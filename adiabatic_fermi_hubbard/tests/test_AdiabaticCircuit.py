'''Unit tests for theAdiabaticCircuit class.'''
import adiabatic_fermi_hubbard as afh

from qiskit import QuantumCircuit, execute
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit_aer import Aer

import numpy as np

import platform
import copy as cp


def test_str():
    """Test method for str dunder method"""
    lattice1 = afh.Lattice(2, pbc=True)
    ham1 = afh.HubbardHamiltonian(lattice1, 2, 2, 2)
    ad_circ = afh.AdiabaticCircuit(ham1, 0.01, 1000)

    test_string = (
        "Time step: 0.01\nNumber of steps: 1000\nNumber of qubits: 4\n\nHubbardHamiltonian:\n"
        + str(ad_circ.get_hubbard_hamiltonian())
    )
    assert str(ad_circ) == test_string


def test_pauli_string_rotation():
    """Test method for pauli_string_rotation method"""
    # do IXYZ rotation (rightmost position is qubit 0), phase pi
    circ = QuantumCircuit(4, 0)
    circ.h(2)
    circ.rx(3 * np.pi / 2, 1)
    circ.i(3)
    circ.cx(0, 2)
    circ.cx(1, 2)
    circ.rz(np.pi, 2)
    circ.cx(1, 2)
    circ.cx(0, 2)
    circ.rx(-3 * np.pi / 2, 1)
    circ.h(2)

    # use method
    lattice1 = afh.Lattice(2, 1)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1)
    paulis = Pauli("IXYZ")
    circ1 = ad_circ.pauli_string_rotation(paulis, np.pi)

    assert str(circ) == str(circ1)

    # test to cover case where all qubits involved in rotation

    # do XYXY rotation (rightmost position is qubit 0), phase pi
    circ2 = QuantumCircuit(4, 0)
    circ2.rx(3 * np.pi / 2, 0)
    circ2.h(1)
    circ2.rx(3 * np.pi / 2, 2)
    circ2.h(3)
    circ2.cx(0, 3)
    circ2.cx(1, 3)
    circ2.cx(2, 3)
    circ2.rz(np.pi, 3)
    circ2.cx(2, 3)
    circ2.cx(1, 3)
    circ2.cx(0, 3)
    circ2.rx(-3 * np.pi / 2, 0)
    circ2.h(1)
    circ2.rx(-3 * np.pi / 2, 2)
    circ2.h(3)

    paulis2 = Pauli("XYXY")
    circ3 = ad_circ.pauli_string_rotation(paulis2, np.pi)

    assert str(circ2) == str(circ3)


def test_evolution_operator():
    """Test method for evolution_operator"""
    lattice1 = afh.Lattice(2, pbc=False)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1, 0.1, 2)

    circ = ad_circ.evolution_operator(0) # evolution operator for step 0 from method

    # build evolution operator for step 0 by hand
    circ1 = QuantumCircuit(4, 0)

    for i in range(4):
        circ1.rx(0.1, i)

    jw_ham = ham1.jw_hamiltonian()
    paulis = jw_ham.paulis

    seq = [0, 1, 2, 3]

    for pauli in paulis:
        circ1.append(ad_circ.pauli_string_rotation(pauli, 0), seq)

    assert str(circ1.decompose()) == str(circ.decompose())


def test_create_circuit():
    """Test method for create_circuit"""
    lattice1 = afh.Lattice(2, pbc=False)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1, 0.1, 2)

    circ = ad_circ.create_circuit() # circuit from method
    
    # build adiabatic state prep circuit by hand
    circ1 = QuantumCircuit(4, 0)

    for i in range(4):
        circ1.x(i)
        circ1.h(i)

    seq = [0, 1, 2, 3]

    circ1.append(ad_circ.evolution_operator(0), seq)
    circ1.append(ad_circ.evolution_operator(1), seq)
    circ1.append(ad_circ.evolution_operator(2), seq)

    assert str(circ.decompose().decompose()) == str(circ1.decompose().decompose())


def test_run():
    """Test method for run method"""
    lattice1 = afh.Lattice(2, pbc=False)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1, 0.1, 2)

    circ = ad_circ.create_circuit()
    result = ad_circ.run(circ)

    simulator = Aer.get_backend("statevector_simulator")
    result1 = execute(circ, backend=simulator).result()

    assert result.get_statevector().equiv(result1.get_statevector())


def test_run_eigensolver_comparison():
    """Test method for eigensolver_comparison"""
    # no PBC
    lattice1 = afh.Lattice(2, pbc=False)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1, 0.01, 100)

    eigensolver_energy = ad_circ.run_eigensolver_comparison()

    assert np.isclose(eigensolver_energy, -11.403124237432852)

    # PBC
    lattice2 = afh.Lattice(3, pbc=True)
    ham2 = afh.HubbardHamiltonian(lattice2)
    ad_circ2 = afh.AdiabaticCircuit(ham2, 0.01, 100)

    eigensolver_energy2 = ad_circ2.run_eigensolver_comparison()

    assert np.isclose(eigensolver_energy2, -17.150070940649563)



def test_calc_energy():
    """Test method for calc_energy"""
    # no PBC
    lattice1 = afh.Lattice(2, pbc=False)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1, 0.01, 100)

    circuit = ad_circ.create_circuit()
    result = ad_circ.run(circuit)
    energy = ad_circ.calc_energy(result)

    assert np.isclose(energy, -6.156190736486687)

    # PBC
    lattice2 = afh.Lattice(3, pbc=True)
    ham2 = afh.HubbardHamiltonian(lattice2)
    ad_circ2 = afh.AdiabaticCircuit(ham2, 0.01, 100)

    circuit2 = ad_circ2.create_circuit()
    result2 = ad_circ2.run(circuit2)
    energy2 = ad_circ2.calc_energy(result2)

    assert np.isclose(energy2, -9.385097650142908)



def test_get_num_qubits():
    """Test method for get_num_qubits"""
    lattice1 = afh.Lattice(2, pbc=False)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1, 0.1, 2)

    assert ad_circ.get_num_qubits() == 4


def test_get_hubbard_hamiltonian():
    """Test method for get_hubbard_hamiltonian"""
    lattice1 = afh.Lattice(2, pbc=False)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1, 0.1, 2)

    assert str(ham1) == str(ad_circ.get_hubbard_hamiltonian())


def test_get_time_step():
    """Test method for get_time_step"""
    lattice1 = afh.Lattice(2, pbc=False)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1, 0.1, 2)

    assert ad_circ.get_time_step() == 0.1


def test_get_step_count():
    """Test method for get_step_count"""
    lattice1 = afh.Lattice(2, pbc=False)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1, 0.1, 2)

    assert ad_circ.get_step_count() == 2


def test_diagonalize_hamiltonian():
    """Test method for diagonalize_hamiltonian"""
    lattice1 = afh.Lattice(2, pbc=False)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1, 0.1, 2)

    # build matrix form of initial Hamiltonian (all X's)
    a = np.array([(0, 1), (1, 0)])  # X gate

    # create array of n = 2N identity gates
    eyes = []
    for i in range(2 * lattice1.get_num_sites()):
        eyes.append(np.eye(2))

    matrices1 = []
    for i in range(2 * lattice1.get_num_sites()):
        list_2 = cp.copy(eyes)
        list_2[i] = a  # replace the ith identity gate with an X gate

        # take tensor product of all gates
        for j in range(len(list_2)):
            if j == 0:
                b = 1
            b = np.kron(list_2[j], b)
        matrices1.append(b)

    # add together matrices in matrices1 list

    init_matrix = np.zeros(2 ** lattice1.get_num_sites())
    for i in range(len(matrices1)):
        if i == 0:
            init_matrix = matrices1[i]
        else:
            init_matrix = np.add(init_matrix, matrices1[i])

    # manual conversion of JW-transformed Fermi-Hubbard Hamiltonian to matrix form
    ham = ham1.jw_hamiltonian()
    ham_paulis = ham.paulis
    ham_coeffs = ham.coeffs

    matrices = []

    for pauli_string in ham_paulis:
        for i in range(
            len(pauli_string)
        ):  # iterate through all gates in the Pauli string
            if str(pauli_string[i]) == "X":
                a = np.array([(0, 1), (1, 0)])  # numpy implementation of X gate
            elif str(pauli_string[i]) == "Y":
                a = np.asarray(
                    [(0, -1.0j), (1.0j, 0)]
                )  # numpy implementation of Y gate
            elif str(pauli_string[i]) == "Z":
                a = np.array([(1, 0), (0, -1)])  # numpy implementation of Z gate
            else:
                a = np.eye(2) # identity gate

            if i == 0:  # if gate i is the first gate in the Pauli string
                b = 1

            b = np.kron(a, b)  # take tensor product of gate i with previous gates
        matrices.append(b)

    # sum up all matrices
    for i in range(len(matrices)):
        if i == 0:
            final_matrix = ham_coeffs[i] * matrices[i]
        else:
            final_matrix = np.add(final_matrix, ham_coeffs[i] * matrices[i])

    # diagonalize for three cases:

    # 1. all H_{initial}
    k = 0
    matrix = final_matrix * (k / ad_circ.get_step_count()) + init_matrix * (
        1 - k / ad_circ.get_step_count()
    )

    if platform.python_version()[:3] == "3.8":
        energies = np.linalg.eig(matrix)
    else:
        energies = np.linalg.eig(matrix).eigenvalues

    energies_from_method = ad_circ.diagonalize_hamiltonian(0)

    # 2. H_{initial}, H_{final} have equal contributions
    k2 = ad_circ.get_step_count() // 2

    matrix2 = final_matrix * (k2 / ad_circ.get_step_count()) + init_matrix * (
        1 - k2 / ad_circ.get_step_count()
    )

    if platform.python_version()[:3] == "3.8":
        energies_2 = np.linalg.eig(matrix2)
    else:
        energies_2 = np.linalg.eig(matrix2).eigenvalues

    energies_from_method_2 = ad_circ.diagonalize_hamiltonian(k2)

    # 3. all H_{final}
    k3 = ad_circ.get_step_count()

    matrix3 = final_matrix * (k3 / ad_circ.get_step_count()) + init_matrix * (
        1 - k3 / ad_circ.get_step_count()
    )

    if platform.python_version()[:3] == "3.8":
        energies_3 = np.linalg.eig(matrix3)
    else:
        energies_3 = np.linalg.eig(matrix3).eigenvalues

    energies_from_method_3 = ad_circ.diagonalize_hamiltonian(k3)

    assert str(energies) == str(energies_from_method)
    assert str(energies_2) == str(energies_from_method_2)
    assert str(energies_3) == str(energies_from_method_3)


def test_initialize_ising():
    """Test method for initialize_ising"""
    lattice1 = afh.Lattice(2, pbc=False)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1, 0.1, 2)

    ad_circ.initialize_ising(1)

    assert ad_circ.ising_ham == SparsePauliOp(
        ["ZZII", "IZZI", "IIZZ"], coeffs=[1 + 0.0j, 1 + 0.0j, 1 + 0.0j]
    )
    assert ad_circ.ising_n == ad_circ.get_num_qubits()

    # PBC case
    lattice2 = afh.Lattice(4, True)
    ham2 = afh.HubbardHamiltonian(lattice2)
    ad_circ1 = afh.AdiabaticCircuit(ham2, 0.1, 2)

    ad_circ1.initialize_ising(1)

    assert ad_circ1.ising_ham == SparsePauliOp(
        ["ZZII", "IZZI", "IIZZ", "ZIIZ"],
        coeffs=[1 + 0.0j, 1 + 0.0j, 1 + 0.0j, 1 + 0.0j],
    )


def test_ising_evolution_operator():
    """Test method for ising_evolution_operator"""
    lattice1 = afh.Lattice(2, pbc=False)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1, 0.1, 2)

    ad_circ.initialize_ising(1)

    circ = ad_circ.ising_evolution_operator(0)
    circ1 = QuantumCircuit(4, 0)

    for i in range(4):
        circ1.rx(0.1, i)

    paulis = ad_circ.ising_ham.paulis

    seq = [0, 1, 2, 3]

    for pauli in paulis:
        circ1.append(ad_circ.pauli_string_rotation(pauli, 0), seq)

    assert str(circ1.decompose()) == str(circ.decompose())


def test_ising_create_circuit():
    """Test method for ising_create_circuit"""
    lattice1 = afh.Lattice(2, pbc=False)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1, 0.1, 2)

    ad_circ.initialize_ising(1)

    circ = ad_circ.ising_create_circuit()
    circ1 = QuantumCircuit(4, 0)

    for i in range(4):
        circ1.x(i)
        circ1.h(i)

    seq = [0, 1, 2, 3]

    circ1.append(ad_circ.ising_evolution_operator(0), seq)
    circ1.append(ad_circ.ising_evolution_operator(1), seq)
    circ1.append(ad_circ.ising_evolution_operator(2), seq)

    assert str(circ.decompose().decompose()) == str(circ1.decompose().decompose())


def test_ising_calc_energy():
    """Test method for ising_calc_energy"""

    # PBC
    lattice1 = afh.Lattice(2, pbc=True)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1, 0.01, 100)

    ad_circ.initialize_ising(1)

    circuit = ad_circ.ising_create_circuit()
    result = ad_circ.run(circuit)
    ising_energy = ad_circ.ising_calc_energy(result)

    assert np.isclose(ising_energy, -0.30494050361526925)

    # no PBC
    lattice2 = afh.Lattice(2, pbc=False)
    ham2 = afh.HubbardHamiltonian(lattice2)
    ad_circ2 = afh.AdiabaticCircuit(ham2, 0.01, 100)

    ad_circ2.initialize_ising(1)

    circuit2 = ad_circ2.ising_create_circuit()
    result2 = ad_circ2.run(circuit2)
    ising_energy2 = ad_circ2.ising_calc_energy(result2)

    assert np.isclose(ising_energy2, -0.23057848327680522)
