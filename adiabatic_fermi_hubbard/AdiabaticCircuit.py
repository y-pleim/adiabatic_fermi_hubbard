from adiabatic_fermi_hubbard import HubbardHamiltonian

from qiskit import QuantumCircuit, execute, QuantumRegister
from qiskit.result import Result
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit_aer import Aer
from qiskit_algorithms import NumPyMinimumEigensolver
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.hamiltonians import FermiHubbardModel
from qiskit_nature.second_q.problems import LatticeModelProblem
from qiskit_nature.second_q.mappers import JordanWignerMapper

import numpy as np
import platform

class AdiabaticCircuit:
    def __init__(
        self, ham: HubbardHamiltonian, time_step: float = 0.01, step_count: int = 1000
    ):
        """A class for building/executing Qiskit circuits for adiabatic state preparation of the Hubbard Hamiltonian ground state through
        interpolating between :math:'H_{init}' and :math:'H_{Fermi-Hubbard}' according to

        .. math:: H(k) = H_{init} (1-k/M) + H_{final} (k/M).

        where :math:'M' is the number of interpolation steps, :math:'k = 0, 1, ... , M', and :math:'H_{init} = \\sum_i X_i'

        Parameters
        ----------
        ham : HubbardHamiltonian
            The HubbardHamiltonian object to find the ground state for.
        time_step : float
            Duration of each interpolating step.
        step_count : float
            The number of steps to take when interpolating between :math:'H_{init}' and :math:'H_{Fermi-Hubbard}'.
        """

        self.hubbard_hamiltonian = ham
        self.n = 2 * ham.get_lattice().get_num_sites()
        self.time_step = time_step
        self.step_count = step_count

    def __str__(self):
        return (
            "Time step: "
            + str(self.time_step)
            + "\nNumber of steps: "
            + str(self.step_count)
            + "\nNumber of qubits: "
            + str(self.n)
            + "\n\nHubbardHamiltonian:\n"
            + str(self.hubbard_hamiltonian)
        )

    def pauli_string_rotation(self, pauli_string: Pauli, argument: float):
        """Performs a rotation about a specified Pauli string through a specified angle. Based on Nielsen and Chuang, Ch 4.

        Parameters
        ----------
        pauli_string : Pauli
            A qiskit Pauli object which contains an n-qubit Pauli string.
        argument : float
            The angle through which to rotate.

        Returns
        -------
        circ : QuantumCircuit
            A qiskit QuantumCircuit object which performs the rotation when executed.

        """
        circ = QuantumCircuit(self.n, 0)

        # extract constituent Pauli gates
        gates_list = []
        for i in range(len(pauli_string)):
            gates_list.append(pauli_string[i])

        # list of qubits which participate in the rotation
        rotation_qubits = []

        # initial step - converting Y, X gates and assigning rotation qubits
        for i in range(len(gates_list)):
            if str(gates_list[i]) == "X":
                circ.h(i)
                rotation_qubits.append(i)
            elif str(gates_list[i]) == "Y":
                circ.rx(3 * np.pi / 2, i)
                rotation_qubits.append(i)
            elif str(gates_list[i]) == "Z":
                rotation_qubits.append(i)
            else:
                circ.id(i)

        # encode parity onto last qubit in rotation_qubits
        j = 0
        for i in rotation_qubits:
            if len(rotation_qubits) == len(
                gates_list
            ):  # if all qubits participate in the rotation
                j = len(pauli_string) - 1
                if i < j:
                    circ.cx(i, j)  # encode parity onto last qubit
            elif i != max(
                rotation_qubits
            ):  # if qubit i isn't the last one that participates in the rotation
                j = max(rotation_qubits)
                circ.cx(i, j)  # encode parity onto last qubit in rotation_qubits

        # perform rotation on last qubit
        if rotation_qubits:
            circ.rz(argument, j)

        # undoes parity encoding
        index = self.n - 2
        while index >= 0:
            if index in rotation_qubits and len(rotation_qubits) == len(gates_list):
                j = len(pauli_string) - 1
                circ.cx(index, j)
            elif index in rotation_qubits and index != max(rotation_qubits):
                j = max(rotation_qubits)
                circ.cx(index, j)
            index -= 1

        # undoes conversion
        for i in range(len(gates_list)):
            if str(gates_list[i]) == "X":
                circ.h(i)
            elif str(gates_list[i]) == "Y":
                circ.rx(-3 * np.pi / 2, i)

        return circ

    def evolution_operator(self, k):
        """Implements the operation

        .. math:: U(k) = \\exp{H_{init} (1-k/M) \\delta t}\\exp{H_{FH} (k/M) \\delta t}

        where :math:'H_{init} = \\sum_i{X_i}', :math:'H_{FH}' is the Jordan-Wigner transformed Fermi-Hubbard Hamiltonian and
        :math:'M' is the number of interpolation steps.

        Parameters
        ----------
        k : integer
            The index k as defined above.

        Returns
        -------
        circ : QuantumCircuit
            A qiskit QuantumCircuit object which implements the evolution.

        """
        delta_s = 1 / self.step_count

        circ = QuantumCircuit(self.n, 0)

        final_ham = (
            self.hubbard_hamiltonian.jw_hamiltonian()
        )  # get SparsePauliOp representation of HubbardHamiltonian
        final_ham_paulis = final_ham.paulis  # paulis list
        final_ham_coeffs = final_ham.coeffs  # weights list

        # stores qubit indices
        seq = []
        for i in range(self.n):
            seq.append(i)

        # first exponential factor: evolution under H_init for step k on qubits in seq
        for i in range(self.n):
            circ.rx(self.time_step * (1 - k * delta_s), i)

        # second exponential factor: evolution under H_HF for step k on qubits in seq
        for i in range(len(final_ham_paulis)):
            circ.append(
                self.pauli_string_rotation(
                    final_ham_paulis[i],
                    self.time_step * k * delta_s * np.real(final_ham_coeffs[i]),
                ),
                seq,
            )

        return circ

    def create_circuit(self):
        """Creates the circuit which performs the adiabatic state preparation of the ground state of the HubbardHamiltonian associated
        with this AdiabaticCircuit.

        Returns
        -------
        circ : QuantumCircuit
            The circuit which will carry out the adiabatic state preparation when run.

        """
        circ = QuantumCircuit(self.n, 0)

        # start system in ground state of H_init (i.e., the all minus state)
        for i in range(self.n):
            circ.x(i)
            circ.h(i)

        # store qubit indices
        seq = []
        for i in range(self.n):
            seq.append(i)

        # build a circuit of M+1 evolution operators, starting from k = 0 and ending at k = M
        for i in range(self.step_count + 1):
            circ.append(self.evolution_operator(i), seq)

        return circ

    def run(self, circ: QuantumCircuit):
        """Executes a quantum circuit using Qiskit's statevector simulator.

        Parameters
        ----------
        circ : QuantumCircuit
            The circuit to be run.

        Returns
        -------
        result : Result
            A qiskit Result object containing the execution results.
        """
        simulator = Aer.get_backend("statevector_simulator")
        result = execute(circ, backend=simulator).result()
        return result

    def run_eigensolver_comparison(self):
        """Computes the ground state energy for the HubbardHamiltonian associated with the AdiabaticCircuit
        through qiskit-nature methods. Based on the example code found in 'LatticeModels'_.

        .. _LatticeModels: https://qiskit-community.github.io/qiskit-nature/tutorials/10_lattice_models.html#The-Fermi-Hubbard-model

        Returns
        -------
        result : float
            The ground state energy.
        """
        map = JordanWignerMapper()

        ham = self.hubbard_hamiltonian
        lattice = ham.get_lattice().get_qiskit_object()
        t = ham.get_t_value()
        u = ham.get_U_value()
        mu = ham.get_mu_value()

        # create FermiHubbardModel with interaction strength t, onsite interaction strength u, chemical potential mu
        fh_model = FermiHubbardModel(
            lattice.uniform_parameters(
                uniform_interaction=-t,
                uniform_onsite_potential=mu,
            ),
            onsite_interaction=u,
        )

        # setting up the problem and solver
        problem = LatticeModelProblem(fh_model)
        solver = NumPyMinimumEigensolver(
            filter_criterion=problem.get_default_filter_criterion()
        )
        gs_solver = GroundStateEigensolver(map, solver)

        result = gs_solver.solve(problem).eigenvalues[
            0
        ]  # extracts ground state energy as float

        return result

    def calc_energy(self, result: Result):
        """Calculates the energy for the state represented by result using the HubbardHamiltonian associated
        with this AdiabaticCircuit object.

        Parameters
        ----------
        result : Result
            A qiskit Result containing the state to calculate the energy for.

        Returns
        -------
        energy : float
            The energy of the state.

        """
        result_statevector = result.get_statevector()
        ham = self.hubbard_hamiltonian
        jw_ham = ham.jw_hamiltonian()
        energy = np.real(result_statevector.expectation_value(jw_ham))

        return energy

    def get_num_qubits(self):
        """Access method to get the number of qubits associated with this AdiabaticCircuit object.

        Returns
        -------
        n : int
            The number of qubits.
        """
        return self.n

    def get_hubbard_hamiltonian(self):
        """Access method to get the HubbardHamiltonian associated with this AdiabaticCircuit object.

        Returns
        -------
        hubbard_hamiltonian : HubbardHamiltonian
            The HubbardHamiltonian.
        """
        return self.hubbard_hamiltonian

    def get_time_step(self):
        """Access method to get the time step associated with this AdiabaticCircuit object.

        Returns
        -------
        time_step : float
            The time step involved in creating the adiabatic state preparation circuit.
        """
        return self.time_step

    def get_step_count(self):
        """Access method to get the step count associated with this AdiabaticCircuit object.

        Returns
        -------
        step_count : int
            The number of steps involved in creating the adiabatic state preparation circuit.
        """
        return self.step_count

    def diagonalize_ham(self):
        """Diagonalize the final Fermi-Hubbard Hamiltonian to validate adiabatic state result and qiskit-nature eigensolver result."""
        ham = self.hubbard_hamiltonian.jw_hamiltonian()
        ham_paulis = ham.paulis
        ham_coeffs = ham.coeffs

        matrices = []

        for pauli_string in ham_paulis:
            for i in range(
                len(pauli_string)
            ):  # iterate through all gates in the Pauli string
                if str(pauli_string[i]) == "X":
                    a = np.asarray([(0, 1), (1, 0)])  # numpy implementation of X gate
                elif str(pauli_string[i]) == "Y":
                    a = np.asarray(
                        [(0, -1.0j), (1.0j, 0)]
                    )  # numpy implementation of Y gate
                elif str(pauli_string[i]) == "Z":
                    a = np.array([(1, 0), (0, -1)])  # numpy implementation of Z gate
                else:
                    a = np.eye(2)

                if i == 0:  # if gate i is the first gate in the Pauli string
                    b = 1

                b = np.kron(a, b)  # take tensor product of gate i with previous gates
            matrices.append(b)

        # sum up all matrices
        for i in range(len(matrices)):
            if i == 0:
                total_matrix = ham_coeffs[i] * matrices[i]
            else:
                total_matrix = np.add(total_matrix, ham_coeffs[i] * matrices[i])

        # find array of eigenvalues

        if platform.python_version()[:3] == '3.8':
            eigs = np.linalg.eig(total_matrix)
        else:
            eigs = np.linalg.eig(total_matrix).eigenvalues

        return eigs

    def ising_setup(self, exchange):
        """Support method that assigns an Ising Hamiltonian (for the same number of qubits) to the AdiabaticCircuit object:

        .. math:: H = J\\sum_{<i,j>} Z_i Z_j

        where :math:'<i,j>' denotes nearest neighbors on the lattice. Implemented with periodic boundary
        conditions if the Lattice associated with the AdiabaticCircuit has periodic boundary conditions.

        Parameters
        ----------
        exchange : float
            value of J parameter

        """
        self.ising_n = self.n
        val = exchange
        if self.get_hubbard_hamiltonian().get_lattice().has_pbc():
            self.ising_ham = SparsePauliOp(
                ["ZZII", "IZZI", "IIZZ", "ZIIZ"], coeffs=[val, val, val, val]
            )
        else:
            self.ising_ham = SparsePauliOp(
                ["ZZII", "IZZI", "IIZZ"], coeffs=[val, val, val]
            )

    def ising_evolution_operator(self, k):
        """Support method that implements the evolution operator

        .. math:: U(k) = \\exp{H_{init} (1-k/M) \\delta t}\\exp{H_{Ising} (k/M) \\delta t}

        where :math:'H_{Ising}' is the Ising Hamiltonian created by running the ising_setup method, :math:'M' is the step count,
        and :math:'k = 0, 1, ..., M".

        Parameters
        ----------
        k : integer
            The index k as defined above.

        Returns
        -------
        circ : QuantumCircuit
            A qiskit QuantumCircuit object which implements the evolution.

        """
        delta_s = 1 / self.step_count

        circ = QuantumCircuit(self.ising_n, 0)

        final_ham = self.ising_ham
        final_ham_paulis = final_ham.paulis  # paulis list
        final_ham_coeffs = final_ham.coeffs  # weights list

        # stores qubit indices
        seq = []
        for i in range(self.ising_n):
            seq.append(i)

        # first exponential factor: evolution under H_init for step k on qubits in seq
        for i in range(self.ising_n):
            circ.rx(self.time_step * (1 - k * delta_s), i)

        # second exponential factor: evolution under H_ising for step k on qubits in seq
        for i in range(len(final_ham_paulis)):
            circ.append(
                self.pauli_string_rotation(
                    final_ham_paulis[i],
                    self.time_step * k * delta_s * np.real(final_ham_coeffs[i]),
                ),
                seq,
            )
        return circ

    def ising_create_circuit(self):
        """Creates the circuit which performs the adiabatic state preparation of the ground state of the Ising Hamiltonian associated
        with this AdiabaticCircuit (provided the ising_setup method has already been run).

        Returns
        -------
        circ : QuantumCircuit
            The circuit which will carry out the adiabatic state preparation when run.

        """
        circ = QuantumCircuit(self.ising_n, 0)

        # start system in ground state of H_init (i.e., the all minus state)
        for i in range(self.ising_n):
            circ.x(i)
            circ.h(i)

        # store qubit indices
        seq = []
        for i in range(self.ising_n):
            seq.append(i)

        # build a circuit of M+1 evolution operators with i = 0, 1, ..., M
        for i in range(self.step_count + 1):
            circ.append(self.ising_evolution_operator(i), seq)

        return circ

    def ising_calc_energy(self, result: Result):
        """Calculates the energy for the state represented by result using the Ising Hamiltonian associated
        with this AdiabaticCircuit object (provided the ising_setup method has already been run).

        Parameters
        ----------
        result : Result
            A qiskit Result containing the state to calculate the energy for.

        Returns
        -------
        energy : float
            The energy of the state.

        """
        result_statevector = result.get_statevector()
        ham = self.ising_ham
        print(ham)
        energy = np.real(result_statevector.expectation_value(ham))

        return energy
