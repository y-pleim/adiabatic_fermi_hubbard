from adiabatic_fermi_hubbard import HubbardHamiltonian

from qiskit import QuantumCircuit, execute
from qiskit.result import Result
from qiskit.quantum_info import Pauli
from qiskit_aer import Aer
from qiskit_algorithms import NumPyMinimumEigensolver
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.hamiltonians import FermiHubbardModel
from qiskit_nature.second_q.problems import LatticeModelProblem
from qiskit_nature.second_q.mappers import JordanWignerMapper

import numpy as np


class AdiabaticCircuit:
    """A class for building/executing Qiskit circuits for adiabatic state preparation of a Hubbard Hamiltonian's ground state"""

    def __init__(
        self, ham: HubbardHamiltonian, time_step: float = 0.1, step_count: int = 10
    ):
        self.hubbard_hamiltonian = ham
        self.n = 2 * ham.get_lattice().get_num_sites()
        self.time_step = time_step
        self.step_count = step_count

    def pauli_string_rotation(self, pauli_string: Pauli, argument: float):
        """Performs a rotation about a specified Pauli string through a specified angle.

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
        circ = QuantumCircuit(self.n, 0)  # 0 qubits on classical register
        gates_list = []
        for i in range(len(pauli_string)):
            gates_list.append(pauli_string[i])

        # initial step - converting Y, X gates
        for i in range(len(gates_list)):
            if str(gates_list[i]) == "X":
                circ.h(i)
            elif str(gates_list[i]) == "Y":
                circ.rx(3 * np.pi / 2, i)

        # encode parity onto final qubit
        for i in range(self.n - 1):
            circ.cx(i, self.n - 1)

        # perform rotation
        circ.rz(argument, self.n - 1)

        # uncompute parity
        index = self.n - 2
        while index >= 0:
            circ.cx(index, self.n - 1)
            index -= 1

        # uncompute conversion
        for i in range(len(gates_list)):
            if str(gates_list[i]) == "X":
                circ.h(i)
            elif str(gates_list[i]) == "Y":
                circ.rx(-3 * np.pi / 2, i)

        return circ

    def evolution_operator(self, k):
        """Implements the operation

        .. math:: U(k) = \exp{H_init (1-k/M) \delta t}\exp{H_{FH} (k/M) \delta t}

        where :math:'H_{init} = \sum_i{X_i}', :math:'H_{FH}' is the Jordan-Wigner transformed Fermi-Hubbard Hamiltonian and
        :math:'M' is the number of steps associated with this AdiabaticCircuit object.

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
        """Creates the circuit which performs the adiabatic preparation of the ground state of the HubbardHamiltonian associated
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

        # build a circuit of M+1 steps out of evolution operators
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

        # create FermiHubbardModel with interaction strength t and onsite interaction strength u
        fh_model = FermiHubbardModel(
            lattice.uniform_parameters(
                uniform_interaction=-t,
                uniform_onsite_potential=0.0,  # no potential offset at each site
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
