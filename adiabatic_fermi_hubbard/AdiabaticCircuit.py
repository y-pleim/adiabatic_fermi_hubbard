from adiabatic_fermi_hubbard import Lattice, HubbardHamiltonian
import qiskit_algorithms as qa
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.hamiltonians import FermiHubbardModel
from qiskit_nature.second_q.hamiltonians.lattices import LineLattice,BoundaryCondition
from qiskit_nature.second_q.problems import LatticeModelProblem
import numpy as np
import qiskit as qi


class AdiabaticCircuit():
    """A class for building/executing Qiskit circuits for adiabatic state preparation of a Hubbard Hamiltonian's ground state"""
    def __init__(self,ham:HubbardHamiltonian,time_step:float=0.1,step_count:int=10):
        self.hamiltonian = ham
        self.n = 2 * ham.get_lattice().get_num_sites()
        self.time_step = time_step
        self.step_count = step_count

    def pauli_string_rotation(self,pauli_string,argument):
        circ = qi.QuantumCircuit(self.n,0) # 0 qubits on classical register
        gates_list = []
        for i in range(len(pauli_string)):
            gates_list.append(pauli_string[i])
        
        #initial step - convert Y, X gates
        print(gates_list)
        print(gates_list[0])
        for i in range(len(gates_list)):
            if str(gates_list[i]) == "X":
                circ.h(i)
            elif str(gates_list[i]) == "Y":
                circ.rx(np.pi/2, i)
    
        for i in range(self.n-1):
            circ.cx(i,self.n-1)
        
        circ.rz(argument,self.n-1)

        index = self.n-2
        while index >= 0:
            circ.cx(index,self.n-1)
            index -= 1
        
        for i in range(len(gates_list)):
            if str(gates_list[i]) == "X":
                circ.h(i)
            elif str(gates_list[i]) == "Y":
                circ.rx(-np.pi/2, i)
        
        return circ
                    

    def get_circuit(self):
        pass

    def get_num_qubits(self):
        return self.n