from adiabatic_fermi_hubbard import Lattice
from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.mappers import JordanWignerMapper

class HubbardHamiltonian():
    """A class for constructing the Fermi-Hubbard Hamiltonian for particular lattice, hopping strength, and interaction strength"""
    def __init__(self, lattice:Lattice, t:float=1.0, U:float=1.0):
        self.t = t
        self.U = U
        self.lattice = lattice

        # get single fermionic operators
        n = lattice.get_num_sites()
        i = 0
        creation_operators = []
        annihilation_operators = []

        for i in range(2*n):
            label = str(i)
            creation_operators.append(FermionicOp({"+_"+label : 1.0}, num_spin_orbitals=2*n))
            annihilation_operators.append(FermionicOp({"-_"+label : 1.0}, num_spin_orbitals=2*n))
        
        # build interaction term
        number_operators = []
        for i in range(0,2*n-1,2):
            number_operators.append(creation_operators[i] @ annihilation_operators[i] @ creation_operators[i+1] @ annihilation_operators[i+1])     
        self.interaction_term = U*sum(number_operators)

        # build hopping term
        hopping_operators = []
        for i in range(0,2*n-3,2):
            hopping_operators.append(-1 * creation_operators[i] @ annihilation_operators[i+2] + annihilation_operators[i] @ creation_operators[i+2])

        for i in range(1,2*n-2,2):
            hopping_operators.append(-1 * creation_operators[i] @ annihilation_operators[i+2] + annihilation_operators[i] @ creation_operators[i+2])

        if self.lattice.has_pbc():
            hopping_operators.append(-1 * creation_operators[2*n-2] @ annihilation_operators[0] + annihilation_operators[2*n-2] @ creation_operators[0])
            hopping_operators.append(-1 * creation_operators[2*n-1] @ annihilation_operators[1] + annihilation_operators[2*n-1] @ creation_operators[1])
        
        self.hopping_term = t*sum(hopping_operators)
        self.hamiltonian = self.hopping_term + self.interaction_term

    def __str__(self):
        return "t = " + str(self.t) + ", U = " + str(self.U) + "\nLattice:\n" + str(self.lattice) + "\nAs FermionicOp:\n" + str(self.hamiltonian)
    
    def jw_interaction_term(self):
        jw = JordanWignerMapper()
        jw_interaction = jw.map(self.interaction_term)
        return jw_interaction

    def jw_hopping_term(self):
        jw = JordanWignerMapper()
        jw_hopping = jw.map(self.hopping_term)
        return jw_hopping

    def jw_hamiltonian(self):
        jw = JordanWignerMapper()
        jw_hamiltonian = jw.map(self.hamiltonian)
        return jw_hamiltonian
    
    def get_t_value(self):
        return self.t

    def get_U_value(self):
        return self.U
    
    def get_lattice(self):
        return self.lattice

    def get_interaction_term(self):
        return self.interaction_term

    def get_hopping_term(self):
        return self.hopping_term
    
    def get_hamiltonian(self):
        return self.hamiltonian
