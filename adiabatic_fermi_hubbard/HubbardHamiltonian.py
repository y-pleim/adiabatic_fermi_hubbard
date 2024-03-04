class HubbardHamiltonian():
    """A class for constructing the Fermi-Hubbard Hamiltonian for particular lattice, hopping strength, and interaction strength"""
    def __init__(self, lattice,t:float=1.0,U:float=1.0):
        self.t = t
        self.U = U
        self.lattice = lattice
