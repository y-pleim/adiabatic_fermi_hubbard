from qiskit_nature.second_q.hamiltonians.lattices import BoundaryCondition, LineLattice

class Lattice():
    """A class for representing a 1D lattice of fermions."""
    def __init__(self,num_sites:int=2,pbc:bool=True):
        self.num_sites = num_sites
        self.pbc = pbc
        if pbc:
            self.lattice = LineLattice(num_sites,1.0,0.0,BoundaryCondition.PERIODIC)
        else:
            self.lattice = LineLattice(num_sites,1.0,0.0,BoundaryCondition.OPEN)

    def __str__(self):
        return "Number of sites: " + str(self.get_num_sites()) + " sites, \nPeriodic boundary conditions: " + str(self.has_pbc())+ "."
    
    def get_num_sites(self):
        return self.num_sites
    
    def has_pbc(self):
        return self.pbc
    
    def get_qiskit_object(self):
        return self.lattice

    
    
