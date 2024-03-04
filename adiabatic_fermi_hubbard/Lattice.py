class Lattice():
    """A class for representing a lattice of fermions."""
    def __init__(self,type:str="1D",num_sites:int=2,pbc:bool=True):
        self.num_sites = num_sites