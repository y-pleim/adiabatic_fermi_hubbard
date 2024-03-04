class AdiabaticCircuit():
    """A class for building/executing Qiskit circuits for adiabatic state preparation of a Hubbard Hamiltonian's ground state"""
    def __init__(self,ham,time_step:float=0.1,step_count:int=10):
        self.ham = ham
        self.time_step = time_step
        self.step_count = step_count
