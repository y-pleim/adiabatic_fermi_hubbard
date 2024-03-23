import adiabatic_fermi_hubbard as afh
import qiskit as qi
import numpy as np

def test_pauli_string_rotation():
    '''Test method for pauli_string_rotation function'''
    #do XYIZ rotation, phase pi
    circ = qi.QuantumCircuit(4,0)
    circ.h(0)
    circ.rx(3*np.pi/2,1)
    circ.cx(0,3)
    circ.cx(1,3)
    circ.cx(2,3)
    circ.rz(np.pi, 3)
    circ.cx(2,3)
    circ.cx(1,3)
    circ.cx(0,3)
    circ.rx(-3*np.pi/2,1)
    circ.h(0)

    circ_drawn = circ.draw(output="text")

    #use method
    lattice1 = afh.Lattice(2,1)
    ham1 = afh.HubbardHamiltonian(lattice1)
    ad_circ = afh.AdiabaticCircuit(ham1)
    paulis = qi.quantum_info.Pauli("XYIZ")
    circ1 = ad_circ.pauli_string_rotation(paulis,np.pi)
    circ1_drawn = circ1.draw(output="text")

    assert str(circ) == str(circ1)