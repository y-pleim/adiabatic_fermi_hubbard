"""Unit tests for the Lattice class."""

import adiabatic_fermi_hubbard as afh
from qiskit_nature.second_q.hamiltonians.lattices import BoundaryCondition, LineLattice


def test_get_num_sites():
    """Test for Lattice get_num_sites method."""
    lattice1 = afh.Lattice(10, 1)
    assert 10 == lattice1.get_num_sites()


def test_has_pbc():
    """Test for Lattice has_pbc method."""
    lattice1 = afh.Lattice(2, 1)
    lattice2 = afh.Lattice(2, 0)
    assert True == lattice1.has_pbc()
    assert False == lattice2.has_pbc()


def test_get_qiskit_object():
    """Test for Lattice get_qiskit_object method."""
    lattice1 = afh.Lattice(2, 1)
    qiskit_from_lattice1 = lattice1.get_qiskit_object()
    qiskit_lattice = LineLattice(2, 1.0, 0.0, BoundaryCondition.PERIODIC)
    assert type(qiskit_lattice) == type(qiskit_from_lattice1)
    assert qiskit_from_lattice1.boundary_condition == qiskit_lattice.boundary_condition
    assert qiskit_from_lattice1.size == qiskit_lattice.size
    assert qiskit_from_lattice1.edge_parameter == qiskit_lattice.edge_parameter
    assert qiskit_from_lattice1.onsite_parameter == qiskit_lattice.onsite_parameter
