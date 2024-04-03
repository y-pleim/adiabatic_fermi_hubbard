"""Unit tests for HubbardHamiltonian. """

import adiabatic_fermi_hubbard as afh
from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.hamiltonians import FermiHubbardModel
from qiskit.quantum_info import SparsePauliOp


def test_str():
    """Test for __str__ method"""
    lattice1 = afh.Lattice(2, 1)
    ham1 = afh.HubbardHamiltonian(lattice1, 2, 2, 2)
    test_string = (
        "t = 2\nU = 2\nmu = 2\n\nLattice:\n"
        + str(lattice1)
        + "\n\n"
        + str(ham1.get_hamiltonian())
    )
    assert str(ham1) == test_string


def test_jw_interaction_term():
    """Test for jw_interaction_term method"""
    lattice1 = afh.Lattice(2, 1)
    ham1 = afh.HubbardHamiltonian(lattice1, 2, 2, 2)
    test_pauli_op = SparsePauliOp(
        ["IIII", "IIZI", "IIIZ", "IIZZ", "ZIII", "IZII", "ZZII"],
        coeffs=[
            1.0 + 0.0j,
            -0.5 + 0.0j,
            -0.5 + 0.0j,
            0.5 + 0.0j,
            -0.5 + 0.0j,
            -0.5 + 0.0j,
            0.5 + 0.0j,
        ],
    )
    assert test_pauli_op == ham1.jw_interaction_term()


def test_jw_hopping_term():
    """Test for jw_hopping_term method"""
    lattice1 = afh.Lattice(2, 1)
    ham1 = afh.HubbardHamiltonian(lattice1, 2, 2, 2)
    test_pauli_op = SparsePauliOp(
        ["IYZY", "IXZX", "YZYI", "XZXI"],
        coeffs=[-1.0 + 0.0j, -1.0 + 0.0j, -1.0 + 0.0j, -1.0 + 0.0j],
    )
    assert test_pauli_op == ham1.jw_hopping_term()


def test_jw_chemical_potential_term():
    """Test for jw_chemical_potential_term method"""
    lattice1 = afh.Lattice(2, 1)
    ham1 = afh.HubbardHamiltonian(lattice1, 2, 2, 2)
    test_pauli_op = SparsePauliOp(
        ["IIII", "IIIZ", "IIZI", "IZII", "ZIII"],
        coeffs=[4 + 0.0j, -1 + 0.0j, -1 + 0.0j, -1 + 0.0j, -1 + 0.0j],
    )
    assert test_pauli_op == ham1.jw_chemical_potential_term()


def test_jw_hamiltonian():
    """Test for jw_hamiltonian method"""
    lattice1 = afh.Lattice(2, 1)
    ham1 = afh.HubbardHamiltonian(lattice1, 2, 2, 2)
    test_pauli_op = SparsePauliOp(
        [
            "IYZY",
            "IXZX",
            "YZYI",
            "XZXI",
            "IIII",
            "ZIII",
            "IZII",
            "ZZII",
            "IIZI",
            "IIIZ",
            "IIZZ",
        ],
        coeffs=[
            -1.0 + 0.0j,
            -1.0 + 0.0j,
            -1.0 + 0.0j,
            -1.0 + 0.0j,
            5.0 + 0.0j,
            -1.5 + 0.0j,
            -1.5 + 0.0j,
            0.5 + 0.0j,
            -1.5 + 0.0j,
            -1.5 + 0.0j,
            0.5 + 0.0j,
        ],
    )

    test_pauli_op2 = SparsePauliOp(
        [
            "IYZY",
            "IXZX",
            "YZYI",
            "XZXI",
            "IIII",
            "IIZI",
            "IIIZ",
            "IIZZ",
            "ZIII",
            "IZII",
            "ZZII",
        ],
        coeffs=[
            -1.0 + 0.0j,
            -1.0 + 0.0j,
            -1.0 + 0.0j,
            -1.0 + 0.0j,
            5.0 + 0.0j,
            -1.5 + 0.0j,
            -1.5 + 0.0j,
            0.5 + 0.0j,
            -1.5 + 0.0j,
            -1.5 + 0.0j,
            0.5 + 0.0j,
        ],
    )

    assert (
        test_pauli_op == ham1.jw_hamiltonian()
        or test_pauli_op2 == ham1.jw_hamiltonian()
    )


def test_get_t_value():
    """Test for get_t_value() access method."""
    lattice1 = afh.Lattice(4, 1)
    ham1 = afh.HubbardHamiltonian(lattice1, 2, 2, 2)
    assert ham1.get_t_value() == 2.0


def test_get_U_value():
    """Test for get_U_value() access method."""
    lattice1 = afh.Lattice(4, 1)
    ham1 = afh.HubbardHamiltonian(lattice1, 2, 2, 2)
    assert ham1.get_U_value() == 2.0


def test_get_mu_value():
    """Test for get_mu_value() access method."""
    lattice1 = afh.Lattice(4, 1)
    ham1 = afh.HubbardHamiltonian(lattice1, 2, 2, 2)
    assert ham1.get_mu_value() == 2.0


def test_get_lattice():
    """Test for get_lattice() access method."""
    lattice1 = afh.Lattice(4, 1)
    ham1 = afh.HubbardHamiltonian(lattice1, 2, 2, 2)
    lattice_from_ham = ham1.get_lattice()

    assert type(lattice1) == type(lattice_from_ham)
    assert lattice1.get_num_sites() == lattice_from_ham.get_num_sites()
    assert lattice1.has_pbc() == lattice_from_ham.has_pbc()
    assert type(lattice1.get_qiskit_object()) == type(
        lattice_from_ham.get_qiskit_object()
    )


def test_get_interaction_term():
    """Test for get_interaction_term() method."""
    lattice1 = afh.Lattice(2, 1)
    ham1 = afh.HubbardHamiltonian(lattice1, 2, 2, 2)
    interaction_ham1 = ham1.get_interaction_term()

    create_0 = FermionicOp({"+_0": 1.0}, num_spin_orbitals=4)
    create_1 = FermionicOp({"+_1": 1.0}, num_spin_orbitals=4)
    create_2 = FermionicOp({"+_2": 1.0}, num_spin_orbitals=4)
    create_3 = FermionicOp({"+_3": 1.0}, num_spin_orbitals=4)

    anni_0 = FermionicOp({"-_0": 1.0}, num_spin_orbitals=4)
    anni_1 = FermionicOp({"-_1": 1.0}, num_spin_orbitals=4)
    anni_2 = FermionicOp({"-_2": 1.0}, num_spin_orbitals=4)
    anni_3 = FermionicOp({"-_3": 1.0}, num_spin_orbitals=4)

    interaction_test = 2.0 * (
        create_0 @ anni_0 @ create_1 @ anni_1 + create_2 @ anni_2 @ create_3 @ anni_3
    )

    assert interaction_ham1 == interaction_test


def test_get_hopping_term():
    """Test for get_hopping_term() method."""
    lattice1 = afh.Lattice(2, 1)
    ham1 = afh.HubbardHamiltonian(lattice1, 2, 2, 2)
    hopping_ham1 = ham1.get_hopping_term()

    # open boundary conditions case
    create_0 = FermionicOp({"+_0": 1.0}, num_spin_orbitals=4)
    create_1 = FermionicOp({"+_1": 1.0}, num_spin_orbitals=4)
    create_2 = FermionicOp({"+_2": 1.0}, num_spin_orbitals=4)
    create_3 = FermionicOp({"+_3": 1.0}, num_spin_orbitals=4)

    anni_0 = FermionicOp({"-_0": 1.0}, num_spin_orbitals=4)
    anni_1 = FermionicOp({"-_1": 1.0}, num_spin_orbitals=4)
    anni_2 = FermionicOp({"-_2": 1.0}, num_spin_orbitals=4)
    anni_3 = FermionicOp({"-_3": 1.0}, num_spin_orbitals=4)

    hopping_test = 2.0 * (
        (-1 * (create_0 @ anni_2) + anni_0 @ create_2)
        + (-1 * (create_1 @ anni_3) + anni_1 @ create_3)
    )

    # periodic boundary conditions case
    lattice2 = afh.Lattice(3, 1)
    ham2 = afh.HubbardHamiltonian(lattice2, 2, 2)
    hopping_ham2 = ham2.get_hopping_term()

    assert hopping_ham1 == hopping_test

    create_0 = FermionicOp({"+_0": 1.0}, num_spin_orbitals=6)
    create_1 = FermionicOp({"+_1": 1.0}, num_spin_orbitals=6)
    create_2 = FermionicOp({"+_2": 1.0}, num_spin_orbitals=6)
    create_3 = FermionicOp({"+_3": 1.0}, num_spin_orbitals=6)
    create_4 = FermionicOp({"+_4": 1.0}, num_spin_orbitals=6)
    create_5 = FermionicOp({"+_5": 1.0}, num_spin_orbitals=6)

    anni_0 = FermionicOp({"-_0": 1.0}, num_spin_orbitals=6)
    anni_1 = FermionicOp({"-_1": 1.0}, num_spin_orbitals=6)
    anni_2 = FermionicOp({"-_2": 1.0}, num_spin_orbitals=6)
    anni_3 = FermionicOp({"-_3": 1.0}, num_spin_orbitals=6)
    anni_4 = FermionicOp({"-_4": 1.0}, num_spin_orbitals=6)
    anni_5 = FermionicOp({"-_5": 1.0}, num_spin_orbitals=6)

    hopping_test_2 = 2.0 * (
        (-1 * (create_0 @ anni_2) + anni_0 @ create_2)
        + (-1 * (create_2 @ anni_4) + anni_2 @ create_4)
        + (-1 * (create_1 @ anni_3) + anni_1 @ create_3)
        + (-1 * (create_3 @ anni_5) + anni_3 @ create_5)
        + (-1 * (create_4 @ anni_0) + anni_4 @ create_0)
        + (-1 * (create_5 @ anni_1) + anni_5 @ create_1)
    )

    assert hopping_ham2 == hopping_test_2


def test_get_chemical_potential_term():
    """Test for get_chemical_potential() method."""
    lattice1 = afh.Lattice(2, 1)
    ham1 = afh.HubbardHamiltonian(lattice1, 2, 2, 2)
    chem_ham1 = ham1.get_chemical_potential_term()

    create_0 = FermionicOp({"+_0": 1.0}, num_spin_orbitals=4)
    create_1 = FermionicOp({"+_1": 1.0}, num_spin_orbitals=4)
    create_2 = FermionicOp({"+_2": 1.0}, num_spin_orbitals=4)
    create_3 = FermionicOp({"+_3": 1.0}, num_spin_orbitals=4)

    anni_0 = FermionicOp({"-_0": 1.0}, num_spin_orbitals=4)
    anni_1 = FermionicOp({"-_1": 1.0}, num_spin_orbitals=4)
    anni_2 = FermionicOp({"-_2": 1.0}, num_spin_orbitals=4)
    anni_3 = FermionicOp({"-_3": 1.0}, num_spin_orbitals=4)

    chem_test = 2.0 * (
        create_0 @ anni_0 + create_1 @ anni_1 + create_2 @ anni_2 + create_3 @ anni_3
    )

    assert chem_ham1 == chem_test


def test_get_hamiltonian():
    """Test for get_hamiltonian() method."""
    lattice1 = afh.Lattice(2, 1)
    ham1 = afh.HubbardHamiltonian(lattice1, 2, 2, 2)
    assert (
        ham1.get_hopping_term()
        + ham1.get_interaction_term()
        + ham1.get_chemical_potential_term()
        == ham1.get_hamiltonian()
    )


def test_get_hamiltonian_2():
    """Additional test for get_hamiltonian() based on the qiskit-nature FermiHubbardModel object"""
    lattice1 = afh.Lattice(2, 1)
    ham1 = afh.HubbardHamiltonian(lattice1, 2, 2, 2)
    qiskit_lattice = lattice1.get_qiskit_object()

    fermi_hubbard_test = FermiHubbardModel(
        qiskit_lattice.uniform_parameters(
            uniform_interaction=-2.0, uniform_onsite_potential=2.0
        ),
        onsite_interaction=2.0,
    )

    assert ham1.get_hamiltonian() == fermi_hubbard_test.second_q_op()
