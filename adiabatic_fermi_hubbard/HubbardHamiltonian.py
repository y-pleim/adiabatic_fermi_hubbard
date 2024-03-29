from adiabatic_fermi_hubbard import Lattice
from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.mappers import JordanWignerMapper


class HubbardHamiltonian:
    """A class for constructing the Fermi-Hubbard Hamiltonian for particular lattice, hopping strength, and interaction strength.

    Parameters
    ----------
    lattice : Lattice
        The Lattice to find the Hamiltonian for.
    t : float
        The strength of the hopping term.
    U : float
        The strength of the onsite interaction term.
    """

    def __init__(self, lattice: Lattice, t: float = 1.0, U: float = 1.0):
        self.t = t
        self.U = U
        self.lattice = lattice

        # get single fermionic operators
        n = lattice.get_num_sites()  # number of lattice sites
        i = 0
        creation_operators = []
        annihilation_operators = []

        for i in range(
            2 * n
        ):  # make all FermionicOp creation/annihilation operators for 2*n spin orbitals
            label = str(i)
            creation_operators.append(
                FermionicOp({"+_" + label: 1.0}, num_spin_orbitals=2 * n)
            )
            annihilation_operators.append(
                FermionicOp({"-_" + label: 1.0}, num_spin_orbitals=2 * n)
            )

        # build interaction term
        number_operators = []
        for i in range(0, 2 * n - 1, 2):
            number_operators.append(
                creation_operators[i]
                @ annihilation_operators[i]
                @ creation_operators[i + 1]
                @ annihilation_operators[i + 1]
            )
        self.interaction_term = U * sum(number_operators)

        # build hopping term
        hopping_operators = []

        # hopping factors for up-spin fermions (operators with even indices)
        for i in range(0, 2 * n - 3, 2):
            hopping_operators.append(
                -1 * creation_operators[i] @ annihilation_operators[i + 2]
                + annihilation_operators[i] @ creation_operators[i + 2]
            )

        # hopping factors for down-spin fermions (operators with odd indices)
        for i in range(1, 2 * n - 2, 2):
            hopping_operators.append(
                -1 * creation_operators[i] @ annihilation_operators[i + 2]
                + annihilation_operators[i] @ creation_operators[i + 2]
            )

        # adds hopping between final and initial sites if lattice has pbc and has num_sites > 2
        if self.lattice.has_pbc() and self.lattice.get_num_sites() > 2:
            hopping_operators.append(
                -1 * creation_operators[2 * n - 2] @ annihilation_operators[0]
                + annihilation_operators[2 * n - 2] @ creation_operators[0]
            )
            hopping_operators.append(
                -1 * creation_operators[2 * n - 1] @ annihilation_operators[1]
                + annihilation_operators[2 * n - 1] @ creation_operators[1]
            )

        self.hopping_term = t * sum(hopping_operators)
        self.hamiltonian = self.hopping_term + self.interaction_term

    def __str__(self):
        return (
            "t = "
            + str(self.t)
            + "\nU = "
            + str(self.U)
            + "\n\nLattice:\n"
            + str(self.lattice)
            + "\n\n"
            + str(self.hamiltonian)
        )

    def jw_interaction_term(self):
        """Applies Jordan-Wigner transformation to HubbardHamiltonian interaction term and returns result.

        Returns
        -------
        jw_interaction : SparsePauliOp
            The qiskit-nature representation of the Jordan-Wigner transformed interaction term.

        """
        jw = JordanWignerMapper()
        jw_interaction = jw.map(self.interaction_term)
        return jw_interaction

    def jw_hopping_term(self):
        """Applies Jordan-Wigner transformation to HubbardHamiltonian hopping term and returns result.

        Returns
        -------
        jw_hopping : SparsePauliOp
            The qiskit-nature representation of the Jordan-Wigner transformed hopping term.

        """
        jw = JordanWignerMapper()
        jw_hopping = jw.map(self.hopping_term)
        return jw_hopping

    def jw_hamiltonian(self):
        """Applies Jordan-Wigner transformation to HubbardHamiltonian and returns result.

        Returns
        -------
        jw_hamiltonian : SparsePauliOp
            The qiskit-nature representation of the Jordan-Wigner transformed Hamiltonian.

        """
        jw = JordanWignerMapper()
        jw_hamiltonian = jw.map(self.hamiltonian)
        return jw_hamiltonian

    def get_t_value(self):
        """Access method for acquiring the hopping strength of a HubbardHamiltonian.

        Returns
        -------
        t : float
            The value of the hopping strength for the HubbardHamiltonian.
        """
        return self.t

    def get_U_value(self):
        """Access method for acquiring the interaction strength of a HubbardHamiltonian.

        Returns
        -------
        U : float
            The value of the interaction strength for the HubbardHamiltonian.
        """
        return self.U

    def get_lattice(self):
        """Access method for acquiring the Lattice associated with a HubbardHamiltonian.

        Returns
        -------
        lattice : Lattice
            The Lattice object described by the HubbardHamiltonian.
        """
        return self.lattice

    def get_interaction_term(self):
        """Access method for acquiring the interaction term of a HubbardHamiltonian in FermionicOp form.

        Returns
        -------
        interaction_term : FermionicOp
            The FermionicOp representation of the interaction term.
        """
        return self.interaction_term

    def get_hopping_term(self):
        """Access method for acquiring the hopping term of a HubbardHamiltonian in FermionicOp form.

        Returns
        -------
        hopping_term : FermionicOp
            The FermionicOp representation of the hopping term.
        """
        return self.hopping_term

    def get_hamiltonian(self):
        """Access method for acquiring the HubbardHamiltonian in FermionicOp form.

        Returns
        -------
        hamiltonian : FermionicOp
            The FermionicOp representation of the HubbardHamiltonian.
        """
        return self.hamiltonian
