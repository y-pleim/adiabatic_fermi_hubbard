from qiskit_nature.second_q.hamiltonians.lattices import BoundaryCondition, LineLattice


class Lattice:
    """A class for representing a 1D lattice of fermions.

    Parameters
    ----------
    num_sites : int, default: 2
        The number of lattice sites.
    bc : int, default: 1
        Specifies boundary conditions; 1 for periodic, 0 for open.
    """

    def __init__(self, num_sites: int = 2, bc: int = 1):
        self.num_sites = num_sites

        if bc == 0:
            self.pbc = False
            self.lattice = LineLattice(num_sites, 1.0, 0.0, BoundaryCondition.OPEN)
        elif bc == 1:
            self.pbc = True
            self.lattice = LineLattice(num_sites, 1.0, 0.0, BoundaryCondition.PERIODIC)
        else:
            raise ValueError(
                "Invalid boundary condition specification. Must be 0 for open boundary conditions or 1 for periodic boundary conditions."
            )

    def __str__(self):
        return (
            "Number of sites: "
            + str(self.get_num_sites())
            + " sites, \nPeriodic boundary conditions: "
            + str(self.has_pbc())
            + "."
        )

    def get_num_sites(self):
        """Access method which returns number of lattice sites.

        Returns
        -------
        num_sites : int
            The number of lattice sites.
        """
        return self.num_sites

    def has_pbc(self):
        """Access method which indicates if Lattice has periodic boundary conditions.

        Returns
        -------
        pbc : bool
            Periodic boundary conditions flag; if false, boundary conditions are open.
        """
        return self.pbc

    def get_qiskit_object(self):
        """Access method which returns the qiskit LineLattice object associated with the Lattice instance.

        Returns
        -------
        lattice : LineLattice
            qiskit-nature LineLattice object with the same number of sites and periodic boundary conditions stored in Lattice instance.
        """
        return self.lattice
