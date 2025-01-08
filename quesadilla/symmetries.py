import numpy as np
from numpy.typing import ArrayLike
from pymatgen.core.structure import Structure

import quesadilla.espresso_symm as espresso_symm


class Symmetrizer:
    def __init__(self, structure: Structure, threshold: float = 1e-5):
        """Initialize a symmetry analysis object for a structure.

        The constructor takes the pymatgen Structure object and sets up attributes in
        a FORTRAN friendly manner. Its various methods call the Quantum ESPRESSO
        symmetry routines to compute the Fourier-transformed force constants in the star of a given q-point, and then symmetrize them using space group symmetries.

        Args:
            structure: The crystal structure to analyze.
            threshold: Numerical threshold for symmetry detection, defaults to 1e-5.
        """

        # Structure
        self.structure = structure

        # Define QE Variables with fortran-ready data types
        # NOTE: we use the same names as the QE Fortran routines
        self.nat = np.intc(len(structure))  # Number of atoms
        # s(:,:,i) is the i-th symmetry operation. Zeros for missing symmetries
        self.s = np.zeros((3, 3, 48), dtype=np.intc, order="F")
        # irt(i,j) is the index of the atom you get from applying symmetry i to atom j
        self.irt = np.zeros((48, self.nat), dtype=np.intc, order="F")
        # Array with index of inverse symmetry of each operation
        self.invs = np.zeros((48), dtype=np.intc, order="F")
        # Array with S.r_a - r_a for each symm op and atom a
        self.rtau = np.zeros((3, 48, self.nat), dtype=np.float64, order="F")
        # ft(:, j) is the fractional translation of symmetry j in FRAC coords
        self.ft = np.zeros((3, 48), dtype=np.float64, order="F")
        # Number of symmetries in the crystal
        self.nsym = np.intc(0)
        # Number of symmetries in the Bravais lattice
        self.nrot = np.intc(0)

        # Whether there is a symm op S_{-} such that S_{-} . q = -q + G
        self.minus_q = np.array([True], dtype=np.intc)
        # s(:, :, irotmq) is the symmetry operation that maps q to -q + G
        self.irotmq = np.intc(0)
        # This is the G vector in S_{-} . q = -q + G
        self.gimq = np.zeros((3), dtype=np.float64, order="F")
        # gi(:, i) is the G vector in S(:,:,i)@q = q + G
        self.gi = np.zeros((3, 48), dtype=np.float64, order="F")
        # Number of symmetries in the star of a q-point
        self.nsymq = np.intc(0)

        # Dictionary of atomic symbol : unique number
        unique_species = list({site.species_string for site in structure})
        symbs = {symb: i + 1 for i, symb in enumerate(unique_species)}
        # Unique species index for each site
        self.ityp = np.zeros(self.nat, dtype=np.intc)
        for i, site in enumerate(structure):
            self.ityp[i] = symbs[site.species_string]
        # Atomic positions in cartesian coordinates
        self.tau = np.array(self.structure.cart_coords.T, dtype=np.float64, order="F")

        # Lattice vectors (QE routines expect columns as basis vectors)
        self.at = np.array(structure.lattice.matrix.T, dtype=np.float64, order="F")
        # Reciprocal lattice vectors (QE routines expect columns as basis vectors)
        self.bg = np.array(
            structure.lattice.reciprocal_lattice.matrix.T / (2 * np.pi),
            dtype=np.float64,
            order="F",
        )

    def get_xq_from_aq(self, aq: ArrayLike):
        """
        Get the q-point in Cartesian coordinates / (2*pi) from a q-point in fractional coordinates.

        Parameters:
            - aq : ndarray(3)
                The q-point in fractional coordinates.

        Returns:
            ndarray(3)
                The q-point in Cartesian coordinates / (2*pi).
        """
        return np.array(
            self.structure.lattice.reciprocal_lattice.get_cartesian_coords(aq)
            / (2 * np.pi),
            dtype=np.float64,
            order="F",
        )

    def get_aq_from_xq(self, xq: ArrayLike):
        """
        Get the q-point in fractional coordinates from a q-point in Cartesian coordinates / (2*pi).

        Parameters:
            - xq : ndarray(3)
                The q-point in Cartesian coordinates / (2*pi).

        Returns:
            ndarray(3)
                The q-point in fractional coordinates.
        """

        return np.array(
            self.structure.lattice.reciprocal_lattice.get_fractional_coords(
                xq * (2 * np.pi)
            ),
            dtype=np.float64,
            order="F",
        )

    def _setup_lattice_symmetries(self, verbose: bool = False):
        """
        Sets up the symmetries of the bravais lattice.

        Parameters:
            verbose : bool
                Whether to print the number of symmetries of the bravais lattice.
        """

        # Sets up the symmetries of the bravais lattice
        espresso_symm.symm_base.set_sym_bl(self.at)

        # TODO: make these return values instead of modifying module vars
        self.s = np.copy(espresso_symm.symm_base.s)
        self.ft = np.copy(espresso_symm.symm_base.ft)
        self.nrot = espresso_symm.symm_base.nrot

        if verbose:
            print("Symmetries of the bravais lattice:", self.nrot)

    def _setup_crystal_symmetries(self, verbose: bool = False):
        """
        Sets up the symmetries of the crystal.

        Parameters:
            verbose : bool
                Whether to print the number of symmetries of the crystal.
        """
        # TODO: implement magnetism (currently just a dummy variable)
        # TODO: some lines in symm_base need to be uncommented/checked
        m_loc = np.zeros((3, self.nat), dtype=np.float64, order="F")
        nspin_mag = 1

        # Find the symmetries of the crystal
        espresso_symm.symm_base.set_sym(
            self.at, self.bg, self.tau, self.ityp, nspin_mag, m_loc
        )

        # TODO: make these return values instead of modifying module vars
        self.s = np.copy(espresso_symm.symm_base.s)
        self.ft = np.copy(espresso_symm.symm_base.ft)
        self.nsym = espresso_symm.symm_base.nsym
        self.crystal_s = np.copy(self.s)
        self.crystal_invs = np.copy(self.invs)
        self.crystal_irt = np.copy(self.irt)

        if verbose:
            print("Symmetries of the crystal:", self.nsym)

    def _setup_little_group(self, xq: ArrayLike, verbose: bool = False):
        """
        Sets up the little group (small group) of a q-point.

        Parameters:
            - xq : ndarray(3), q point in CARTESIAN coordinates / (2 * pi)
        """
        # Set up the symmetries of the lattice
        self._setup_lattice_symmetries(verbose)
        # Set up the symmetries of the crystal
        self._setup_crystal_symmetries(verbose)
        # ~~~~~~~~~~PART 1~~~~~~~~~~
        # Set up the symmetries of the small group of q
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~
        syms = np.zeros((48), dtype=np.intc)
        syms[: self.nsym] = np.intc(1)

        # TODO: make these return values instead of modifying module vars
        espresso_symm.smallg_q(xq, 0, self.at, self.nsym, self.s, syms, self.minus_q)
        # Copy the symmetries of the small group of q
        self.nsymq = espresso_symm.symm_base.copy_sym(
            espresso_symm.symm_base.nsym, syms
        )
        # Recompute the inverses as the order of sym.ops. has changed
        espresso_symm.symm_base.inverse_s()
        # Copy stuff into python
        self.s = np.copy(espresso_symm.symm_base.s)
        self.invs = np.copy(espresso_symm.symm_base.invs)
        self.ft = np.copy(espresso_symm.symm_base.ft)
        self.irt = np.copy(espresso_symm.symm_base.irt)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # PART 2: figure out the q -> -q+G symmetries
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        lgamma = np.allclose(xq, [0, 0, 0], rtol=0, atol=1e-5)
        self.irotmq, self.gi, self.gimq = espresso_symm.set_giq(
            xq, lgamma, self.bg, self.at, self.s, self.nsymq, self.nsym, self.minus_q
        )
        # Part 3: Compute rtau = S . r_a - r_a
        self.rtau = espresso_symm.sgam_lr(
            self.at, self.bg, self.nsym, self.s, self.irt, self.tau
        )

        # Print some info
        if verbose:
            q = self.structure.lattice.reciprocal_lattice.get_fractional_coords(
                xq * 2 * np.pi
            )
            print(f"Symmetries of the small group of q = {np.round(q, 5)}:", self.nsymq)
            if self.minus_q:
                gimq = self.structure.lattice.reciprocal_lattice.get_fractional_coords(
                    self.gimq * 2 * np.pi
                )
                print(
                    (
                        f"in addition sym. q -> -q+G, with irotmq = {self.irotmq}"
                        f" and G = {np.round(gimq, 3)}"
                    )
                )

    def get_star_q(self, xq: ArrayLike, verbose: bool = False, debug: bool = False):
        # sourcery skip: extract-method
        """
        Get the star of a q-point.

        Parameters:
            - xq : ndarray(3), q point in CARTESIAN coordinates / (2 * pi)
            - verbose : bool
                Whether to print the number of q-points in the star and the list of q-points in the star.
            - debug : bool
                Whether to print debug information from Fortran routines.

        Returns:
            ndarray(nq, 3), the q-points in the star in CARTESIAN coordinates / (2 * pi)
        """
        # To get the star we need the symmetries of the WHOLE CRYSTAL
        # So that's just the litle group of the Gamma point
        self._setup_little_group(np.array([0, 0, 0]), verbose)
        # Return values are:
        # 1. number of q-points in the star
        # 2. The vector that q is mapped to under each symmop
        # 3. Index of q in sxq
        # 4. Index of -q in sxq, 0 if not present
        self.nqs, self.sxq, self.isq, self.imq = espresso_symm.star_q(
            xq,  # q point in cartesian coordinates/2pi
            self.at,  # lattice vectors (see __init__ for format)
            self.bg,  # rec. lattice vectors (see __init__ for format)
            self.nsym,  # Number of symmetries in the small group of q
            self.s,  # Array of ALL symmetries of the crystal
            self.invs,  # Index of inverse of s in self.s
            False,  # Debug flag
        )
        star_xq = self.sxq[:, : self.nqs].T
        if verbose:
            q_star = np.array([self.get_aq_from_xq(xqs) for xqs in star_xq])
            print("Number of q in the star:", self.nqs)
            print("List of q in the star:")
            for i, qq in enumerate(q_star):
                print(f"{i+1}    {np.round(qq, 8)}")
            if self.imq == 0:
                print("In addition, there is the -q list, which is NOT in the star:")
                for i, qq in enumerate(q_star):
                    print(f"{i+1+self.nqs}    {np.round(qq, 8)}")
            else:
                mq = self.get_aq_from_xq(self.sxq[:, self.imq - 1])  # -q
                gmq = self.get_aq_from_xq(self.gimq)  # G
                q = self.get_aq_from_xq(xq)  # q
                print("-q is also in the star: ", np.round(mq, 8))
                print("With G = ", np.round(gmq, 8))
                print("So that S_ @ q - (-q + G)", np.round(mq + q - gmq, 5))

        return star_xq

    def symmetrize_fcq(self, fcq: ArrayLike, xq: ArrayLike, verbose: bool = False):
        """
        Symmetrize the force constants at a q-point.

        Parameters:
            - fcq : ndarray(3 * nat, 3 * nat), in FRACTIONAL coordinates
            - q : ndarray(3)
                The q vector, in CARTESIAN coordinates / (2 * pi)

        Returns:
            ndarray(3 * nat, 3* nat): The symmetrized force constants at q in FRACTIONAL coordinates.
        """
        self._setup_little_group(xq, verbose)
        return espresso_symm.symdynph_gq_new(
            xq,
            self.at,
            self.bg,
            fcq,
            self.s,
            self.invs,
            self.rtau,
            self.irt,
            self.nsymq,
            self.irotmq,
            self.minus_q,
        )

    def get_fcq_in_star(self, fcq, aq, verbose=True):
        """
        Get the force constants in the star of a q-point.

        Parameters:
            - fcq : ndarray(3 * nat, 3 * nat), Fourier-transformed force constants
                        in FRACTIONAL coordinates
            - q : ndarray(3)
                The q vector, in FRACTIONAL coordinates

        Returns:
            dict: The force constants in the star of q. Keys are the q-points in fractional coordinates. Values are the force constants in the star in the form of a 3*nat x 3*nat matrix in FRACTIONAL coordinates.
        """

        # Symmetrize the force constants at q
        aq = np.array(aq, dtype=np.float64, order="F")
        xq = self.get_xq_from_aq(aq)
        fcq = np.array(fcq, dtype=np.complex128, order="F")
        fcq_symm = self.symmetrize_fcq(fcq, xq, verbose)

        # Get star of the q-point
        star_xq = self.get_star_q(xq, verbose)
        # In the special case that -q is NOT in the star,
        # We compute the FCQ at -q anyway using TRS and include it
        if self.imq == 0:
            nq_tot = 2 * self.nqs
            star_xq = np.concatenate((star_xq, -star_xq), axis=0)
        else:
            nq_tot = self.nqs

        # Get the FC(q) in the star
        fcq_star = espresso_symm.q2qstar_ph(
            fcq_symm,
            self.at,
            self.bg,
            self.nsym,
            self.s,
            self.invs,
            self.irt,
            self.rtau,
            self.nqs,
            self.sxq,
            self.isq,
            self.imq,
            nq_tot,
        )

        # Turn into a dictionary
        final_fcq = {}
        for i, xqq in enumerate(star_xq):
            qq = self.get_aq_from_xq(xqq)
            final_fcq[tuple(qq)] = fcq_star[i, :, :]

        return final_fcq
