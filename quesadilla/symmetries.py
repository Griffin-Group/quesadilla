import numpy as np
from numpy.typing import ArrayLike
from pymatgen.core.structure import Structure

# from pymatgen.util.coord import pbc_diff
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

        # Setup the threshold
        self.threshold = threshold
        espresso_symm.symm_base.set_accep_threshold(self.threshold)

        # Define QE Variables with fortran-ready data types
        # NOTE: we use the same names as the QE Fortran routines
        self.nat = np.intc(len(structure))  # Number of atoms
        # Array with symmetry matrices in fractional coordinates
        self.s = np.zeros((3, 3, 48), dtype=np.intc, order="F")
        # TODO: array with ??
        self.irt = np.zeros((48, self.nat), dtype=np.intc, order="F")
        # Array with index of inverse symmetry of each operation
        self.invs = np.zeros((48), dtype=np.intc, order="F")
        # Array with S.r_a - r_a for each symm op and atom a
        self.rtau = np.zeros((3, 48, self.nat), dtype=np.float64, order="F")
        # TODO: Array with ??
        self.ft = np.zeros((3, 48), dtype=np.float64, order="F")
        # Arrays for G vectors
        self.gi = np.zeros((3, 48), dtype=np.float64, order="F")
        self.gimq = np.zeros((3), dtype=np.float64, order="F")

        self.minus_q = np.array([True], dtype=np.intc)
        self.irotmq = np.intc(0)
        self.nsymq = np.intc(0)
        self.nsym = np.intc(0)

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

    def setup_lattice_symmetries(self, verbose: bool = False):
        # TODO: make these input values instead of modifying module vars
        espresso_symm.symm_base.set_at_bg(self.at, self.bg)

        # Prepare the symmetries
        espresso_symm.symm_base.set_sym_bl()

        # TODO: make these return values instead of modifying module vars
        self.s = np.copy(espresso_symm.symm_base.s)
        self.ft = np.copy(espresso_symm.symm_base.ft)
        self.nsym = espresso_symm.symm_base.nrot

        if verbose:
            print("Symmetries of the bravais lattice:", self.nsym)

    def setup_crystal_symmetries(self, verbose: bool = False):
        # TODO: implement magnetism (currently just a dummy variable)
        # TODO: some lines in symm_base need to be uncommented/checked
        m_loc = np.zeros((3, self.nat), dtype=np.float64, order="F")

        # Find the symmetries of the crystal
        espresso_symm.symm_base.find_sym(self.tau, self.ityp, False, m_loc, False)

        # Copy into python
        # TODO: make these return values instead of modifying module vars
        self.s = np.copy(espresso_symm.symm_base.s)
        self.ft = np.copy(espresso_symm.symm_base.ft)
        self.nsym = espresso_symm.symm_base.nsym
        self.crystal_s = np.copy(self.s)
        self.crystal_invs = np.copy(self.invs)
        self.crystal_irt = np.copy(self.irt)

        if verbose:
            print("Symmetries of the crystal:", self.nsym)

    def setup_little_cogroup(self, xq: ArrayLike, verbose: bool = False):
        """
        Sets up the little co-group (small group) of a q-point.

        Parameters:
            - xq : ndarray(3), q point in CARTESIAN coordinates / (2 * pi)
        """
        # ~~~~~~~~~~PART 1~~~~~~~~~~
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
        # set_giq (xq,lgamma,bg,at,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
        self.irotmq, self.gi, self.gimq = espresso_symm.set_giq(
            xq, lgamma, self.bg, self.at, self.s, self.nsymq, self.nsym, self.minus_q
        )
        if verbose:
            q = self.structure.lattice.reciprocal_lattice.get_fractional_coords(
                xq * 2 * np.pi
            )
            print(f"Symmetries of the small group of q = {np.round(q, 5)}:", self.nsymq)
            if self.minus_q:
                print("in addition sym. q -> -q+G:")
                print(f"irotmq = {self.irotmq}")
                gimq = self.structure.lattice.reciprocal_lattice.get_fractional_coords(
                    self.gimq * 2 * np.pi
                )
                print(f"gi = {np.round(gimq, 3)}")

    def get_star_q(self, xq: ArrayLike, verbose: bool = False, debug: bool = False):
        # sourcery skip: extract-method
        """
        Get the star of a q-point.
        """
        # Returns:
        # 1. number of q-points in the star
        # 2. The vector that q is mapped to under each symmop
        # 3. Index of q in sxq
        # 4. Index of -q in sxq, 0 if not present
        self.nqs, self.sxq, self.isq, self.imq = espresso_symm.star_q(
            xq,  # q point in cartesian coordinates/2pi
            self.at,  # lattice vectors (see __init__ for format)
            self.bg,  # rec. lattice vectors (see __init__ for format)
            self.nsym,  # Number of symmetries in the small group of q
            self.crystal_s,  # Array of ALL symmetries of the crystal
            self.crystal_invs,  # Index of inverse of s in self.s
            debug,  # Verbosity flag
        )
        star_xq = self.sxq[:, : self.nqs].T
        if verbose:
            q_star = np.array(
                [
                    self.structure.lattice.reciprocal_lattice.get_fractional_coords(
                        xqs * (2 * np.pi)
                    )
                    for xqs in star_xq
                ]
            )
            print("Number of q in the star:", self.nqs)
            print("List of q in the star:")
            for i, qq in enumerate(q_star):
                print(f"{i+1}    {np.round(qq, 8)}")
            if self.imq == 0:
                print("In addition, there is the -q list, which is NOT in the star:")
                for i, qq in enumerate(q_star):
                    print(f"{i+1+self.nqs}    {np.round(qq, 8)}")
            else:
                mxq = self.sxq[:, self.imq - 1]
                mq = self.structure.lattice.reciprocal_lattice.get_fractional_coords(
                    mxq * (2 * np.pi)
                )
                gmq = self.structure.lattice.reciprocal_lattice.get_fractional_coords(
                    self.gimq * (2 * np.pi)
                )
                q = self.structure.lattice.reciprocal_lattice.get_fractional_coords(
                    xq * (2 * np.pi)
                )
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
            ndarray(3, 3, nat, nat): The symmetrized force constants at q.
        """
        # Compute rtau = S . r_a - r_a
        fcq_symm = np.zeros((3, 3, self.nat, self.nat), dtype=np.complex128, order="F")
        fcq_symm = espresso_symm.symdynph_gq_new(
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
        return fcq_symm

    def get_fcq_in_star(self, fcq, q, verbose=True):
        """
        TODO: docstring
        """

        xq = np.array(
            self.structure.lattice.reciprocal_lattice.get_cartesian_coords(q)
            / (2 * np.pi),
            dtype=np.float64,
            order="F",
        )
        fcq = np.array(fcq, dtype=np.complex128, order="F")
        self.setup_lattice_symmetries(verbose)
        self.setup_crystal_symmetries(verbose)
        self.setup_little_cogroup(np.array([0, 0, 0]), verbose)
        # -----------
        # fcq_symm = self.symmetrize_fcq(fcq, xq, verbose)
        fcq_symm = fcq
        self.rtau = espresso_symm.sgam_lr(
            self.at, self.bg, self.nsym, self.s, self.irt, self.tau
        )
        # self.rtau *= 0 # Should work for CsCl
        # -----------
        # star_xq = self.get_star_q(xq, verbose)
        self.nqs, self.sxq, self.isq, self.imq = espresso_symm.star_q(
            xq,  # q point in cartesian coordinates/2pi
            self.at,  # lattice vectors (see __init__ for format)
            self.bg,  # rec. lattice vectors (see __init__ for format)
            self.nsym,  # Number of symmetries in the small group of q
            self.s,  # Array of ALL symmetries of the crystal
            self.invs,  # Index of inverse of s in self.s
            False,  # Verbosity flag
        )
        star_xq = self.sxq[:, : self.nqs].T
        # print("Got total star:", q_star)

        # subroutine q2qstar_ph(fcq, at, bg, nat, nsym, s, invs, irt, rtau, &
        #              nq, sxq, isq, imq, nq_tot, fcqstar)
        nq_tot = 2 * self.nqs if self.imq == 0 else self.nqs
        fcq_star = np.zeros(
            (nq_tot, 3, 3, self.nat, self.nat), dtype=np.complex128, order="F"
        )
        # print("self.nsym", self.nsym)
        # print(self.crystal_symmetries)
        # self.setup_lattice_symmetries()
        # self.setup_crystal_symmetries()
        # self.setup_little_cogroup([0, 0, 0])
        # self.imq = 0
        # FIXME: this is broken :(
        print(f"I have {self.nqs} q-points in the star with imq {self.imq}")
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
        for i, f in enumerate(fcq_star):
            print(f"for q point {i} in the star, we have:")
            print(f)

        # Get dynmat in the whole star
        final_fcq = {}
        for i, xqq in enumerate(star_xq):
            qq = np.round(
                self.structure.lattice.reciprocal_lattice.get_fractional_coords(
                    xqq * (2 * np.pi)
                ),
                decimals=6,
            )
            D_blocks = np.zeros((self.nat, self.nat, 3, 3), dtype=np.complex128)
            for na in range(self.nat):
                for nb in range(self.nat):
                    D_blocks[na, nb] = fcq_star[i, :, :, na, nb]
            final_fcq[tuple(qq)] = D_blocks.swapaxes(1, 2).reshape(
                3 * self.nat, 3 * self.nat
            )
            # final_fcq[tuple(qq)] = fcq_star[i]

        return final_fcq
