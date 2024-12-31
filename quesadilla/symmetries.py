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
        espresso_symm.symm_base.set_at_bg(self.at, self.bg)

        # Prepare the symmetries
        espresso_symm.symm_base.set_sym_bl()

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
        self.s = np.copy(espresso_symm.symm_base.s)
        self.ft = np.copy(espresso_symm.symm_base.ft)
        self.nsym = espresso_symm.symm_base.nsym

        if verbose:
            print("Symmetries of the crystal:", self.nsym)

    def setup_little_cogroup(self, q: ArrayLike, verbose: bool = False):
        """
        Get symmetries of the small group of q

        Setup the symmetries in the small group of Q.

        Parameters
        ----------
            q_point : ndarray
            q vector in FRACTIONAL COORDINATES
        """
        q = np.array(q, dtype=np.float64, order="F")

        # Setup the bravais lattice
        # self.setup_lattice_symmetries(verbose)

        # ------- CRYSTAL SYMMETRIES -------
        # self.setup_crystal_symmetries(verbose)

        # -----------------------------------

        # ------- SMALL GROUP OF Q -------

        #! part 2: this computes gi, gimq
        #!
        #! finally this does some of the above again and also computes rtau...
        #! TODO: provide this routine
        # ALLOCATE(rtau( 3, 48, nat))
        # CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
        # -----------------------------
        # ~~~~~~~~~~PART 1~~~~~~~~~~
        # minus_q = .true.
        # sym = .false.
        # sym(1:nsym) = .true.
        # CALL smallg_q(xq, 0, at, bg, nsym, s, sym, minus_q)
        # nsymq = copy_sym(nsym, sym) ! TODO: provide this routine
        #! recompute the inverses as the order of sym.ops. has changed
        # CALL inverse_s ( )
        syms = np.zeros((48), dtype=np.intc)
        syms[: self.nsym] = np.intc(1)

        xq = self.structure.lattice.reciprocal_lattice.get_cartesian_coords(q) / (
            2 * np.pi
        )
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
        # ~~~~~~~~~~~~~~~~~~~~
        # PART 2: figure out the q -> -q+G symmetries
        # call set_giq (xq,s,nsymq,nsym,minus_q,gi,gimq,lgamma) ! TODO: provide this routine
        # WRITE(*, '(5x,a,i3)') "Symmetries of small group of q:", nsymq
        # IF(minus_q) WRITE(*, '(10x,a)') "in addition sym. q -> -q+G:"
        lgamma = np.allclose(q, [0, 0, 0], rtol=0, atol=1e-5)
        # set_giq (xq,lgamma,bg,at,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
        self.irotmq, self.gi, self.gimq = espresso_symm.set_giq(
            xq, lgamma, self.bg, self.at, self.s, self.nsymq, self.nsym, self.minus_q
        )
        if verbose:
            print(f"Symmetries of the small group of q = {np.round(q, 5)}:", self.nsymq)
            if self.minus_q:
                print("in addition sym. q -> -q+G:")
                print(f"irotmq = {self.irotmq}")
                gimq = self.structure.lattice.reciprocal_lattice.get_fractional_coords(
                    self.gimq * 2 * np.pi
                )
                print(f"gi = {np.round(gimq, 3)}")
        # ~~~~~~~~~~~~~~~~~~~~
        # PART 3: re-set up the symmetries and compute rtau
        # --------------------------------
        # sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
        self.rtau = espresso_symm.sgam_lr(
            self.at, self.bg, self.nsym, self.s, self.irt, self.tau
        )

    def get_star_q(self, q: ArrayLike, verbose: bool = False, debug: bool = False):
        """
        GET THE Q STAR
        ==============

        Given a vector in q space, get the whole star.
        We use the quantum espresso subrouitine.

        Parameters
        ----------
            q_vector : ndarray(size= 3, dtype = np.float64)
                The q vector

        Results
        -------
            q_star : ndarray(size = (nq_star, 3), dtype = np.float64)
                The complete q star
        """
        #!
        # CALL star_q(xq, at, bg, nsym, s, invs, nqs, sxq, isq, imq, .true. ) ! TODO: provide
        q = np.array(q, dtype=np.float64, order="F")
        self.setup_lattice_symmetries(verbose)
        self.setup_crystal_symmetries(verbose)
        crystal_symmetries = np.copy(self.s)
        crystal_invs = np.copy(self.invs)

        self.setup_little_cogroup(q, verbose)

        xq = np.array(
            self.structure.lattice.reciprocal_lattice.get_cartesian_coords(q)
            / (2 * np.pi),
            dtype=np.float64,
            order="F",
        )

        # Returns:
        # 1. number of q-points in the star
        # 2. The vector that q is mapped to under each symmop
        # 3. Index of q in sxq
        # 4. Index of -q in sxq, 0 if not present
        nqs, sxq, isq, imq = espresso_symm.star_q(
            xq,  # q point in cartesian coordinates/2pi
            self.at,  # lattice vectors (see __init__ for format)
            self.bg,  # rec. lattice vectors (see __init__ for format)
            self.nsym,  # Number of symmetries in the small group of q
            crystal_symmetries,  # Array of ALL symmetries of the crystal
            crystal_invs,  # Index of inverse of s in self.s
            debug,  # Verbosity flag
        )
        xq_star = sxq[:, :nqs].T
        q_star = np.array(
            [
                self.structure.lattice.reciprocal_lattice.get_fractional_coords(
                    xqs * (2 * np.pi)
                )
                for xqs in xq_star
            ]
        )
        if verbose:
            print("Number of q in the star:", nqs)
            print("List of q in the star:")
            for i, qq in enumerate(q_star):
                print(f"{i+1}    {np.round(qq, 8)}")
            if imq == 0:
                print("In addition, there is the -q list, which is NOT in the star:")
                for i, qq in enumerate(q_star):
                    print(f"{i+1+nqs}    {np.round(qq, 8)}")
            else:
                mxq = sxq[:, imq - 1]
                mq = self.structure.lattice.reciprocal_lattice.get_fractional_coords(
                    mxq * (2 * np.pi)
                )
                gmq = self.structure.lattice.reciprocal_lattice.get_fractional_coords(
                    self.gimq * (2 * np.pi)
                )
                print("-q is also in the star: ", np.round(mq, 8))
                print("With G = ", np.round(gmq, 8))
                print("So that S_ @ q - (-q + G)", np.round(mq + q - gmq, 5))

        return q_star

    def get_fcq_in_star(self, fcq, q):
        """
        APPLY THE Q STAR SYMMETRY
        =========================

        Given the fc matrix at each q in the star, it applies the symmetries in between them.

        Parameters
        ----------
            - fcq : ndarray(nq, 3xnat, 3xnat)
                The dynamical matrices for each q point in the star
            - q_star : ndarray(nq, 3)
                The q vectors that belongs to the same star
        """

        # TODO: REWRITE like this
        # DO i = 1, 3 * nat
        # na = (i - 1) / 3 + 1
        # icar = i - 3 * (na - 1)
        # DO j = 1, 3 * nat
        #    nb = (j - 1) / 3 + 1
        #    jcar = j - 3 * (nb - 1)
        #    d2 (i, j) = phi(icar, jcar, na, nb)
        # ENDDO
        # ENDDO
        #!
        #!if (imq == 0) then
        #!    nq_tot = 2 * nq
        #! else
        #!    nq_tot = nq
        #!end if
        # CALL q2qstar_ph (d2, at, bg, nat, nsym, s, invs, irt, rtau, &
        #                nqs, sxq, isq, imq, 1) ! TODO: provide this routine
        # Setup all the symmetries
        q = np.array(q, dtype=np.float64, order="F")
        q_star = self.get_star_q(q, verbose=False)
        nq_total = np.shape(q_star)[0]
        # print("Got total star:", q_star)

        # q = np.array(q, dtype=np.float64, order="F")
        # Convert to cartesian coordinates
        # q = np.array(
        #    self.structure.lattice.reciprocal_lattice.get_cartesian_coords(q)
        #    / (2 * np.pi),
        #    dtype=np.float64,
        #    order="F",
        # )
        # self.setup_little_cogroup(q, verbose=True)
        nq_no_mq, sxq, isq, imq = espresso_symm.star_q(
            q_star[0],
            self.at,
            self.bg,
            self.nsymq,
            self.s,
            self.invs,
            1,
        )
        # print("I am getting nq_new vs nq:", nq_new, nq)
        # print("I am getting isq:", isq)
        # print("I am getting imq:", imq)
        # print("I am getting sxq:", sxq.T)
        # assert nq_new == nq

        fcq = np.array(fcq, dtype=np.complex128, order="F")
        dyn_star = np.zeros(
            (nq_total, 3, 3, self.nat, self.nat), dtype=np.complex128, order="F"
        )
        print("self.nsym:", self.nsym)
        dyn_star = espresso_symm.q2qstar_out(
            fcq,  # FC matrix at q
            self.at,  # lattice vectors
            self.bg,  # reciprocal lattice vectors
            self.nsym,  # Number of symmetries in whole crystal
            self.s,  # All symmetries in the crystal
            self.invs,  # Index of inverse of s in self.s
            self.irt,  # Index of atom you get from apply S to atom a
            self.rtau,  # Array with S.r_a - r_a for each symm op and atom a
            nq_no_mq,  # Number of q points in the star
            sxq,  # S.q for each symmetry in the crystal
            isq,  # Index of S.q in sxq for each S
            imq,  # Index of -q in sxq for each S
            nq_total,  # Will be 2*nq_no_mq only if -q is NOT in the star
            nat=self.nat,  # Number of atoms in the crystal
        )

        dynmats = {}
        for i, q in enumerate(q_star):
            q = np.round(
                self.structure.lattice.reciprocal_lattice.get_fractional_coords(
                    q * (2 * np.pi)
                ),
                decimals=6,
            )
            D_blocks = np.zeros((self.nat, self.nat, 3, 3), dtype=np.complex128)
            for na in range(self.nat):
                for nb in range(self.nat):
                    D_blocks[na, nb] = dyn_star[i, :, :, na, nb]
            dynmats[tuple(q)] = D_blocks.swapaxes(1, 2).reshape(
                3 * self.nat, 3 * self.nat
            )

        return dynmats

    # def SymmetrizeDynQ(self, dyn_matrix, q_point):
    #    """
    #    DYNAMICAL MATRIX SYMMETRIZATION
    #    ===============================

    #    Use the Quantum ESPRESSO fortran code to symmetrize the dynamical matrix
    #    at the given q point.

    #    NOTE: the symmetries must be already initialized.

    #    Parameters
    #    ----------
    #        dyn_matrix : ndarray (3nat x 3nat)
    #            The dynamical matrix associated to the specific q point (cartesian coordinates)
    #        q_point : ndarray 3
    #            The q point related to the dyn_matrix.

    #    The input dynamical matrix will be modified by the current code.
    #    """
    ######################### star of q #########################
    # TODO: replace with this
    # do na = 1, nat
    # do nb = 1, nat
    #    call trntnsc (phi (1, 1, na, nb), at, bg, - 1) ! TODO: provide this routine
    # enddo
    # enddo
    # CALL symdynph_gq_new (xq, phi, s, invs, rtau, irt, nsymq, nat, &
    #    irotmq, minus_q) ! TODO: provide this routine
    # do na = 1, nat
    # do nb = 1, nat
    #    call trntnsc (phi (1, 1, na, nb), at, bg, + 1) ! TODO: provide this routine
    # enddo
    # enddo

    #    # TODO: implement hermitianity to speedup the conversion

    #    # Prepare the array to be passed to the fortran code
    #    QE_dyn = np.zeros(
    #        (3, 3, self.QE_nat, self.QE_nat), dtype=np.complex128, order="F"
    #    )

    #    # Get the crystal coordinates for the matrix
    #    for na in range(self.QE_nat):
    #        for nb in range(self.QE_nat):
    #            fc = dyn_matrix[3 * na : 3 * na + 3, 3 * nb : 3 * nb + 3]
    #            QE_dyn[:, :, na, nb] = Methods.convert_matrix_cart_cryst(
    #                fc, self.structure.unit_cell, False
    #            )

    #    # Prepare the xq variable
    #    # xq = np.ones(3, dtype = np.float64)
    #    xq = np.array(q_point, dtype=np.float64)

    #    # USE THE QE library to perform the symmetrization
    #    espresso_symm.symdynph_gq_new(
    #        xq,
    #        QE_dyn,
    #        self.QE_s,
    #        self.QE_invs,
    #        self.QE_rtau,
    #        self.QE_irt,
    #        self.QE_irotmq,
    #        self.QE_minus_q,
    #        self.QE_nsymq,
    #        self.QE_nat,
    #    )

    #    # Return to cartesian coordinates
    #    for na in range(self.QE_nat):
    #        for nb in range(self.QE_nat):
    #            fc = QE_dyn[:, :, na, nb]
    #            dyn_matrix[3 * na : 3 * na + 3, 3 * nb : 3 * nb + 3] = (
    #                Methods.convert_matrix_cart_cryst(
    #                    fc, self.structure.unit_cell, True
    #                )
    #            )
