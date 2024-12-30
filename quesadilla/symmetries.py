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

        self.minus_q = False
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

        # ------- BRAVAIS LATTICE SYMMETRIES -------
        # Setup the bravais lattice
        espresso_symm.symm_base.set_at_bg(self.at, self.bg)

        # Prepare the symmetries
        espresso_symm.symm_base.set_sym_bl()

        if verbose:
            print("Symmetries of the bravais lattice:", espresso_symm.symm_base.nrot)
        # Now copy all the work initialized on the symmetries inside python
        self.s = np.copy(espresso_symm.symm_base.s)
        self.ft = np.copy(espresso_symm.symm_base.ft)
        self.nsym = espresso_symm.symm_base.nrot
        # ---------------------------

        # ------- CRYSTAL SYMMETRIES -------
        # TODO: implement magnetism (currently just a dummy variable)
        m_loc = np.zeros((3, self.nat), dtype=np.float64, order="F")

        # Find the symmetries of the crystal
        # TODO: this doesn't work for nonsymmorphic SGs?
        espresso_symm.symm_base.find_sym(self.tau, self.ityp, 6, 6, 6, False, m_loc)

        if verbose:
            print("Symmetries of the crystal:", espresso_symm.symm_base.nsym)

        # Now copy all the work initialized on the symmetries inside python
        self.s = np.copy(espresso_symm.symm_base.s)
        self.ft = np.copy(espresso_symm.symm_base.ft)
        # -----------------------------------

        # ------- SMALL GROUP OF Q -------
        syms = np.zeros((48), dtype=np.intc)

        # Initialize to true the symmetry of the crystal
        syms[: espresso_symm.symm_base.nsym] = np.intc(1)

        self.minus_q = espresso_symm.symm_base.smallg_q(q, 0, syms)
        self.nsymq = espresso_symm.symm_base.copy_sym(
            espresso_symm.symm_base.nsym, syms
        )
        self.nsym = espresso_symm.symm_base.nsym

        # Recompute the inverses
        espresso_symm.symm_base.inverse_s()

        if verbose:
            print(f"Symmetries of the small group of q = {np.round(q, 5)}:", self.nsymq)
        # --------------------------------

        # Assign symmetries
        self.s = np.copy(espresso_symm.symm_base.s)
        self.invs = np.copy(espresso_symm.symm_base.invs)
        self.ft = np.copy(espresso_symm.symm_base.ft)
        self.irt = np.copy(espresso_symm.symm_base.irt)

        # Compute the additional shift caused by fractional translations
        self.rtau = espresso_symm.sgam_ph_new(
            self.at,
            self.bg,
            espresso_symm.symm_base.nsym,
            self.s,
            self.irt,
            self.tau,
            self.nat,
        )

        # Convert q to cartesian coordinates
        # lgamma = np.allclose(q, [0, 0, 0])
        # self.irotmq = 0
        # if self.minus_q:
        #    xq = np.array(
        #        self.structure.lattice.reciprocal_lattice.get_cartesian_coords(q)
        #        / (2 * np.pi),
        #        dtype=np.float64,
        #        order="F",
        #    )
        #    self.irotmq = espresso_symm.set_irotmq(
        #        xq,
        #        self.s,
        #        self.nsymq,
        #        self.nsym,
        #        self.minus_q,
        #        self.bg,
        #        self.at,
        #        lgamma,
        #    )
        #    print("IROTMQ from QE:", self.irotmq)
        # self.irotmq = 0
        # if self.minus_q:
        #   # Get the first symmetry:
        #   for k in range(self.nsym):
        #       # Position feels the symmetries with S (fortran S is transposed)
        #       # While q vector feels the symmetries with S^t (so no .T required for fortran matrix)
        #       new_q = self.s[:,:, k].dot(q)
        #       # Compare new_q with aq
        #       new_q = pbc_diff(new_q, [0, 0, 0])
        #       if  np.allclose(q, -new_q, atol=1e-3, rtol=0):
        #           #print("Found a symmetry that maps q to -q")
        #           #print("Symmetry:", self.s[:,:, k])
        #           #print("q:", np.round(q, 5))
        #           #print("new_q:", np.round(new_q, 5))
        #           self.irotmq = k + 1
        #           #print("IROTMQ from python:", self.irotmq)
        #           break
        #   if self.irotmq == 0:
        #       #print ("Error, the fortran code tells me there is S so that Sq = -q + G")
        #       #print ("But I did not find such a symmetry!")
        #       raise ValueError("Error in the symmetrization. See stdout")

    def setup_sg_symmetries(self, verbose=False):
        self.setup_little_cogroup([0, 0, 0], verbose=verbose)

    def get_star_q(self, q: ArrayLike, verbose=False):
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
        q = np.array(q, dtype=np.float64, order="F")
        self.setup_sg_symmetries(verbose=False)
        full_symmetries = np.copy(self.s)
        full_invs = np.copy(self.invs)
        # self.setup_little_cogroup(q, verbose=False)
        q = np.array(
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
        nq_new, sxq, isq, imq = espresso_symm.star_q(
            q,  # q point in cartesian coordinates/2pi
            self.at,  # lattice vectors (see __init__ for format)
            self.bg,  # rec. lattice vectors (see __init__ for format)
            self.nsymq,  # Number of symmetries in the small group of q
            full_symmetries,  # Array of ALL symmetries of the crystal
            full_invs,  # Index of inverse of s in self.s
            verbose,  # Verbosity flag
        )

        # print("----------Inside routine get_star_q")
        # print("I am getting nq_new:", nq_new)
        # print("I am getting sxq:", sxq.T[:nq_new])
        # print("I am getting isq:", isq)
        # print("I am getting imq:", imq)
        # do while (isq (isym) /= imq)
        #    isym = isym + 1
        # enddo
        # irotmq = np.where(isq == imq)[0][0]
        # print("I am getting irotmq:", irotmq)
        ##if imq > 0 and not np.allclose(sxq[:, irotmq], -q, atol=1e-3, rtol=0):
        ##    print("WARNING: sxq[imq] is not the inverse of -q")
        # print(f"sxq[{irotmq}]:", sxq[:, irotmq])
        # print("-q:", -q)
        # print("----------End of Inside routine get_star_q")

        # TODO: this implementation is quite confusing
        # TODO: it is only really necessary because the imq thing doesn't
        # TODO: work properly for nonsymmorphic space groups
        if imq != 0:
            total_star = np.zeros((nq_new, 3), dtype=np.float64)
        else:
            # If -q is not in the star we stick it in there
            total_star = np.zeros((2 * nq_new, 3), dtype=np.float64)

        total_star[:nq_new, :] = sxq[:, :nq_new].T

        if imq == 0:
            # Stick -q into the star
            total_star[nq_new:, :] = -sxq[:, :nq_new].transpose()

        return total_star

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
