import math
from fractions import Fraction
from typing import Tuple

import numpy as np
import pulp
import spglib
import tomli
import tomlkit
from numpy.typing import ArrayLike
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.structure.cells import Primitive, Supercell


class SupercellGenerator:
    def __init__(
        self,
        atoms: PhonopyAtoms,
        grid: ArrayLike,
        q_ibz: np.ndarray = None,
        sc_matrices: np.ndarray = None,
        sc_sizes: np.ndarray = None,
        q_comm: list = None,
    ):
        # Setup primitive structure and supercell
        grid = np.array(grid)
        self.supercell = Supercell(atoms, np.diag(grid))
        self.primitive = Primitive(self.supercell, np.diag(1 / grid))
        # Setup BZ data
        self.grid = grid
        if q_ibz is None:
            self.q_ibz = self.get_ibz()
        else:
            q_ibz = np.array(q_ibz)
            q_ibz = q_ibz[np.lexsort(q_ibz.T)]
            try:
                assert np.allclose(q_ibz, self.get_ibz())
            except AssertionError as e:
                raise ValueError(
                    "IBZ q-points passed do not match what we expect from the "
                    "structure's symmetry."
                ) from e

        # Setup supercell data
        self.sc_matrices = np.array(sc_matrices) if sc_matrices is not None else None
        self.sc_sizes = np.array(sc_sizes) if sc_sizes is not None else None
        self.q_comm = q_comm

    def get_ibz(
        self,
    ) -> np.ndarray:
        """
        Gets the q points in a structures IBZ. The output is sorted for consistency.

        Parameters:
        - primitive: phonopy.structure.cells.Primitive
            The input structure (primitive cell, assumed to be standardized).
        - grid: array-like
            The grid of q points (e.g., [4, 4, 4] for a 4x4x4 grid).

        Returns:
        - qpoints: numpy.ndarray
            The q points in the IBZ (fractional coordinates, sorted).
        """
        grid = np.array(self.grid)
        # Contruct the SPG-style cell tuple
        # https://spglib.readthedocs.io/en/stable/python-interface.html#crystal-structure-cell
        mapping = {
            element: i + 1 for i, element in enumerate(set(self.primitive.symbols))
        }
        zs = [mapping[element] for element in self.primitive.symbols]
        cell = (
            tuple(map(tuple, self.primitive.cell.tolist())),
            tuple(map(tuple, self.primitive.scaled_positions.tolist())),
            tuple(zs),
        )
        # Get irr q-points
        # TODO: cell = primitive.totuple()
        mapping, all_qpoints = spglib.get_ir_reciprocal_mesh(
            grid, cell, is_shift=np.zeros(3), symprec=1e-5
        )
        irr_qpoints = np.array([all_qpoints[idx] / grid for idx in np.unique(mapping)])

        # Sort by (z, y, x) for consistency
        return irr_qpoints[np.lexsort(irr_qpoints.T)]

    def _atoms_to_toml(self) -> tomlkit.table:
        """
        Convert the primitive structure to a TOML table.
        """
        if self.primitive is None:
            raise ValueError("Primitive structure must be set.")

        primitive_toml = tomlkit.table()
        primitive_toml.add("lattice", np.round(self.primitive.cell, 12).tolist())
        primitive_toml.add("species", self.primitive.symbols)
        primitive_toml.add(
            "frac_coords", np.round(self.primitive.scaled_positions, 12).tolist()
        )
        primitive_toml["lattice"].multiline(True)
        primitive_toml["frac_coords"].multiline(True)
        return primitive_toml

    def _bz_to_toml(self) -> tomlkit.table:
        """
        Convert the Brillouin zone data to a TOML table.
        """
        if self.grid is None or self.q_ibz is None:
            raise ValueError("Grid and irreducible q-points must be set.")
        bz_toml = tomlkit.table()
        bz_toml.add("grid", self.grid.tolist())
        bz_toml.add("irreducible_q", self.q_ibz.tolist())
        bz_toml["irreducible_q"].multiline(True)
        return bz_toml

    def _supercells_to_toml(self) -> tomlkit.aot:
        """
        Convert the supercells data to a TOML array of tables.
        """
        if self.sc_matrices is None or self.sc_sizes is None or self.q_comm is None:
            raise ValueError("Supercell data must be set.")

        supercells_toml = tomlkit.aot()
        for i, (T, sz, q) in enumerate(
            zip(self.sc_matrices, self.sc_sizes, self.q_comm)
        ):
            sc_table = tomlkit.table()
            sc_table.add("index", i + 1)
            sc_table.add("size", int(sz))
            sc_table.add("commensurate_q", q.tolist())
            sc_table["commensurate_q"].multiline(True)
            sc_table.add("matrix", T.tolist())
            sc_table["matrix"].multiline(True)
            supercells_toml.append(sc_table)
        return supercells_toml

    def to_toml(self, output_file):
        """
        Write the supercell data to a TOML file.

        Parameters:
        - output_file: str
            The output file path.
        """
        doc = tomlkit.document()
        doc.add("primitive", self._atoms_to_toml())
        doc.add("brillouin_zone", self._bz_to_toml())
        doc.add("supercells", self._supercells_to_toml())

        with open(output_file, "w") as f:
            f.write(doc.as_string())
        print(f"Supercells written to {output_file}")

    @classmethod
    def from_toml(cls, input_file: str):
        """
        Read the supercell data from a TOML file.

        Parameters:
        - input_file: str
            The input file path.
        """

        with open(input_file, "rb") as f:
            data = tomli.load(f)

        primitive = cls.atoms_from_toml(data["primitive"])
        bz = data["brillouin_zone"]
        grid = bz["grid"]
        irr_q = np.array(bz["irreducible_q"])

        if "supercells" in data:
            supercells = data["supercells"]
            T_matrices = np.array([sc["matrix"] for sc in supercells])
            sc_size = np.array([sc["size"] for sc in supercells])
            comm_q = [np.array(sc["commensurate_q"]) for sc in supercells]
        else:
            T_matrices = None
            sc_size = None
            comm_q = None
        return cls(primitive, grid, irr_q, T_matrices, sc_size, comm_q)

    @staticmethod
    def atoms_from_toml(data: dict) -> PhonopyAtoms:
        """Reconstructs a PhonopyAtoms object from TOML data."""
        return PhonopyAtoms(
            symbols=data["species"],
            scaled_positions=np.array(data["frac_coords"]),
            cell=np.array(data["lattice"]),
        )

    def generate_supercells(
        self, minkowski_reduce: bool = True, trim: bool = False
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate nondiagonal supercells commensurate with the IBZ.

        Parameters:
        - reduce: bool
            Whether to Minkowski reduce the supercells (default: True, recommended).
        - trim: bool
            Whether to trim the supercells using ILP (default: False, broken).
        """

        self.sc_matrices = self._get_ndsc_matrices()
        if minkowski_reduce:
            self.sc_matrices = np.array(
                [minkowski_reduce_sc(T, self.primitive.cell) for T in self.sc_matrices]
            )

        self.sc_sizes = np.array([np.linalg.det(T) for T in self.sc_matrices])
        # TODO: recompute all q-points commensurate with a supercell, not just one
        # self.q_comm = np.array([[q] for q in self.q_ibz])
        self.q_comm = [
            np.array([q for q in self.q_ibz if np.allclose(np.rint(T @ q), T @ q)])
            for T in self.sc_matrices
        ]

    def _get_ndsc_matrices(self) -> np.ndarray:
        qpoints_frac = convert_to_fraction_array(self.q_ibz)
        sc_matrices = np.zeros((self.q_ibz.shape[0], 3, 3), dtype=int)

        for i, Q in enumerate(qpoints_frac):
            sc_matrices[i, :] = find_nondiagonal(Q)

        return sc_matrices


# Utility functions for supercell generation
def convert_to_fraction_array(arr: np.ndarray) -> np.ndarray:
    """
    Takes a numpy array of floats and converts them to fractions.

    The output array has shape (N, 2, M) where N is the number of rows in the input array and M is the number of columns. The second dimension is used to store the numerator and denominator of the fraction, respectively.
    """
    result = np.empty((arr.shape[0], 2, arr.shape[1]), dtype=int)

    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            frac = Fraction(arr[i, j]).limit_denominator()
            result[i, 0, j] = frac.numerator
            result[i, 1, j] = frac.denominator

    return result


def find_integers(nums, g23, g12, g31, g123):
    """
    Compute integers for off-diagonal supercell matrix elements
    Called by find_nondiagonal()

    This function is copied from QEPlayground
    """
    # Compute p (it's a modulo equation)
    if g23 == 1:
        p = 0
    else:
        for i in range(1, g23):
            if (nums[1] + i * nums[2]) % g23 == 0:
                p = i
                break
        # Compute q
    g12_r = int(g12 / g123)
    g23_r = int(g23 / g123)
    g31_r = int(g31 / g123)
    if g12_r == 1:
        q = 0
    else:
        for i in range(1, g12_r):
            if (g23_r * nums[0] + i * g31_r * nums[1]) % g12_r == 0:
                q = i
                break
    # Compute r
    gg_r = int(g31 * g23 / g123)
    z = g23 * nums[0] / g12 + g31 * q * nums[1] / g12
    if gg_r == 1:
        r = 0
    else:
        for i in range(1, gg_r):
            if (z + i * nums[2]) % gg_r == 0:
                r = i
                break
    return p, q, r


def find_nondiagonal(Q: np.ndarray) -> np.ndarray:
    """
    Nondiagonal supercell, based on [Phys. Rev. B 92, 184301]
    This function is copied from QEPlayground

    Parameters:
    """
    # Take care of components already at Gamma
    Q[1, np.where(Q[0] == 0)] = 1
    # Shift the q-point into the positive quadrant of the reciprocal unit cell
    Q[0, np.where(Q[0] < 0)] += Q[1, np.where(Q[0] < 0)]
    # GCDs of Q[1] (in the logical order of the derivation)
    g23 = math.gcd(Q[1, 1], Q[1, 2])
    g12 = math.gcd(Q[1, 0], Q[1, 1])
    g31 = math.gcd(Q[1, 2], Q[1, 0])
    g123 = math.gcd(Q[1, 0], math.gcd(Q[1, 1], Q[1, 2]))
    # Integers needed to solve the supercell matrix equation
    p, q, r = find_integers(Q[0], g23, g12, g31, g123)
    # Matrix elements (in order of derivation) and supercell matrix
    S_33 = Q[1, 2]
    S_22 = Q[1, 1] / g23
    S_23 = p * Q[1, 2] / g23
    S_11 = g123 * Q[1, 0] / (g12 * g31)
    S_12 = q * g123 * Q[1, 1] / (g12 * g23)
    S_13 = r * g123 * Q[1, 2] / (g31 * g23)
    return np.array([[S_11, S_12, S_13], [0, S_22, S_23], [0, 0, S_33]])


# Functions for minimizing number of cells
def pick_smallest_supercells(
    commensurate: np.ndarray[bool], sc_sizes: np.ndarray[int], verbose: bool = False
) -> np.ndarray[bool]:
    """
    Solves the set cover problem with associated supercell sizes, selecting the
    smallest set of supercells that cover all q points.

    Parameters:
    -----------
    commensurate : numpy.ndarray
        A boolean NxN array where each row corresponds to a supercell and
        each column corresponds to a q points. An entry commensurate[i, j] is True if supercell `i` is commensurate with q-point `j`.

    sc_sizes : numpy.ndarray
        A 1D array of length N where each entry is the size of the corresponding
        supercell.

    verbose : bool
        Whether to print the solver output.

    Returns:
    --------
    selected_cells : list of int
        A list of indices of the selected supercells.

    total_size : int or float
        The total size of the selected supercells.

    Notes:
    ------
    This function uses integer linear programming (ILP) to ensure an optimal selection
    of supercells with the smallest total size while covering all q-points.
    The function requires the `pulp` library to solve the ILP problem.
    """
    N = len(sc_sizes)

    # Create a problem instance
    prob = pulp.LpProblem("PickSmallestSupercells", pulp.LpMinimize)

    # Create binary variables for each supercell
    x = [pulp.LpVariable(f"x_{i}", cat="Binary") for i in range(N)]

    # Objective function: Minimize the total supercell size
    prob += pulp.lpSum(sc_sizes[i] * x[i] for i in range(N))

    # Constraints: Ensure each q-point is covered by at least one bin
    for j in range(N):
        prob += pulp.lpSum(commensurate[i, j] * x[i] for i in range(N)) >= 1

    # Solve the problem
    prob.solve(pulp.PULP_CBC_CMD(msg=False))

    # Get the selected supercells and total size
    selected_cells = [int(pulp.value(x[i])) for i in range(N)]
    selected_cells = [i for i in range(N) if selected_cells[i] == 1]

    return selected_cells


# Minkowski reduction functions
def minkowski_reduce_sc(T, lattice):
    """
    Reduce a supercell matrix using Minkowski reduction.

    Args:
        T (numpy.ndarray): Supercell matrix.
        lattice (numpy.ndarray): Primitive lattice.

    Returns:
        numpy.ndarray: Minkowski-reduced supercell matrix with positive determinant.
    """
    ndsc_lattice = np.dot(T, lattice)
    ndsc_lattice = mink_reduce(ndsc_lattice)
    T = np.dot(
        ndsc_lattice,
        np.linalg.inv(lattice),
    )
    return make_positive_det(np.rint(T).astype(int))


def make_positive_det(matrix: np.ndarray) -> np.ndarray:
    """
    If the matrix has a negative determinant, this function flips the sign of the row with the most negative entries. Phonopy requires the supercell matrix to have a positive determinant. This doesn't change the q-point that the supercell is commensurate with.

    Args:
        matrix (numpy.ndarray): Input square matrix.

    Returns:
        numpy.ndarray: Adjusted matrix.
    """
    if np.linalg.det(matrix) < 0:
        # Find the row index with the most negative entries
        negative_sums = np.sum(np.minimum(matrix, 0), axis=1)
        row_to_flip = np.argmin(negative_sums)
        matrix[row_to_flip] *= -1

    return matrix


def mink_reduce(vecs: np.ndarray, tol: float = 1e-7, max_iter: int = 100) -> np.ndarray:
    """
    Perform Minkowski reduction on a set of 3 vectors in 3D space.

    Parameters
    ----------
    vecs : np.ndarray of shape (3, 3)
        Input array where each row represents a 3D vector.
    tol : float, optional
        Tolerance for floating-point comparisons, default is 1e-7.
    max_iter : int, optional
        Maximum number of iterations, default is 100.

    Returns
    -------
    np.ndarray
        Minkowski-reduced vectors.
    """
    if vecs.shape != (3, 3):
        raise ValueError("Input must have shape (3, 3).")

    i = 0
    while True:
        # Keep track of whether any reduction occurred
        changed = False

        for i in range(3):
            temp_vecs = vecs.copy()
            temp_vecs[i] = 0.0  # Temporarily zero out the i-th vector
            reduced_vecs = reduce_vectors(temp_vecs, tol)
            reduced_vecs[i] = vecs[i]  # Restore the i-th vector

            if not np.allclose(reduced_vecs, vecs):
                vecs = reduced_vecs
                changed = True
                break

        # Check combinations involving all three vectors
        if not changed:
            reduced_vecs = reduce_vectors(vecs, tol)
            if not np.allclose(reduced_vecs, vecs, atol=tol, rtol=0):
                vecs = reduced_vecs
                continue

        # Stop if no changes occurred in this iteration
        if not changed:
            break

        i += 1
        if i > max_iter:
            raise RuntimeError("Too many iterations in Minkowski reduction.")

    return vecs


def reduce_vectors(vecs: np.ndarray, tol: float) -> np.ndarray:
    """
    Reduce three 3D vectors by replacing the longest vector with a linear combination that is shorter.

    Parameters
    ----------
    vecs : np.ndarray of shape (3, 3)
        Input array where each row is a 3D vector.
    tol : float
        Tolerance for floating-point comparisons.

    Returns
    -------
    np.ndarray
        A new array with the reduced vectors if a reduction occurs, or the original array unchanged.
    """
    lengths_sq = np.sum(vecs**2, axis=1)
    longest_idx = np.argmax(lengths_sq)
    max_len_sq = lengths_sq[longest_idx]

    a, b, c = vecs
    candidates = [
        a + b - c,
        a - b + c,
        -a + b + c,
        a + b + c,
    ]

    for candidate in candidates:
        if np.dot(candidate, candidate) < max_len_sq * (1 - tol):
            new_vecs = vecs.copy()
            new_vecs[longest_idx] = candidate
            return new_vecs

    return vecs
