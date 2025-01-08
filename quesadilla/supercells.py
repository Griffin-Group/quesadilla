import math
from fractions import Fraction
from typing import Tuple

import numpy as np
import pulp
from numpy.typing import ArrayLike
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.structure.cells import Primitive, Supercell
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def struct_to_phonopy(prim: Structure, grid: ArrayLike) -> Tuple[Supercell, Primitive]:
    """
    Convert a pymatgen Structure to Phonopy's Supercell and Primitive objects.

    Parameters:
    -----------
    prim : pymatgen.Structure
        The primitive structure.

    grid : array-like
        The supercell transformation matrix
        TODO: define convention
    """
    atoms = PhonopyAtoms(
        symbols=[s.species_string for s in prim],
        masses=[s.species.elements[0].atomic_mass for s in prim],
        scaled_positions=prim.frac_coords,
        cell=prim.lattice.matrix,
    )
    grid = np.array(grid)
    supercell = Supercell(atoms, np.diag(grid))
    primitive = Primitive(supercell, np.diag(1 / grid))
    return primitive, supercell


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
    # FIXME: Added by Omar, needs a proper check why loop ends without finding r sometimes
    r = 0
    if gg_r == 1:
        r = 0
    else:
        for i in range(1, gg_r):
            if (z + i * nums[2]) % gg_r == 0:
                r = i
                break
    return p, q, r


def find_nondiagonal(Q):
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


def convert_to_fraction_array(arr):
    result = np.empty((arr.shape[0], 2, arr.shape[1]), dtype=int)

    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            frac = Fraction(arr[i, j]).limit_denominator()
            result[i, 0, j] = frac.numerator
            result[i, 1, j] = frac.denominator

    return result


def get_qpoints(
    structure: Structure,
    grid: ArrayLike,
    ibz: bool = True,
) -> np.ndarray:
    """
    Get the necessary q points for a given structure and grid.

    Parameters:
    - structure: pymatgen.Structure
        The input structure (primitive cell).
    - grid: array-like
        The grid of q points (e.g., [4, 4, 4] for a 4x4x4 grid).
    - ibz: bool
        Whether to return the irreducible q points. Always use True.

    Returns:
    - qpoints: numpy.ndarray
        The q points in the IBZ (fractional coordinates).
    """
    sga = SpacegroupAnalyzer(structure)
    if ibz:
        qpoints = np.array([q[0] for q in sga.get_ir_reciprocal_mesh(grid)])
    else:
        qpoints, _ = sga.get_ir_reciprocal_mesh_map(grid)
    return qpoints[np.lexsort(qpoints.T)]


def get_T_matrices(qpoints: np.ndarray) -> np.ndarray:
    qpoints_frac = convert_to_fraction_array(qpoints)
    T_matrices = np.zeros((qpoints.shape[0], 3, 3), dtype=int)

    # pymatgen uses different convention, need transpose
    for i, Q in enumerate(qpoints_frac):
        T_matrices[i, :] = find_nondiagonal(Q)

    return T_matrices


def get_supercells(
    prim: Structure, irr_q, reduce: bool = True, trim: bool = False
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generate nondiagonal supercells commensurate with the IBZ.

    Parameters:
    - prim: pymatgen.Structure
        The input structure (primitive cell).
    - grid: array-like
        The grid of q points (e.g., [4, 4, 4] for a 4x4x4 grid).
    - reduce: bool
        Whether to Minkowski reduce the supercells (default: True, recommended).
    - trim: bool
        Whether to trim the supercells using ILP (default: False, broken).

    Returns:
    - T_mat: list of numpy.ndarray
        The transformation matrices.
    - q_comm: list of numpy.ndarray
        q_comm[i] is the q-vector (frac coords) commensurate with T_mat[i]
    """

    T_matrices = get_T_matrices(irr_q)
    # TODO: recompute all q-points commensurate with a supercell, not just one
    q_comm = irr_q
    for i, T in enumerate(T_matrices):
        ndsc_lattice = np.dot(T, prim.lattice.matrix)
        ndsc_lattice = minkowski_reduce(ndsc_lattice)
        T = np.dot(ndsc_lattice, prim.lattice.reciprocal_lattice.matrix.T / 2 / np.pi)
        T_matrices[i] = ensure_positive_det(np.rint(T).astype(int))
    T_matrices = np.array(T_matrices)
    sc_size = np.array([np.linalg.det(T) for T in T_matrices])

    return T_matrices, sc_size, q_comm


def pick_smallest_supercells(commensurate, sc_sizes, verbose=False):
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


# FIXME: How to do Niggli reduction without messing up
# the commensurate q-points? Here, the niggli-reduced cells
# actually no longer span the full IBZ...?


def ensure_positive_det(matrix: np.ndarray) -> np.ndarray:
    """
    If the matrix has a negative determinant, this function flips the sign of the row with the most negative entries. Phonopy requires the supercell matrix to have a positive determinant.

    Args:
        matrix (numpy.ndarray): Input square matrix.

    Returns:
        numpy.ndarray: Adjusted matrix.
    """
    if np.linalg.det(matrix) < 0:
        # Calculate the sum of negative entries in each row
        negative_sums = np.sum(np.minimum(matrix, 0), axis=1)

        # Find the row index with the most negative entries
        row_to_flip = np.argmin(negative_sums)

        # Flip the sign of the selected row
        matrix[row_to_flip] *= -1

    return matrix


# TODO: these two functions are written very similarly to
# the FORTRAN code so it's they're very difficult to understand.
# Needs a complete rewrite in a pythonic and understandable way.
def minkowski_reduce(vecs: np.ndarray) -> np.ndarray:
    """
    Given a set of 3 vectors in 3D space (rows of `vecs`) that might not be
    Minkowski-reduced, iteratively attempt to Minkowski-reduce them until no
    further changes occur.

    Minkowski reduction (in 3D) implies:
      - The first vector is the shortest non-zero vector in the lattice.
      - Each subsequent vector is the shortest possible that still ensures
        linear independence with previously chosen vectors.

    This function modifies `vecs` in place until it is Minkowski-reduced.

    Parameters
    ----------
    vecs : np.ndarray of shape (3, 3)
        Each row is a 3D vector. This array is modified in place during the
        iterative process.

    Returns
    -------
    np.ndarray
        The Minkowski-reduced vectors (same reference as `vecs`).
    """
    # Sanity check
    if vecs.shape != (3, 3):
        raise ValueError("Input array 'vecs' must have shape (3, 3).")

    # We'll keep iterating until no more changes happen
    while True:
        # Save a copy of the current vectors so we can restore each row after testing
        # vecs_snapshot = vecs.copy()

        # 1) Check linear combinations involving two vectors
        #    We do this by zeroing out each row in turn, calling reduce_vec, then restoring.
        changed_in_this_pass = False
        for i in range(3):
            saved_row = vecs[i].copy()  # backup
            vecs[i] = 0.0  # zero out
            changed = reduce_vec(vecs)  # see if a reduction is triggered
            vecs[i] = saved_row  # restore the row

            if changed:
                changed_in_this_pass = True
                break  # break out of this for-loop to restart the outer loop

        if changed_in_this_pass:
            continue  # restart from scratch

        # 2) Check linear combinations involving all three vectors directly
        changed = reduce_vec(vecs)
        if changed:
            continue  # if it changed, restart

        # 3) If we got here, no changes occurred in either step -> done
        break

    return vecs


def reduce_vec(vecs: np.ndarray, tol: float = 1e-7) -> bool:
    """
    Attempt to reduce three 3D vectors (rows of vecs) by checking specific
    linear combinations. If any combination is shorter than the current longest
    vector, replace that longest vector with the shorter one.

    Parameters
    ----------
    vecs : np.ndarray of shape (3, 3)
        Each row is a 3D vector: [v1, v2, v3].
        This array is modified in place if a reduction occurs.
    tol : float, optional
        Relative tolerance for floating-point comparisons.

    Returns
    -------
    bool
        True if a replacement (and thus a reduction) happened.
        False otherwise.
    """
    # Compute the squared lengths of each of the three rows (vectors)
    lengths_sq = np.sum(vecs**2, axis=1)
    # Identify the index of the longest vector
    longest_idx = np.argmax(lengths_sq)
    max_len_sq = lengths_sq[longest_idx]

    # Construct the 4 new candidate vectors
    # Using Python’s array slicing to keep it readable:
    a, b, c = vecs
    new_vectors = np.array(
        [
            a + b - c,
            a - b + c,
            -a + b + c,
            a + b + c,
        ]
    )

    # Check whether any candidate is strictly shorter (within tolerance)
    for new_v in new_vectors:
        new_len_sq = np.dot(new_v, new_v)
        # Compare with a relative tolerance
        if new_len_sq < max_len_sq - tol * max_len_sq:
            # Replace the longest vector with the new shorter one
            vecs[longest_idx] = new_v
            return True
    return False
