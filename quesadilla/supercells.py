import math
from copy import deepcopy
from fractions import Fraction
from typing import Tuple

import numpy as np
import pulp
from numpy.typing import ArrayLike
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.structure.cells import Primitive, Supercell
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints
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


# Copied from QEPlayground
def find_integers(nums, g23, g12, g31, g123):
    """Compute integers for off-diagonal supercell matrix elements
    Called by find_nondiagonal()
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


# Copied from QEPlayground
def find_nondiagonal(Q):
    """Nondiagonal supercell, based on [Phys. Rev. B 92, 184301]"""
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


def get_qpoints(s, grid, shift=(0, 0, 0), ibz=True):
    sga = SpacegroupAnalyzer(s)
    if ibz:
        qpoints = np.array([q[0] for q in sga.get_ir_reciprocal_mesh(grid, shift)])
    else:
        qpoints, _ = sga.get_ir_reciprocal_mesh_map(grid, shift)
    return qpoints[~np.all(qpoints == [0, 0, 0], axis=1)]


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


def get_T_matrices(qpoints):
    qpoints_frac = convert_to_fraction_array(qpoints)
    T_matrices = np.zeros((qpoints.shape[0], 3, 3), dtype=int)

    # pymatgen uses different convention, need transpose
    for i, Q in enumerate(qpoints_frac):
        T_matrices[i, :] = find_nondiagonal(Q)

    # Eliminate identity if it's there
    # mask = ~np.all(T_matrices == np.eye(3), axis=(1,2))
    return T_matrices


# FIXME: How to do Niggli reduction without messing up
# the commensurate q-points? Here, the niggli-reduced cells
# actually no longer span the full IBZ...?
def get_supercells(prim, grid, shift=(0, 0, 0), reduce=True, ibz=True, trim=True):
    qpoints = get_qpoints(prim, grid, shift, ibz)
    T_matrices = get_T_matrices(qpoints)
    nq = qpoints.shape[0]
    commensurate = np.full((nq, nq), False, dtype=bool)
    for i, T in enumerate(T_matrices):
        commensurate[i, :] = [np.all(q @ T.T == np.round(q @ T.T)) for q in qpoints]
    sc_sizes = np.array([np.linalg.det(T) for T in T_matrices])

    # Trim to find the smallest number of supercells that span
    # the entire IBZ while minimizing the total supercell size
    if trim:
        selected_cells = pick_smallest_supercells(commensurate, sc_sizes)
    else:
        selected_cells = np.arange(len(sc_sizes))
    print(
        (
            f"Need {len(selected_cells)} supercells "
            f"with sizes {sc_sizes[selected_cells]}"
        )
    )
    supercells = []
    T_new = []
    for T in T_matrices[selected_cells]:
        sc = deepcopy(prim)
        sc.make_supercell(T.T)
        if reduce:
            sc_reduced = sc.get_reduced_structure(reduction_algo="niggli")
        else:
            sc_reduced = sc
        supercells.append(sc_reduced)

        # New transformation matrix
        S = prim.lattice.matrix.T
        D = sc_reduced.lattice.matrix.T
        T_new.append(np.round(np.matmul(np.linalg.inv(S), D), 10))

    # Get commensurate q points for the selected supercells
    q_comm = [qpoints[c] for c in commensurate[selected_cells, :]]
    return supercells, T_new, q_comm, T_matrices[selected_cells]


def minkowski_reduce(matrix):
    """
    TODO: implement this function
    """
    pass


def get_KPOINTS(struct, kspace):
    """
    TODO: switch to autoGR
    Prepares a KPOINTS object with a mesh of given spacing in 1/Å.
    Parameters:
    -----------
    struct : pymatgen.Structure
        The structure object.
    kspace : float
        The spacing of the k-point mesh in 1/Å.
    Returns:
    --------
    kpoints : pymatgen.io.vasp.inputs.Kpoints
        The KPOINTS object.
    """
    assert kspace > 0, "argument kspace must not be negative"
    b = np.array(struct.lattice.reciprocal_lattice.matrix)
    N = np.maximum(
        np.array([1, 1, 1]), np.round(np.sqrt(np.sum(b**2, axis=1)) / kspace)
    ).astype(np.int64)
    k_dict = {
        "nkpoints": 0,
        "generation_style": "Gamma",
        "kpoints": [[N[0], N[1], N[2]]],
        "usershift": [0, 0, 0],
        "comment": f"Mesh with spacing {kspace} 1/Å",
    }
    return Kpoints.from_dict(k_dict)


def ensure_positive_det(matrix):
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
