import itertools
import os

import numpy as np
import phonopy
from phonopy import Phonopy
from phonopy.harmonic.dynmat_to_fc import (
    DynmatToForceConstants,
)
from pymatgen.core.structure import Structure

import quesadilla.symmetries as symmetries
from quesadilla.supercells import struct_to_phonopy


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


def read_monserrat(path):
    """
    Reads kpoint_to_supercell.dat and associated supercell files.

    Args:
        path (str): Path to the directory containing kpoint_to_supercell.dat and supercell.<i>.dat files.

    Returns:
        tuple: A NumPy array of vectors (n x 3) and a list of associated matrices (each 3 x 3).
    """
    # File paths
    kpoint_file = os.path.join(path, "kpoint_to_supercell.dat")

    # Read kpoint_to_supercell.dat using NumPy
    kpoint_data = np.loadtxt(kpoint_file)

    # Extract vectors (first 3 columns) and indices (last column)
    q = kpoint_data[:, :3]
    indices = kpoint_data[:, 3].astype(int)

    # Read the associated supercell matrices
    matrices = []
    for i in indices:
        supercell_file = os.path.join(path, f"supercell.{i}.dat")
        T = np.loadtxt(supercell_file, dtype=int)  # Read the 3x3 matrix
        matrices.append(ensure_positive_det(T))
        # for qq in q:
        qq = q[i - 1]
        assert np.allclose(
            T @ qq, np.round(T @ qq)
        ), "Supercell {i} is not commensurate with q = {qq}"

    return matrices, q


def get_nd_phonopy(path, grid, T_sc, q_comm):
    prim = Structure.from_file(os.path.join(path, "POSCAR"))
    primitive, supercell = struct_to_phonopy(prim, grid)

    # Get fourier transformered FCs in IBZ and unfold to full BZ
    irr_fcq = get_fcq(path, T_sc, q_comm)
    full_fcq = {}
    symmetrizer = symmetries.Symmetrizer(prim)
    for q, D in irr_fcq.items():
        q = np.array(q)
        print(f"Applying q-star to q = {q.astype(float)}...")
        full_fcq |= symmetrizer.get_fcq_in_star(D, q)

    # Create the phonon object and symmetrize the FCs
    print("Creating the Phonopy object...")
    nd_phonon = Phonopy(primitive, supercell_matrix=np.diag(grid))
    print("Fourier transforming the force constants back to real space...")
    nd_phonon.force_constants = fcq_to_fcr(full_fcq, prim, grid)
    print("Symmetrizing the force constants...")
    nd_phonon.symmetrize_force_constants_by_space_group()
    nd_phonon.symmetrize_force_constants()

    return nd_phonon


def fcq_to_fcr(fcq, prim, grid):
    """
    Fourier transform the force constants on the full q-grid to real
    space.
    """
    primitive, supercell = struct_to_phonopy(prim, grid)
    all_q = np.array(list(fcq.keys()))
    print(f"Found {len(all_q)} q-points in the full BZ")
    dynmat = np.array([fcq_to_dynmat(prim, fcq[tuple(q)], q) for q in all_q])
    d2f = DynmatToForceConstants(
        primitive,
        supercell,
        is_full_fc=True,
    )
    d2f.commensurate_points = all_q
    d2f.dynamical_matrices = dynmat
    print("Running DynmatToForceConstants...")
    d2f.run()
    return d2f.force_constants


def get_fcq(path, T_sc, q_comm, average_gamma=False, average_non_gamma=False):
    """
    Gets the Fourier transformed force constantm atrices
    TODO: should take a list of force constants not path to the files!
    """
    prim = Structure.from_file(os.path.join(path, "POSCAR"))
    dynmat = {}
    n_sc = len(T_sc)

    for i in range(n_sc):
        phonon = phonopy.load(
            os.path.join(path, f"phonopy-{i + 1}.yaml"),
            is_symmetry=False,
            symmetrize_fc=False,
            is_compact_fc=False,
        )
        dyn_mat = phonon.dynamical_matrix

        q = q_comm[i]
        P = phonon.supercell_matrix
        # TODO: double check transpose crap again...
        assert np.allclose(P.T @ q, (P.T @ q).astype(int)), "q is not commensurate"
        dyn_mat.run(q)
        dynmat.setdefault(tuple(q), []).append(dyn_mat.dynamical_matrix)

    for q, dyn_mats in dynmat.items():
        print(f"Found {len(dyn_mats)} dynamical matrices at q = {q}")
        if q == (0, 0, 0):
            dynmat[q] = np.mean(dyn_mats, axis=0) if average_gamma else dyn_mats[0]
        elif len(dyn_mats) > 1:
            dynmat[q] = np.mean(dyn_mats, axis=0) if average_non_gamma else dyn_mats[0]
        else:
            dynmat[q] = dyn_mats[0]

    for q, D in dynmat.items():
        dynmat[q] = dynmat_to_fcq(prim, D, q)

    return dynmat


def dynmat_to_fcq(prim, D, q):
    N = len(prim)
    masses = np.array([s.species.elements[0].atomic_mass for s in prim])
    D_blocks = D.reshape(N, 3, N, 3).swapaxes(1, 2)
    for i, j in itertools.product(range(N), range(N)):
        r_i = prim[i].frac_coords
        r_j = prim[j].frac_coords
        D_blocks[i, j] *= (
            np.exp(1j * 2 * np.pi * np.dot(q, r_i))
            * np.exp(-1j * 2 * np.pi * np.dot(q, r_j))
            * np.sqrt(masses[i] * masses[j])
        )
    return D_blocks.swapaxes(1, 2).reshape(3 * N, 3 * N)


def fcq_to_dynmat(prim, fcq, q):
    N = len(prim)
    masses = np.array([s.species.elements[0].atomic_mass for s in prim])
    D_blocks = fcq.reshape(N, 3, N, 3).swapaxes(1, 2)
    for i, j in itertools.product(range(N), range(N)):
        r_i = prim[i].frac_coords
        r_j = prim[j].frac_coords
        D_blocks[i, j] *= (
            np.exp(-1j * 2 * np.pi * np.dot(q, r_i))
            * np.exp(1j * 2 * np.pi * np.dot(q, r_j))
            / np.sqrt(masses[i] * masses[j])
        )
    return D_blocks.swapaxes(1, 2).reshape(3 * N, 3 * N)
