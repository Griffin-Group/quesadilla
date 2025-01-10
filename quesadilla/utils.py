import os

import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from quesadilla.supercells import SupercellGenerator, ensure_positive_det

GET_YAML = """
# Read force sets from text file
FORCE_SETS = READ
# Write force_constants.hdf5
FORCE_CONSTANTS = WRITE
FC_FORMAT = HDF5
# Do not ymmetrize force constants
FC_SYMMETRY = .FALSE.
#SYMMETRY = .FALSE.
# Add force constants to .yaml file
INCLUDE_FC = .TRUE.
# Use primitive cell as is
PRIMITIVE_AXES = 1 0 0 0 1 0 0 0 1
# Supercell matrix
DIM = {}
"""

MAKE_DISP = """
# Displace atoms
CREATE_DISPLACEMENTS = .TRUE.
# Use primitive cell as is
PRIMITIVE_AXES = 1 0 0 0 1 0 0 0 1
# Supercell matrix
DIM = {}
"""


def write_lwm(prim, grid, path):
    """
    Writes the lattice matrix, irreducible q-points, and grid to files
    suitable for input to Lloyd-Williams and Monserrat's code.

    Args:
        prim (pymatgen.Structure): The primitive structure.
        grid (list): The grid of q-points.
        path (str): The directory to write the files to.
    """
    sga = SpacegroupAnalyzer(prim)
    q_irr = sga.get_ir_reciprocal_mesh(grid)
    q_irr = np.array([q[0] for q in q_irr])
    # Write lattice matrix to path/prim.dat
    np.savetxt(f"{path}/prim.dat", prim.lattice.matrix, fmt="%.10f")
    # Write irreducible q-points to path/ibz.dat
    np.savetxt(f"{path}/ibz.dat", q_irr, fmt="%.10f")
    # Write grid to path/grid.dat
    np.savetxt(f"{path}/grid.dat", np.array(grid), fmt="%d", newline=" ")


def read_lwm(path):
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

    return np.array(matrices), np.array(q)


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


def generate_files(
    prime_filename: str, grid: list, k_spacing: float = 0.15, verbose: bool = False
):
    """
    Generates the files for an NDSC phonon calculation.
    """

    root = os.path.dirname(prime_filename)
    structure = Structure.from_file(prime_filename)
    sc_gen = SupercellGenerator(structure, grid)
    sc_gen.generate_supercells()

    T_matrices, sc_size, comm_q = sc_gen.sc_matrices, sc_gen.sc_sizes, sc_gen.q_comm
    prim = sc_gen._pmg_prim
    for i, (T, sz, q) in enumerate(zip(T_matrices, sc_size, comm_q)):
        # Print the supercell information
        if verbose:
            print(
                (
                    f"Supercell {i+1} with size {sz} is commensurate with q = "
                    f"{np.round(q, 3)}"
                )
            )
            print("Supercell matrix:")
            print(T)
        # Generate the supercell
        sc_path = os.path.join(root, f"sc-{i+1}")
        os.makedirs(sc_path, exist_ok=True)
        sc = prim.copy()
        sc.make_supercell(T)
        prim.to(filename=os.path.join(sc_path, "POSCAR"), fmt="poscar")
        sc.to(filename=os.path.join(sc_path, "quesadilla_ndsc.vasp"), fmt="poscar")
        get_KPOINTS(sc, k_spacing).write_file(os.path.join(sc_path, "KPOINTS"))

        # DIM string uses supercell matrix, flattened. Phonopy uses transpose
        T_str = " ".join([f"{int(x):d}" for x in T.T.flatten()])
        with open(os.path.join(sc_path, "make_disp.conf"), "w") as f:
            f.write(MAKE_DISP.format(T_str))
        with open(os.path.join(sc_path, "get_yaml.conf"), "w") as f:
            f.write(GET_YAML.format(T_str))
    sc_gen.to_toml(os.path.join(root, "quesadilla.toml"))
