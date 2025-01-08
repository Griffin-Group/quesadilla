import os

import numpy as np
import tomli
import tomlkit
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from quesadilla.supercells import ensure_positive_det, get_qpoints, get_supercells

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


def write_monserrat(prim, grid, path):
    sga = SpacegroupAnalyzer(prim)
    q_irr = sga.get_ir_reciprocal_mesh(grid)
    q_irr = np.array([q[0] for q in q_irr])
    # Write lattice matrix to path/prim.dat
    np.savetxt(f"{path}/prim.dat", prim.lattice.matrix, fmt="%.10f")
    # Write irreducible q-points to path/ibz.dat
    np.savetxt(f"{path}/ibz.dat", q_irr, fmt="%.10f")
    # Write grid to path/grid.dat
    np.savetxt(f"{path}/grid.dat", np.array(grid), fmt="%d", newline=" ")


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
        # matrices.append(T)
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


def generate_files(prime_filename: str, grid: list, k_spacing: float = 0.15):
    """
    Generates the files for an NDSC phonon calculation.
    """
    root = os.path.dirname(prime_filename)
    prim_original = Structure.from_file(prime_filename)
    # Standardize the primitive structure
    prim = standardize_prim(prim_original, root)
    irr_q = get_qpoints(prim, grid)
    T_matrices, sc_size, comm_q = get_supercells(prim, irr_q)
    for i, (T, sz, q) in enumerate(zip(T_matrices, sc_size, comm_q)):
        # Print the supercell information
        print(
            f"Supercell {i+1} with size {sz} is commensurate with q = {np.round(q, 3)}"
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
    # write_quesadilla_yaml(T_matrices, sc_size, comm_q, root)
    write_quesadilla_toml(prim, grid, irr_q, T_matrices, sc_size, comm_q, root)

    return prim, T_matrices, sc_size, comm_q


def standardize_prim(prim_original, root) -> Structure:
    prim = SpacegroupAnalyzer(prim_original).get_primitive_standard_structure()
    prim.translate_sites(np.arange(len(prim)), [0.0, 0.0, 0.0], to_unit_cell=True)
    if prim_original != prim:
        prime_filename = os.path.join(root, "POSCAR_symm")
        prim.to(filename=prime_filename, fmt="poscar")
        print("NOTE: The primitive structure has been standardized and dumped to")
        print(f"      {prime_filename}")
    return prim


def structure_to_toml(structure: Structure) -> tomlkit.table:
    primitive = tomlkit.table()
    primitive.add("lattice", np.round(structure.lattice.matrix, 10).tolist())
    primitive.add("species", [site.species_string for site in structure])
    primitive.add("frac_coords", np.round(structure.frac_coords, 10).tolist())
    primitive["lattice"].multiline(True)
    primitive["frac_coords"].multiline(True)

    return primitive


def bz_to_toml(grid: list, irr_q: np.ndarray) -> tomlkit.table:
    bz = tomlkit.table()
    bz.add("grid", grid)
    bz.add("irreducible_q", irr_q.tolist())
    bz["irreducible_q"].multiline(True)
    return bz


def structure_from_toml(data: dict) -> Structure:
    """Reconstructs a pymatgen Structure from TOML data."""
    lattice = Lattice(np.array(data["lattice"]))
    species = data["species"]
    frac_coords = np.array(data["frac_coords"])
    return Structure(lattice, species, frac_coords)


def write_quesadilla_toml(prim, grid, irr_q, T_matrices, sc_size, comm_q, output_path):
    doc = tomlkit.document()
    doc.add("primitive", structure_to_toml(prim))
    doc.add("brillouin_zone", bz_to_toml(grid, irr_q))

    supercells = tomlkit.aot()  # Array of Tables
    output_file = os.path.join(output_path, "quesadilla.toml")
    for i, (T, sz, q) in enumerate(zip(T_matrices, sc_size, comm_q)):
        sc_table = tomlkit.table()
        sc_table.add("index", i + 1)
        sc_table.add("size", int(sz))
        sc_table.add("commensurate_q", q.tolist())
        sc_table.add("matrix", T.tolist())
        sc_table["matrix"].multiline(True)
        supercells.append(sc_table)

    doc.add("supercells", supercells)
    with open(output_file, "w") as f:
        f.write(doc.as_string())
    print(f"Supercells written to {output_file}")


def read_quesadilla_toml(input_file: str):
    with open(input_file, "rb") as f:
        data = tomli.load(f)

    # Read primitive structure
    prim = structure_from_toml(data["primitive"])
    # Read BZ Data
    bz = data["brillouin_zone"]
    grid = bz["grid"]
    irr_q = np.array(bz["irreducible_q"])
    # Read supercells
    supercells = data["supercells"]
    T_matrices = np.array([sc["matrix"] for sc in supercells])
    sc_size = np.array([sc["size"] for sc in supercells])
    comm_q = np.array([sc["commensurate_q"] for sc in supercells])
    return prim, grid, irr_q, T_matrices, sc_size, comm_q
