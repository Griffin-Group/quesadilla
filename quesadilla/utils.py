import os
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import spglib
from phonopy import Phonopy
from phonopy.cui.collect_cell_info import collect_cell_info
from phonopy.interface.calculator import (
    write_crystal_structure,
)
from phonopy.structure.atoms import PhonopyAtoms, atom_data
from phonopy.structure.cells import Primitive, get_primitive, guess_primitive_matrix
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from quesadilla.supercells import SupercellGenerator, make_positive_det

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
        matrices.append(make_positive_det(T))
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
    primitive_filename: str, grid: list, k_spacing: float = 0.15, verbose: bool = False
):
    """
    Generates the files for an NDSC phonon calculation.
    """

    root = os.path.dirname(primitive_filename)
    primitive = get_phonopy_prim(primitive_filename)
    sc_gen = SupercellGenerator(primitive, grid)
    sc_gen.generate_supercells()

    # pymatgen structure object
    new_prim_filename = os.path.join(root, "POSCAR_standard_primitive")
    prim = Structure.from_file(new_prim_filename)

    T_matrices, sc_size, comm_q = sc_gen.sc_matrices, sc_gen.sc_sizes, sc_gen.q_comm
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


def get_phonopy_prim(
    cell_filename: Path,
    calculator: str = "vasp",
    symprec: float = 1e-5,
    magmoms: Optional[np.ndarray] = None,
):
    # Get the cell info dict from the file
    cell_info = _get_cell_info(cell_filename, calculator, magmoms)
    # Find the standard primitive cell
    primitive = _find_standard_primitive(cell_info, symprec)
    # Write the standard primitive cell to a file
    fname = f"{cell_filename}_standard_primitive"
    write_crystal_structure(
        fname,
        primitive,
        interface_mode=calculator,
        optional_structure_info=cell_info["optional_structure_info"],
    )
    print(f"Standard primitive cell written to {fname}")

    return primitive


def _get_cell_info(
    cell_filename: Path, calculator: str, magmoms: Optional[np.ndarray] = None
):
    """
    Get the cell info from a file.

    Args:
        cell_filename (Path): The path to the cell file.
        calculator (str): The calculator to use.
        magmoms (np.ndarray): The magnetic moments.

    Returns:
        dict: The cell info dictionary.
    """
    # Get cell info
    cell_info = collect_cell_info(
        interface_mode=calculator,
        cell_filename=cell_filename,
        supercell_matrix=np.diag(np.ones(3, dtype=int)),
    )
    if "error_message" in cell_info:
        print("Phonopy returned this error while reading the cell:")
        print(cell_info["error_message"])
        raise PhonopyError(f"Error reading the cell file {cell_filename}")
    # Set magnetic moments
    if magmoms is not None:
        unitcell = cell_info["unitcell"]
        try:
            assert len(magmoms) in (len(unitcell), len(unitcell) * 3)
            unitcell.magnetic_moments = magmoms
        except AssertionError as e:
            raise PhonopyError(
                "Number of magnetic moments does not match the number of atoms or "
                "number of atoms times 3."
            ) from e

    return cell_info


def _find_standard_primitive(cell_info: dict, symprec: float) -> Primitive:
    """
    Finds the standard primitive cell of a structure. This should find a cell
    similar to what you get from running `phonopy --symmetry`

    Only supports non-magnetic systems. For magnetic systems it will just
    return the input cell as is.

    Args:
        cell_info (dict): The cell info dictionary.
        symprec (float): The symmetry precision.
    """
    phonon = Phonopy(
        cell_info["unitcell"],
        np.eye(3, dtype=int),
        primitive_matrix=cell_info["primitive_matrix"],
        symprec=symprec,
        calculator=cell_info["interface_mode"],
        log_level=0,
    )

    # Phonopy (and maybe spglib?) can't do this for magnetic systems
    if phonon.unitcell.magnetic_moments is not None:
        warnings.warn(
            (
                "Warning: phonopy cannot handle finding standard primitive cell for "
                "magnetic systems yet. I will pretend your input structure is "
                "in primitive in the right setting and proceed with the calculation. "
                "If not, this _may_ cause issues with quesadilla's nondiagonal "
                "supercell calculations. Proceed with caution."
            )
        )
        return phonon.primitive

    (bravais_lattice, bravais_pos, bravais_numbers) = spglib.refine_cell(
        phonon.primitive.totuple(), symprec
    )
    bravais_symbols = [atom_data[n][1] for n in bravais_numbers]
    bravais = PhonopyAtoms(
        symbols=bravais_symbols, scaled_positions=bravais_pos, cell=bravais_lattice
    )
    # Find the primitive cell
    trans_mat = guess_primitive_matrix(bravais, symprec=symprec)
    return get_primitive(bravais, trans_mat, symprec=symprec)


# Phonopy error class
class PhonopyError(Exception):
    pass
