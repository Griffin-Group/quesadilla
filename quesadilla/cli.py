"""
Quesadilla command line interface
"""

import argparse
import os
import sys
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import spglib
import yaml
from phonopy import Phonopy
from phonopy.cui.collect_cell_info import collect_cell_info
from phonopy.cui.phonopy_script import (
    _create_FORCE_SETS_from_settings,
    _init_phonopy,
    _produce_force_constants,
    _read_phonopy_settings,
)
from phonopy.cui.settings import PhonopySettings
from phonopy.interface.calculator import (
    add_arguments_of_calculators,
    calculator_info,
    get_default_cell_filename,
    get_default_displacement_distance,
    get_default_physical_units,
    write_crystal_structure,
    write_supercells_with_displacements,
)
from phonopy.interface.phonopy_yaml import PhonopyYaml
from phonopy.structure.atoms import PhonopyAtoms, atom_data
from phonopy.structure.cells import Primitive, get_primitive, guess_primitive_matrix

from quesadilla.dynmat import NondiagonalPhononCalculator
from quesadilla.supercells import SupercellGenerator


def get_parser():
    """Return ArgumentParser instance."""

    parser = argparse.ArgumentParser(
        description="Phonopy command-line-tool", allow_abbrev=False
    )

    add_arguments_of_calculators(parser, calculator_info)

    # GENERAL
    parser.add_argument(
        "--config",
        dest="conf_filename",
        metavar="FILE",
        default=None,
        help="Phonopy configuration file",
    )
    parser.add_argument(
        "--loglevel", dest="loglevel", type=int, default=None, help="Log level"
    )
    # TODO: implement magnetism
    parser.add_argument(
        "--magmom", nargs="+", dest="magmoms", default=None, help="Same as MAGMOM tag"
    )
    parser.add_argument(
        "-q",
        "--quiet",
        dest="quiet",
        action="store_true",
        default=None,
        help="Print out smallest information",
    )
    parser.add_argument(
        "--random-seed",
        dest="random_seed",
        type=int,
        default=None,
        help="Random seed by a 32 bit unsigned integer",
    )
    parser.add_argument(
        "--tolerance",
        dest="symmetry_tolerance",
        type=float,
        default=None,
        help="Symmetry tolerance to search",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=None,
        help="Detailed information is shown.",
    )
    parser.add_argument(
        "filename",
        nargs="*",
        help=(
            "Phonopy configure file. However if the file is recognized as "
            "phonopy.yaml like file, this file is read as phonopy.yaml like file."
        ),
    )
    ###################################
    # MODE 1 : pre-process
    ###################################
    # -----------------------------
    # Required in this mode
    parser.add_argument(
        "-d",
        "--displacement",
        dest="is_displacement",
        action="store_true",
        default=None,
        help="Create supercells with displacements",
    )
    parser.add_argument(
        "--dim",
        nargs="+",
        dest="supercell_dimension",
        default=None,
        help="Same behavior as DIM tag",
    )
    # ----------------------
    # Good default is guessed depending on the calculator
    parser.add_argument(
        "-c",
        "--cell",
        dest="cell_filename",
        metavar="FILE",
        default=None,
        help="Read unit cell",
    )
    # Control of displacement pattern
    parser.add_argument(
        "--amplitude",
        "--amin",
        dest="displacement_distance",
        type=float,
        default=None,
        help=(
            "Distance of displacements and also minimum distance of displacements "
            "in random displacements"
        ),
    )
    parser.add_argument(
        "--amax",
        dest="displacement_distance_max",
        type=float,
        default=None,
        help="Minimum distance of displacements in random displacements",
    )
    parser.add_argument(
        "--nodiag",
        dest="is_nodiag",
        action="store_true",
        default=None,
        help="Set displacements parallel to axes",
    )
    parser.add_argument(
        "--nosym",
        dest="is_nosym",
        action="store_true",
        default=None,
        help="Symmetry is not imposed.",
    )
    parser.add_argument(
        "--trigonal",
        dest="is_trigonal_displacements",
        action="store_true",
        default=None,
        help="Set displacements of all trigonal axes ",
    )
    parser.add_argument(
        "--pm",
        dest="is_plusminus_displacements",
        action="store_true",
        default=None,
        help="Set plus minus displacements",
    )
    # TODO: not sure how to deal with this one?
    parser.add_argument(
        "--wien2k-p1",
        dest="is_wien2k_p1",
        action="store_true",
        default=None,
        help="Assume Wien2k structs with displacements are P1",
    )

    # MODE 2: post-process
    # -----------------------------
    # MODE 2.1: generate force sets
    # -----------------------------
    parser.add_argument(
        "-f",
        "--force-sets",
        nargs="+",
        dest="create_force_sets",
        default=None,
        help="Create FORCE_SETS",
    )
    parser.add_argument(
        "--fz",
        "--force-sets-zero",
        nargs="+",
        dest="create_force_sets_zero",
        default=None,
        help=(
            "Create FORCE_SETS. disp.yaml in the current directory and "
            "vapsrun.xml's for VASP or case.scf(m) for Wien2k as arguments "
            "are required. The first argument is that of the perfect "
            "supercell to subtract residual forces"
        ),
    )
    # -----------------------------
    # MODE 2.2: generate force constants
    # -----------------------------
    # NOTE: this tag works differently from phonopy
    parser.add_argument(
        "--fc",
        "--force-constants",
        dest="create_force_constants",
        action="store_true",
        default=None,
        help=(
            "Create FORCE_CONSTANTS from vaspurn.xml. "
            "vasprun.xml has to be passed as argument."
        ),
    )
    # Force constant symmetry
    # Should be hijaked and applied only to last step (ND D(q) -> DSC D(q))
    # Individual steps should not use symmetry
    parser.add_argument(
        "--no-fc-symmetry",
        "--no-sym-fc",
        dest="fc_symmetry",
        action="store_false",
        default=None,
        help="Do not symmetrize force constants",
    )
    parser.add_argument(
        "--fc-spg-symmetry",
        dest="fc_spg_symmetry",
        action="store_true",
        default=None,
        help="Enforce space group symmetry to force constants",
    )
    parser.add_argument(
        "--fc-symmetry",
        "--sym-fc",
        dest="fc_symmetry",
        action="store_true",
        default=None,
        help="Symmetrize force constants",
    )

    return parser


def create_displacements(
    sc_gen: SupercellGenerator,
    settings: PhonopySettings,
    cell_info: dict,
    log_level: int,
):
    """
    Create displacements.
    """
    if settings.displacement_distance is None:
        displacement_distance = get_default_displacement_distance(settings.calculator)
    else:
        displacement_distance = settings.displacement_distance

    nd_phonons = []
    cell_info_ndsc = cell_info.copy()
    for i, sc_matrix in enumerate(sc_gen.sc_matrices):
        # NOTE: phonopy convention is transposed compared to PRB 92, 184301 (2015)
        cell_info_ndsc["supercell_matrix"] = sc_matrix.T
        phonon = _init_phonopy(
            settings, cell_info_ndsc, settings.symmetry_tolerance, log_level
        )
        phonon.generate_displacements(
            distance=displacement_distance,
            is_plusminus=settings.is_plusminus_displacement,
            is_diagonal=settings.is_diagonal_displacement,
            is_trigonal=settings.is_trigonal_displacement,
            number_of_snapshots=settings.random_displacements,
            random_seed=settings.random_seed,
            max_distance=settings.displacement_distance_max,
        )
        if log_level:
            print(
                f"Nondiagonal supercell {i+1} needs {len(phonon.supercells_with_displacements)} displacements."
            )
        nd_phonons.append(phonon)

    return nd_phonons


def write_ndsc_with_displacement(
    nd_phonons: list[Phonopy],
    calculator: str,
    confs: dict,
    optional_structure_info: tuple,
):
    """
    Write supercells with displacements as structure files

    Args:
        nd_phonons: list of Phonopy objects for each NDSC
        calculator: force calculator (e.g., "vasp")
        confs: dictionary of configuration as strings
        optional_structure_info: output of read_cell_info
    """
    units = get_default_physical_units(calculator)
    for i, phonon in enumerate(nd_phonons):
        os.makedirs(f"sc-{i+1:03d}", exist_ok=True)
        os.chdir(f"sc-{i+1:03d}")
        write_crystal_structure(
            get_default_cell_filename(calculator),
            phonon.primitive,
            interface_mode=calculator,
            optional_structure_info=optional_structure_info,
        )
        # Write supercells with displacements as structure files
        write_supercells_with_displacements(
            phonon.calculator,
            phonon.supercell,
            phonon.supercells_with_displacements,
            optional_structure_info=optional_structure_info,
            additional_info={"supercell_matrix": phonon.supercell_matrix},
        )
        # Write phonopy_disp.yaml
        _write_phonopy_disp_yaml(phonon, confs, units)
        os.chdir("..")


def _write_phonopy_disp_yaml(phonon: Phonopy, confs: dict, units: dict):
    """
    Write phonopy_disp.yaml file

    Args:
        phonon: Phonopy object
        confs: dict
        units: dict
    """
    yaml_settings = {
        "force_sets": False,
        "force_constants": False,
        "born_effective_charge": False,
        "dielectric_constant": False,
        "displacements": True,
    }
    confs["dim"] = " ".join(map(str, phonon.supercell_matrix.flatten().tolist()))
    phpy_yaml = PhonopyYaml(
        configuration=confs, physical_units=units, settings=yaml_settings
    )
    phpy_yaml.set_phonon_info(phonon)
    with open("phonopy_disp.yaml", "w") as w:
        w.write(str(phpy_yaml))


def _write_final_phonopy_yaml(phonon: Phonopy, confs: dict, units: dict):
    """
    Write the final phonopy.yaml file after force constants have been
    computed from all the NDSCs.
    """
    yaml_settings = {
        "force_sets": False,
        "force_constants": True,
        "born_effective_charge": False,
        "dielectric_constant": False,
        "displacements": False,
    }
    confs["dim"] = " ".join(map(str, phonon.supercell_matrix.flatten().tolist()))
    phpy_yaml = PhonopyYaml(
        configuration=confs, physical_units=units, settings=yaml_settings
    )
    phpy_yaml.set_phonon_info(phonon)
    with open("phonopy.yaml", "w") as w:
        w.write(str(phpy_yaml))


def initialize_quesadilla():
    """
    Initialize Quesadilla command line interface.

    Returns:
        PhonopySettings: Phonopy settings
        dict: Configuration settings
        str: Cell filename
        int: Log level
    """
    parser = get_parser()
    args = parser.parse_args()
    # Set log level
    log_level = 1
    if args.verbose:
        log_level = 2
    if args.quiet:
        log_level = 0
    if args.loglevel is not None:
        log_level = args.loglevel

    argparse_control = {
        "fc_symmetry": False,
        "is_nac": False,
        "load_phonopy_yaml": False,
    }
    settings, confs, cell_filename = _read_phonopy_settings(
        args, argparse_control, log_level
    )

    # These defaults are set internally inside multiple phonopy routines
    # But I want to have them here from the get-go
    if settings.calculator is None:
        settings.calculator = "vasp"
    if cell_filename is None:
        cell_filename = get_default_cell_filename(settings.calculator)
    if settings.symmetry_tolerance is None:
        settings.symmetry_tolerance = 1e-5

    # Ensure settings are correct
    if settings.create_displacements or settings.random_displacements:
        if settings.supercell_matrix is None:
            raise ValueError("Supercell matrix is required for displacements")
        if not np.all(
            np.diag(np.diag(settings.supercell_matrix)) == settings.supercell_matrix
        ):
            raise ValueError(
                "Supercell matrix (--dim or DIM) must be diagonal."
                "This is the size of your q-grid."
            )

    return settings, confs, cell_filename, log_level


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
    fname = f"P{cell_filename}"
    write_crystal_structure(
        fname,
        primitive,
        interface_mode=calculator,
        optional_structure_info=cell_info["optional_structure_info"],
    )
    print(f"Standard primitive cell written to {fname}")
    # Get cell_info for the new file
    # TODO: how to deal with magmoms here?
    cell_info = _get_cell_info(fname, calculator, magmoms)

    return primitive, cell_info


def create_force_sets(
    settings: PhonopySettings,
    log_level: int,
    forces_filename: str,
):
    """
    Create force sets from the displacements for each nondiagonal supercell.

    Args:
        settings (PhonopySettings): Phonopy settings
        log_level (int): Log level
        forces_filename (str): The name of the file to read the forces from. This is the single argument of --force-sets or --force-sets-zero. Then for each supercell in the directory sc-i, the files to read forces from are disp-i/forces_filename.
    """
    sc_gen = SupercellGenerator.from_toml("quesadilla.toml")
    for i in range(len(sc_gen.sc_matrices)):
        os.chdir(f"sc-{i+1:03d}")
        if log_level > 0:
            print(f"Working on supercell {i+1} in directory sc-{i+1:03d}")
        _set_filenames_for_forcesets(settings, log_level, forces_filename)
        _create_FORCE_SETS_from_settings(
            settings,
            "phonopy_disp.yaml",
            settings.symmetry_tolerance,
            log_level,
        )
        if log_level > 0:
            print("------------------------------------")
        os.chdir("..")


def _set_filenames_for_forcesets(
    settings: PhonopySettings, log_level: int, forces_filename: str
):
    """
    Sets the filenames for the forcesets. This is a helper function for create_forcesets. Sets settings.create_force_sets or settings.create_force_sets_zero to a list of files to read the forces from based on the `phonopy_disp.yaml` file.
    """
    with open("phonopy_disp.yaml", "r") as f:
        data = yaml.safe_load(f)
        num_disps = len(data["displacements"])
        if log_level > 0:
            print(f"Found {num_disps} displacements.")
    # Change filenames
    filenames = [
        os.path.join(f"disp-{i+1:03d}", forces_filename) for i in range(num_disps)
    ]

    if settings.create_force_sets_zero:
        filenames.insert(0, os.path.join("disp-000", forces_filename))
        settings.create_force_sets_zero = filenames
    else:
        settings.create_force_sets = filenames

    if log_level > 0:
        print("Files to extract forces from:")
        print(filenames)


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


def get_sc_gen(settings, cell_filename):
    """
    Prepares a supercell generate object from the settings and cell filename.

    This function dumps a standard primitive cell from the cell file and
    proceeds to create a SupercellGenerator object from it, and generate the supercells
    necessary to span the q-grid.
    """
    primitive, cell_info = get_phonopy_prim(
        cell_filename, calculator=settings.calculator
    )
    grid = np.diag(settings.supercell_matrix)
    sc_gen = SupercellGenerator(primitive, grid)
    # TODO: add arguments for minkowski and minimize
    sc_gen.generate_supercells(minimize_supercells=True)
    sc_gen.to_toml("quesadilla.toml")
    return cell_info, sc_gen
def create_force_constants(settings: Phonopy, confs: dict, log_level: int):
    """
    Create force constants from the force sets for each nondiagonal supercell.

    Reads the quesadilla.toml file to setup the SupercellGenerator,
    reads the force sets from each of the sc-i directories for each NDSC and
    computes the force constants for each NDSC, then uses a NondiagonalPhononCalculator
    to compute the final force constants for the full diagonal supercell.

    Args:
        settings (PhonopySettings): Phonopy settings
        log_level (int): Log level
        confs (dict): Configuration settings
    """
    sc_gen = SupercellGenerator.from_toml("quesadilla.toml")
    nd_phonons = _get_nd_phonons(sc_gen, settings, log_level)
    ndsc_calc = NondiagonalPhononCalculator(sc_gen, nd_phonons)
    ndsc_calc.run()
    units = get_default_physical_units(settings.calculator)
    _write_final_phonopy_yaml(ndsc_calc.phonons, confs, units)


def _get_nd_phonons(sc_gen, settings, log_level):
    """
    Get a list of Phonopy objects for each nondiagonal supercell.

    Args:
        sc_gen: SupercellGenerator object
        settings: PhonopySettings object
        log_level: int
    """

    # All nondiagonal supercells have the same primitive cell
    cell_filename = get_default_cell_filename(settings.calculator)
    cell_info = _get_cell_info(
        os.path.join("sc-001", cell_filename), settings.calculator
    )
    nd_phonons = []
    for i, sc_matrix in enumerate(sc_gen.sc_matrices):
        if log_level > 0:
            print(f"Reading force sets for supercell {i+1} in directory sc-{i+1:03d}")
        os.chdir(f"sc-{i+1:03d}")
        # TODO: need to deal with magmoms here?
        cell_info["supercell_matrix"] = sc_matrix.T
        # Initializes a phonopy object ready for force constant computation
        phonon = _init_phonopy(
            settings, cell_info, settings.symmetry_tolerance, log_level
        )
        # Read force sets and compute force constants in the phonon object
        _produce_force_constants(phonon, settings, None, None, log_level)

        nd_phonons.append(phonon)
        if log_level > 0:
            print("------------------------------------")
        units = get_default_physical_units(settings.calculator)
        _write_final_phonopy_yaml(phonon, {}, units)
        os.chdir("..")

    return nd_phonons


# Phonopy error class
class PhonopyError(Exception):
    pass


def main():
    """Main function of Quesadilla command line interface."""
    settings, confs, cell_filename, log_level = initialize_quesadilla()

    if settings.create_displacements or settings.random_displacements:
        cell_info, sc_gen = get_sc_gen(settings, cell_filename)
        nd_phonons = create_displacements(sc_gen, settings, cell_info, log_level)
        write_ndsc_with_displacement(
            nd_phonons, settings.calculator, confs, cell_info["optional_structure_info"]
        )
        sys.exit(0)

    if settings.create_force_sets or settings.create_force_sets_zero:
        forces_filename = (
            settings.create_force_sets[0]
            if settings.create_force_sets
            else settings.create_force_sets_zero[0]
        )
        create_force_sets(settings, log_level, forces_filename)
        sys.exit(0)

    if settings.create_force_constants:
        create_force_constants(settings, confs, log_level)
        sys.exit(0)


