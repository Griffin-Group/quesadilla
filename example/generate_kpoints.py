"""
Simple script to generate a KPOINTS file for a given structure and k-point spacing.

Usage:
    python generate_kpoints.py <path/to/POSCAR> <k_spacing>

Example:
    python generate_kpoints.py Si/sc-001/SPOSCAR 0.15
"""


import os
import sys

import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints


def get_KPOINTS(struct, kspace):
    """
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


sc_path = sys.argv[1]
root = os.path.dirname(sc_path)
k_spacing = float(sys.argv[2])

sc = Structure.from_file(sc_path)
print("I got structure")
print(sc)

KPOINTS = get_KPOINTS(sc, k_spacing)
kpoints_path = os.path.join(root, "KPOINTS")
KPOINTS.write_file(kpoints_path)

print(f"I wrote KPOINTS to {kpoints_path}")
print(KPOINTS)
