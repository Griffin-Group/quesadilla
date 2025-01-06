import os

import numpy as np
import phonopy
import pytest

import quesadilla.dynmat as dynmat

TEST_DIR = os.path.dirname(os.path.abspath(__file__))


# This unit test is just meant as a temporary solution for now
@pytest.mark.parametrize("material, threshold", [("Si", 5e-2), ("CsCl", 5e-2)])
def test_phonon_band_comparison(material, threshold):
    """
    Test to numerically compare phonon band structures between reference and computed data.
    Args:
        material (str): Material name for the test.
        threshold (float): Maximum allowed absolute difference in phonon frequencies.
    """
    # Load structures and phonon data
    root = f"{TEST_DIR}/data/{material}/"
    # prim = Structure.from_file(f"{root}/POSCAR")
    T_sc, q_comm = dynmat.read_monserrat(f"{root}/monserrat")
    nd_phonon = dynmat.get_nd_phonopy(f"{root}", [4, 4, 4], T_sc, q_comm)
    ref_phonon = phonopy.load(f"{root}/phonopy-diag.yaml")

    # Generate band structures
    ref_phonon.auto_band_structure(npoints=31)
    nd_phonon.auto_band_structure(npoints=31)

    # Extract frequencies
    ref_bands = ref_phonon.get_band_structure_dict()["frequencies"]
    nd_bands = nd_phonon.get_band_structure_dict()["frequencies"]

    # Compute absolute differences and assert they are within the threshold
    for b_ref, b_nd in zip(ref_bands, nd_bands):
        diff = np.abs(b_ref - b_nd)
        max_diff = np.max(diff)
        assert (
            max_diff < threshold
        ), f"Maximum difference {max_diff:.3e} exceeds threshold {threshold:.3e}"
