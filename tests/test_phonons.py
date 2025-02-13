import os

import numpy as np
import phonopy
import pytest

from quesadilla.dynmat import NondiagonalPhononCalculator

TEST_DIR = os.path.dirname(os.path.abspath(__file__))


# This unit test is just meant as a temporary solution for now
@pytest.mark.parametrize(
    "material, threshold", [("Si", 0.005), ("CsCl", 0.04), ("hcp_He_28GPa", 0.03)]
)
def test_phonon_band_comparison(material, threshold):
    """
    Test to numerically compare phonon band structures between reference and computed data.
    Args:
        material (str): Material name for the test.
        threshold (float): Maximum allowed absolute difference in phonon frequencies.
    """
    # Load structures and phonon data
    root = os.path.join(TEST_DIR, "data", material)
    ndsc_calc = NondiagonalPhononCalculator.from_toml(f"{root}/quesadilla.toml")
    ndsc_calc.run()
    nd_phonon = ndsc_calc.phonons
    ref_phonon = phonopy.load(os.path.join(root, "phonopy-diag.yaml"))

    # Generate band structures
    ref_phonon.auto_band_structure(npoints=24)
    nd_phonon.auto_band_structure(npoints=24)

    # Extract frequencies
    ref_bands = ref_phonon.get_band_structure_dict()["frequencies"]
    nd_bands = nd_phonon.get_band_structure_dict()["frequencies"]

    # Compute max relative error and assert they are within the threshold
    # FIXME: This is nasty...
    for b_ref, b_nd in zip(ref_bands, nd_bands):
        err = np.abs(b_ref - b_nd) / b_ref
        # Exclude outliers
        # (some isolated points can have weird errors due to how band indices are assigned)
        for e in err.T:
            idx = np.abs(e - np.mean(e)) < 3 * np.std(e)
            max_err = np.max(e[idx])
            assert (
                max_err < threshold
            ), f"Maximum relative error {max_err:.3e} exceeds threshold {threshold:.3e}"
