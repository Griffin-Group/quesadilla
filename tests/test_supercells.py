import os

import numpy as np
import pytest
from pymatgen.core.structure import Structure

import quesadilla
from quesadilla.utils import read_monserrat

TEST_DIR = os.path.dirname(os.path.abspath(__file__))


@pytest.mark.parametrize("material", ["Si", "CsCl"])
def test_phonon_band_comparison(material):
    root = os.path.join(TEST_DIR, "data", material)
    prim = Structure.from_file(os.path.join(root, "POSCAR"))

    # Generate the supercells using Quesadilla
    T_sc2, _, q_comm2 = quesadilla.supercells.get_supercells(prim, [4, 4, 4])

    # Read the stuff generated by Lloyd-Williams and Monserrat's code
    T_sc, q_comm = read_monserrat(os.path.join(root, "monserrat"))
    # Need to sort to have the same order as Quesadilla
    sort_idx = np.lexsort(q_comm.T)
    T_sc = T_sc[sort_idx]
    q_comm = q_comm[sort_idx]


    for i, (T1, T2) in enumerate(zip(T_sc, T_sc2)):
        assert np.allclose(q_comm[i], q_comm2[i])
        assert np.allclose(T1, T2)
