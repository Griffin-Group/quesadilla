import itertools
import os
from pathlib import Path

import numpy as np
import phonopy
from numpy.typing import ArrayLike
from phonopy import Phonopy
from phonopy.harmonic.dynmat_to_fc import (
    DynmatToForceConstants,
)

import quesadilla.symmetries as symmetries
from quesadilla.supercells import SupercellGenerator


class NondiagonalPhononCalculator:
    def __init__(self, sc_gen: SupercellGenerator, nd_phonons: list[Phonopy]):
        self.sc_gen = sc_gen
        self.nd_phonons = nd_phonons

    @classmethod
    def from_toml(cls, toml_file: Path):
        """
        Create a NondiagonalPhononCalculator object from a toml file.

        Generally the toml file is called quesadilla.toml and will be
        in the same directory as the sc-i directories, each of which
        should contain a phonopy.yaml file THAT MUST INCLUDE THE
        FORCE CONSTANTS computed from this supercell.

        Args:
            toml_file: Path to the toml file

        Returns:
            A NondiagonalPhononCalculator object
        """
        # TODO: move to CLI (maybe?)
        sc_gen = SupercellGenerator.from_toml(toml_file)
        root = os.path.dirname(toml_file)
        nd_phonons = cls._parse_ndsc_phonons(root, sc_gen)

        return cls(sc_gen, nd_phonons)

    def run(self, verbose: bool = False):
        # Get fourier transformered FCs in IBZ
        irr_dynmat = self._get_dynmat_in_ibz()
        irr_fcq = {q: self._dynmat_to_fcq(dyn, q) for q, dyn in irr_dynmat.items()}

        ## Unfold the FCs to the full star of each q point
        full_fcq = {}
        symmetrizer = symmetries.Symmetrizer(self.sc_gen.primitive)
        for q, D in irr_fcq.items():
            full_fcq |= symmetrizer.get_fcq_in_star(D, q, verbose=verbose)

        self.phonons = self._get_phonopy_from_fcq(full_fcq)

    @staticmethod
    def _parse_ndsc_phonons(root: Path, sc_gen: SupercellGenerator) -> list[Phonopy]:
        """
        Reads the phonopy.yaml files from the supercells generated by Quesadilla

        Args:
            root: Path to the directory containing the supercells
            sc_gen: SupercellGenerator object
        """
        # TODO: move to CLI (maybe?)
        n_sc = len(
            [d for d in os.listdir(root) if d.startswith("sc-") and d[3:].isdigit()]
        )
        assert n_sc == len(sc_gen.sc_matrices), "Number of supercells does not match"
        phonons = [
            phonopy.load(
                os.path.join(root, f"sc-{i+1:03d}", "phonopy.yaml"),
                is_symmetry=False,
                symmetrize_fc=False,
            )
            for i in range(n_sc)
        ]
        # Now assert that the supercell matrices are the same
        for i, (p, T) in enumerate(zip(phonons, sc_gen.sc_matrices)):
            assert np.all(
                p.supercell_matrix.T == T
            ), f"{i+1}: Supercell matrices do not match, {p.supercell_matrix.T} vs {T}"

        return phonons

    def _get_dynmat_in_ibz(
        self,
        average_gamma: bool = False,
        average_non_gamma: bool = False,
    ) -> dict[tuple[float, float, float], np.ndarray]:
        """
        Gets the Fourier transformed force constant matrices
        """
        dynmats = {}
        for i, p in enumerate(self.nd_phonons):
            for q in self.sc_gen.q_comm[i]:
                T = p.supercell_matrix
                assert np.allclose(
                    T.T @ q, (T.T @ q).astype(int)
                ), "q is not commensurate?"
                p.dynamical_matrix.run(q)
                dynmats.setdefault(tuple(q), []).append(
                    p.dynamical_matrix.dynamical_matrix
                )

        for q, dyn_mats in dynmats.items():
            print(f"Found {len(dyn_mats)} dynamical matrices at q = {np.round(q, 5)}")
            if q == (0, 0, 0):
                dynmats[q] = np.mean(dyn_mats, axis=0) if average_gamma else dyn_mats[0]
            elif len(dyn_mats) > 1:
                dynmats[q] = (
                    np.mean(dyn_mats, axis=0) if average_non_gamma else dyn_mats[0]
                )
            else:
                dynmats[q] = dyn_mats[0]

        return dynmats

    def _get_phonopy_from_fcq(
        self, full_fcq: dict[tuple[float, float, float], np.ndarray]
    ) -> Phonopy:
        """
        Create a Phonopy object from the full force constant matrix.

        Args:
            sc_gen: SupercellGenerator object
            full_fcq: The full force constant matrix
        """
        print("Creating the Phonopy object...")
        nd_phonon = Phonopy(
            self.sc_gen.primitive, supercell_matrix=np.diag(self.sc_gen.grid)
        )
        print("Fourier transforming the force constants back to real space...")
        nd_phonon.force_constants = self._fcq_to_fcr(
            full_fcq,
        )
        print("Applying acoustic sum rules...")
        nd_phonon.symmetrize_force_constants()
        return nd_phonon

    def _fcq_to_fcr(
        self, fcq: dict[tuple[float, float, float], np.ndarray]
    ) -> np.ndarray:
        """
        Fourier transform the force constants on the full q-grid to real
        space.
        """
        all_q = np.array(list(fcq.keys()))
        print(f"Found {len(all_q)} q-points in the full BZ")
        dynmat = np.array([self._fcq_to_dynmat(fcq[tuple(q)], q) for q in all_q])
        d2f = DynmatToForceConstants(
            self.sc_gen.primitive,
            self.sc_gen.supercell,
            is_full_fc=True,
        )
        d2f.commensurate_points = all_q
        d2f.dynamical_matrices = dynmat
        print("Running DynmatToForceConstants...")
        d2f.run()
        return d2f.force_constants

    def _dynmat_to_fcq(self, D: np.ndarray, q: ArrayLike) -> np.ndarray:
        """
        Converts the dynamical matrix as defined by phonopy into
        the Fourier transformer force constant matrix as defined
        by Quantum ESPRESSO and read by Quesadilla's symmetry routines.

        Args:
            prim: The primitive structure
            D: The dynamical matrix (3*nat x 3*nat)
            q: The q-point in FRAC coords

        Returns:
            The Fourier transformed force constant matrix (3*nat x 3*nat)
        """
        N = len(self.sc_gen.primitive)
        D_blocks = D.reshape(N, 3, N, 3).swapaxes(1, 2)
        masses = self.sc_gen.primitive.masses
        for i, j in itertools.product(range(N), range(N)):
            # Get fractional coordinates
            r_i = self.sc_gen.primitive.scaled_positions[i]
            r_j = self.sc_gen.primitive.scaled_positions[j]
            D_blocks[i, j] *= (
                np.exp(1j * 2 * np.pi * np.dot(q, r_i))
                * np.exp(-1j * 2 * np.pi * np.dot(q, r_j))
                * np.sqrt(masses[i] * masses[j])
            )
        return D_blocks.swapaxes(1, 2).reshape(3 * N, 3 * N)

    def _fcq_to_dynmat(self, fcq: np.ndarray, q: np.ndarray):
        """
        Converts the Fourier transformed force constant matrix as defined
        by Quantum ESPRESSO and produced by Quesadilla's symmetry routines into
        the dynamical matrix as defined by phonopy.

        Args:
            prim: The primitive structure
            fcq: The Fourier transformed force constant matrix (3*nat x 3*nat)
            q: The q-point in FRAC coords

        Returns:
            The dynamical matrix (3*nat x 3*nat)
        """
        N = len(self.sc_gen.primitive)
        masses = self.sc_gen.primitive.masses
        D_blocks = fcq.reshape(N, 3, N, 3).swapaxes(1, 2)
        for i, j in itertools.product(range(N), range(N)):
            # Get fractional coordinates
            r_i = self.sc_gen.primitive.scaled_positions[i]
            r_j = self.sc_gen.primitive.scaled_positions[j]
            D_blocks[i, j] *= (
                np.exp(-1j * 2 * np.pi * np.dot(q, r_i))
                * np.exp(1j * 2 * np.pi * np.dot(q, r_j))
                / np.sqrt(masses[i] * masses[j])
            )
        return D_blocks.swapaxes(1, 2).reshape(3 * N, 3 * N)
