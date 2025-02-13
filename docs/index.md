# Quesadilla
Quesadilla (**Que**ue-**sa**ving non**di**agona**l** super**la**ttices) is a python package for nondiagonal supercell (NDSC) phonon calculations, using the approach of [Lloyd-Williams and Monserrat](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.184301). 

The nondiagonal supercell method can dramatically reduce the computational cost of phonon calculations. In the standard "diagonal supercell" approach (as implemented in, e.g., phonopy), sampling the dynamical matrix on an $N \times N \times N$ $q$-grid requires calculating the forces in several $N \times N \times N$ supercells, each with $N^3$ as many atoms as the primitive cell. In the nondiagonal supercell approach, you only need the forces from supercells with $N$-times as many atoms as the primitive cell, leading to immense computational savings even for small systems.

Quesadilla is designed to significantly simplify this process. It can generate all the NDSCs and generate the necessary displacements for each supercell. Once you have used your favorite force calculator to compute the forces (e.g., VASP), quesadilla can compute the force sets in the nondiagonal supercells and then construct the full force constant matrix in real space for the corresponding diagonal supercell. Quesadilla uses symmetry and other tricks to figure out the absolute minimum number of NDSCs that covers the irreducible part of your $q$-grid, then uses little-group symmetries to generate the dynamical matrix in the entire Brillouin zone and symmetrize it, giving you the full force constant matrix in the equivalent diagonal supercell at a fraction of the computational cost.

Quesadilla is built as a layer on top of phonopy, so if you know how to use phonopy, you already know how to use Quesadilla. Further, quesadilla outputs `phonopy.yaml` files, so any workflows designed for phonopy will work flawlessly with quesadilla without any modification! Every force calculator supported by phonopy is supported by quesadilla as well. 

Under the hood, quesadilla uses `spglib` to handle much of the symmetry analysis, and a performant FORTRAN backend---adapted from Quantum ESPRESSO's PHonon routines---to deal with the symmetries of the dynamical matrix.

## Features
* A simple, phonopy-like interface for generating NDSCs and computing force constants.
* Supports all force calculators supported by phonopy (VASP, QE, Wien2k, etc.).
* Generates standard `phonopy.yaml` files that can be used with phonopy or any other package that reads phonopy output files.
* Heavily uses symmetry to minimize the computational cost of the calculations.
* Solves the set cover problem to extract the minimum number of NDSCs needed to cover the irreducible part of the $q$-grid.

### To be implemented
* Magnetic systems.
* Expanded documentation and tutorials.
* Developer documentation.

## Installation
Quesadilla requires a FORTRAN compiler and a BLAS/LAPACK installation on your path. It has only been tested on Linux and MacOS. The package is not yet listed on PyPI, so you need to install it directly from the GitHub repository via `pip`:
```bash
pip install git+git@github.com:oashour/quesadilla.git
``` 

If you have any issues with the installation, please open an issue.

On an Apple Silicon Mac, you may need to reinstall `phonopy` because it sometimes builds for the wrong architecture for some reason (not sure why...?).
```bash
CMAKE_ARGS="-DCMAKE_OSX_ARCHITECTURES=arm64" pip install --upgrade --verbose --force-reinstall --no-cache-dir phonopy
```

If you want to install the package in editable mode, it's slightly different from the usual approach because of the fortran extension.

```bash
pip install numpy meson meson-python
python -m pip install --no-build-isolation --editable . 
```

## Usage
Quesadilla has a simple CLI that is nearly identical to phonopy. If you know how to use phonopy, you already know how to use Quesadilla! See the full [tutorial](https://oashour.github.io/quesadilla/latest/tutorial.html) for more details.

To generate the NDSCs, you can run

```bash
quesadilla -d --dim 4 4 4
```
which will find the smallest number of NDSCs that cover the irreducible part of a 4x4x4 $q$-grid. You can then use your favorite force calculator to compute the forces in these supercells, and then run

```bash
quesadilla -f vasprun.xml
```
to create the `FORCE_SETS` file for each NDSC. Finally, you can compute the force constants in the full 4x4x4 supercell by running

```bash
quesadilla --fc
```

which will produce a `phonopy.yaml` file that can be used with phonopy or any other package that reads phonopy output files.