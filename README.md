# Quesadilla
[![Pre-Alpha](https://img.shields.io/badge/Status-Pre--Alpha-red)](https://oashour.github.io/quesadilla/develop/)
[![GitHub Release](https://img.shields.io/github/v/release/oashour/quesadilla?include_prereleases)](https://github.com/oashour/quesadilla/releases)
[![Tests](https://github.com/oashour/quesadilla/actions/workflows/run_tests.yaml/badge.svg)](https://github.com/oashour/quesadilla/actions)
[![License](https://img.shields.io/badge/License-GPL-blue)](#license "Go to license section")
[![Stable Docs](https://img.shields.io/badge/Docs-Stable-blue)](https://oashour.github.io/quesadilla/latest/)
[![Develop Docs](https://img.shields.io/badge/Docs-Develop-purple)](https://oashour.github.io/quesadilla/develop/)


Quesadilla (**Que**ue-**sa**ving non**di**agonal super**la**ttices) is a python package for nondiagonal supercell phonon calculations, using the approach of [Lloyd-Williams and Monserrat](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.184301)

## Features
* A simple, phonopy-like interface for generating NDSCs and computing force constants.
* Supports all force calculators supported by phonopy (VASP, QE, Wien2k, etc.).
* Generates standard `phonopy.yaml` files that can be used with phonopy or any other package that reads phonopy output files.
* Heavily uses symmetry to minimize the computational cost of the calculations.
* Solves the set cover problem to extract the minimum number of NDSCs needed to cover the irreducible part of the $q$-grid.

## Installation

Quesadilla requires a FORTRAN compiler and a BLAS/LAPACK installation on your path. It has only been tested on Linux and MacOS. The package is not yet listed on PyPI, so you need to install it directly from the GitHub repository via `pip`:
```bash
pip install git+git@github.com:oashour/quesadilla.git
``` 

If you have any issues with the installation, please open an issue.


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