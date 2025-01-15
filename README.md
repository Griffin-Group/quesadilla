# Quesadilla
[![Pre-Alpha](https://img.shields.io/badge/Status-Pre--Alpha-red)](https://oashour.github.io/quesadilla/develop/)
[![GitHub Release](https://img.shields.io/github/v/release/oashour/quesadilla?include_prereleases)](https://github.com/oashour/quesadilla/releases)
[![Tests](https://github.com/oashour/quesadilla/actions/workflows/run_tests.yaml/badge.svg)](https://github.com/oashour/quesadilla/actions)
[![License](https://img.shields.io/badge/License-GPL-blue)](#license "Go to license section")
[![Stable Docs](https://img.shields.io/badge/Docs-Stable-blue)](https://oashour.github.io/quesadilla/latest/)
[![Develop Docs](https://img.shields.io/badge/Docs-Develop-purple)](https://oashour.github.io/quesadilla/develop/)


Quesadilla (**Que**ue-**sa**ving non**di**agonal super**la**ttices) is a python package for nondiagonal supercell phonon calculations, using the approach of [Lloyd-Williams and Monserrat](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.184301)



# Installation

In principle, just running this (in a `venv`) should work, assuming a fortran compiler and `BLAS`/`LAPACK` implementation are available:
```bash
pip install git+git@github.com:oashour/quesadilla.git
```

This works fine on Perlmutter. On an Apple Silicon Mac, you may need to reinstall `phonopy` because it builds for the wrong architecture for some reason (I think a bug with their build system).

```bash
CMAKE_ARGS="-DCMAKE_OSX_ARCHITECTURES=arm64" pip install --upgrade --verbose --force-reinstall --no-cache-dir phonopy
```

If you want to install the package in editable mode, it's slightly different from the usual approach because of the fortran extension.

```bash
python -m pip install --no-build-isolation --editable . 
```

# Usage
Currently the CLI does not exist, but the package is fully usable via the python API (which needs quite a bit of work, especially to support IO for all the codes phonopy supports besides VASP). See the `examples` directory for a Silicon example on a 4x4x4 $q$-grid. You only need to read the Jupyter notebook and follow the instructions.
