# Quesadilla
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

# TODO:

- [x] Meson build system for the Fortran extension
- [x] Documentaton with Mkdocs
- [x] Automatic testing GitHub actions
- [x] Fix issue with nonsymmorphic space groups
- [x] Symmetrize the dynamical matrix using QE instead of spglib
- [x] Implement your Minkowski reduction and ensure it replicates Monserrat's results
- [x] Add unit tests
- [x] Refactor FORTRAN modules
- [x] Write calculation data to a log file of some kind
- [ ] Fix LAPACK issue with GitHub actions
- [ ] Significant code cleanup and refactoring
- [ ] Implement autoGR functionality for k-grids
- [ ] Implement CLI
- [ ] Documentation
- [ ] Implement magnetism
- [ ] Remove pmg as a dependency (use only phonopy + spglib)
- [ ] Test integer linear programming reduction of necessary supercells