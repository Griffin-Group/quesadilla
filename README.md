# quesadilla
Quesadilla (**Que**ue-**Sa**ving non**di**agonal super**la**ttices) is a python package for nondiagonal supercell phonon calculations, using the approach of [Lloyd-Williams and Monserrat](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.184301)

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
