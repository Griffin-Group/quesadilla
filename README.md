# quesadilla
Quesadilla (*Que*ue-*Sa*ving non*di*agonal super*la*ttices) is a python package for nondiagonal supercell phonon calculations, using the approach of [Lloyd-Williams and Monserrat](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.184301)

- [x] Meson build system for the Fortran extension
- [x] Documentaton with Mkdocs
- [x] Automatic testing GitHub actions
- [ ] Fix issue with nonsymmorphic space groups
- [ ] Symmetrize the dynamical matrix using QE instead of spglib
- [ ] Implement your Minkowski reduction and ensure it replicates Monserrat's results
- [ ] Significant code cleanup and refactoring
- [ ] Implement autoGR functionality for k-grids
- [ ] Implement CLI
- [ ] Add unit tests
- [ ] Documentation
- [ ] Implement magnetism
- [ ] Remove pmg as a dependency (use only phonopy + spglib)
- [ ] Test integer linear programming reduction of necessary supercells