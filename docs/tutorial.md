# Tutorial

The quesadilla command line interface is very similar to phonopy's. Regardless of what DFT code you're using, the steps are more or less the same.

## Step 1: Generate supercells with displacements

With phonopy, you would run something like `phonopy -d --dim 4 4 4`, perhaps with extra arguments if you don't use VASP (e.g., `phonopy --qe -d --dim 4 4 4 -c pw.in`). The first step with quesadilla is to use the exact same command, but with `phonopy` replaced with `quesadilla`:

```bash
# Assumes you use VASP and POSCAR is in current directory
quesadilla -d --dim 4 4 4
```
or 
```bash
# Quantum espresso example
quesadilla --qe -d --dim 4 4 4 -c pw.in
```

quesadilla also supports all the other `phonopy` arguments you may want to use when generating the displacements, such as `--pm`, `--trigonal`, etc.

In a conventional `phonopy` calculation, this would have produced a pristine 4x4x4 supercell file (e.g., `SPOSCAR`) as well as some supercells with displacements (e.g., `POSCAR-00{1,2,3...}`) and a `phonopy_disp.yaml`. With quesadilla, the situation is quite similar but with some extra steps:

1. First, quesadilla will find a standard primitive cell. This is similar to running `phonopy --symmetry`. This standard primitive cell is what will be used to build the nondiagonal supercell, and will be saved in a file `PPOSCAR`.
2. Next, quesadilla will figure out the minimum number of nondiagonal supercells necessary to span the 4x4x4 q-grid, build them, and then Minkowski reduce them. These files will be placed in the directories `sc-{001,002,...}`, each corresponding to a nondiagonal supercell. Each NDSC will be in, e.g., the `sc-*/SPOSCAR` file for VASP. The standard primitive `POSCAR` file will also be placed in these directories.
3. Finally, quesadilla will generate the displacemnts in each of these supercells, producing the files `sc-*/POSCAR-*`. This is similar to running `phonopy -d --dim <ndsc matrix>` on every single supercell. 

Each `sc-*` folder will also have a `phonopy_disp.yaml` file, which can be used for debugging if necessary, and will be read later when creating the force sets. In the root directory where the primitive cell is located, quesadilla will also write a file `quesadilla.toml` with a log of the number of supercells, their sizes, the q-points they are commensurate with, etc.

## Step 2: run the DFT calculations

As usual with `phonopy`, we now have to run the DFT calculations on the supercells with the displacements. We will place each `POSCAR-x` in a directory `disp-x` and run the calculation there (e.g., `POSCAR-002` should go in a directory `disp-002`). While this directory structure is optional for phonopy, it is mandatory for quesadilla. This needs to be repeated for every single supercell (i.e., for each `sc-{i}` directory). Optionally, the forces in the pristine supercell can also be calculated so they can be subtracted from the force sets in the next step. If this is desired, the pristine supercell calculation (e.g., the `SPOSCAR` file) should be run in the `disp-000` folder.

Note that, since nondiagonal supercells have lower symmetry, you may require more displacements than you would for a diagonal cell. For example, a 4x4x4 diagonal supercell of Si will only require the calculation of one supercell with displacements, whereas with nondiagonal supercells we will need 4 supercells with 2 displacements each and 1 supercell with 4 displacements.

## Step 3: generate the force sets (`FORCE_SETS`)

The next step is to generate the `FORCE_SETS` file. In phonopy, we would use something like `phonopy disp-00{1..4}/vasprun.xml` to read the forces from the `vasprun.xml` for the four supercells with displacements. Here, we have five nondiagonal supercells with varying numbers of displacements each. As mentioned earlier, your calculations should always be structured so that `POSCAR-001` (and the resultant `vasprun.xml` or similar output file) go into the folder `disp-001`, etc. To generate the `FORCE_SETS` for all supercells with `quesadilla`, we can run from the root directory

```bash
quesadilla -f vasprun.xml
```

which is similar to running `phonopy -f disp-*/vasprun.xml` in every `sc-*` directory, generating one `FORCE_SETS` file per supercell. If you want to use the `--force-sets-zero` option, subtracting the residual forces from the pristine supercell, then with quesadilla we'll run

```bash
quesadilla --force-sets-zero vasprun.xml
```

assuming the pristine supercell calculation was done in the `sc-*/disp-000` directories as explained in the previous section.

## Step 4: generate the force constants

Finally, we need to generate the force constants. This is a little bit different from `phonopy` where often you won't generate them explicitly. The command is (from the root directory)

```bash
quesadilla --fc
```

which will do the following:
1. From the force sets (`FORCE_SETS`) file in each `sc-*` directory, it will compute the real-space force constants for this specific NDSC.
2. Next, for each NDSC, it will compute the dynamical matrix at every q-point that supercell it is commensurate with.
3. If a q-point is commensurate with multiple cells (such as the $\Gamma$ point, commensurate with all of them), the dynamical matrix is averaged. We now have the dynamical matrix at every $q$ point in the iredducible Brillouin zone.
4. Next, quesadilla will symmetrize the reciprocal-space dynamical matrix using the symmetries of the little group of each q-point. This can get rid of small symmetry breaking problems that can happen in non-diagonal supercells.
5. Using the symmetrized dynamical matrices, quesadilla will calculate the star of each q-point in the iredducible brillouin zone, then use space group symmetries, and time reversal symmetry if applicable, to compute the dynamical matrix at every point in the star of q.
6. Once this is done, we now have the dynamical matrix in the full Brillouin zone. We can now Fourier transform it to real space to obtain the force constant matrix of a diagonal 4x4x4 supercell.
7. Finally, quesadilla will apply the acoustic sum rules (ASR) to enforce translational symmetry.

As a result of this command, a `phonopy.yaml` file containing all the calculation information and the force constants will be produced in the root directory. This file can now be used with `phonopy`. `phono3py`, DarkMAGIC, or any other code that reads `phonopy.yaml` files. Since it contains the full force constants, we don't need a `FORCE_CONSTANTS` or `force_constants.hdf5` file, but they can be dumped explicitly from this `phonopy.yaml` by using `phonopy` itself.
