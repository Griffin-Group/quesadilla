#!/bin/bash
#SBATCH -A m3349
#SBATCH -J jobname
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH -N 4
#SBATCH --exclusive
#SBATCH -t 02:30:00
#SBATCH --mail-user=me@email.com
#SBATCH --mail-type=ALL

# This file was written by Omar A. Ashour (OAA) for Perlmutter
# Change log:
# 20231114, added file, (OAA)

# Some notes about the Slurm flags above
# -N is the number of nodes, Perlmutter CPU has 64 cores/socket, 2 sockets/node

# OpenMP Options (For advanced users, I recommend not changing this)
# I generally recommend disabling OpenMP (OMP_NUM_THREADS=1) unless you benchmark your
# system carefully. Using OpenMP means you can't use band group parallelism 
# (NCORE fixed to 1, ignoring what's in INCAR) and in most cases, you can get 
# better performance using an appropriate NCORE than with OpenMP.
# If you want to use OMP, use with 8 threads per core and don't use hyper threading/SMT.
# I suggest reading the NERSC white paper benchmarking VASP on Perlmutter.
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# VASP Module
# Use ISYM=-1 for magnetic systems if using VASP <= 6.3.2 on Perlmutter
module load vasp/6.4.3-cpu
# If you don't have access to the `vasp{5,6}` groups on NERSC, you can't use their 
# binaries from the module. Uncomment this line to use the Griffin group version
# Please read $CFS/m3349/codes/vasp/README.md if you use this version
# export PATH=$CFS/m3349/codes/vasp/vasp.6.4.1-cpu/bin/:$PATH

# --unbuffered ensures that VASP prints to stdout (output.log) and OUTCAR frequently
# -c should always be 2*OMP_NUM_THREADS
# -n should almost always be (Num_nodes * 128 / OMP_NUM_THREADS)
# For exceptions to setting -n this way, see below.
srun -n 64 -c 2 --cpu-bind=cores --unbuffered vasp_std &> job.out

# If you want to use fewer than 128 ranks per node (e.g., to give each MPI rank more memory, 
# or for better workload distribution if your numbers don’t add up properly), 
# you’ll also want to set --ntasks-per-node . For example, if you only want to use 64 MPI ranks per 
# node to double the memory available per rank (rarely necessary), 
# then you should do srun -n 256 -c 2 --ntasks-per-node=64 (assuming 4 nodes and OMP disabled)
# srun -n 256 -c 2 --ntasks-per-node=128 --cpu-bind=cores --unbuffered vasp_std &> output.log 
