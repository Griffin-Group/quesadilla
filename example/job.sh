#!/bin/bash

# You can run this in an interactive session on a cluster

export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

module load vasp/6.4.3-cpu

srun -n 64 -c 2 --cpu-bind=cores --unbuffered vasp_std &> job.out
