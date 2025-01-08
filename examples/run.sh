#!/bin/bash

alias vaspVeryClean='rm -f OSZICAR PCDAT OUTCAR IBZKPT EIGENVAL CONTCAR DOSCAR CHG CHGCAR XDATCAR WAVECAR REPORT job.out job.err vasprun.xml vaspout.h5 PROCAR WAVEDER SOURCEPOT XCPOT dg.md job_details.out &> /dev/null'

# Loop over directories sc-1 to sc-5
for i in {1..8}; do
    dir="sc-$i"
    
    cd "$dir"
    echo "Working on $dir"
        
    for disp_dir in disp-*; do
        cd $disp_dir || exit
        ./job.sh &
        cd ..
    done
    wait
    echo "Directory $disp_dir is done"
    phonopy -f disp-*/vasprun.xml
    phonopy get_yaml.conf
    cd ..
done
