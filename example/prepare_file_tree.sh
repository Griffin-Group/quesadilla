#!/bin/bash

# This bash script creates the necessary file tree for Quesadilla calcs

for i in {1..5}; do
    dir="sc-00$i"
    
    if [ -d "$dir" ]; then
        cd "$dir" || exit

        # Second argument is k-point spacing in 1/A
        python ../generate_kpoints.py SPOSCAR 0.15

        # Pristine supercell
        disp_dir="disp-000"
        mkdir -p $disp_dir
        echo "I made directory $disp_dir"
        cp SPOSCAR "$disp_dir/POSCAR"
        cp ../INCAR ../POTCAR ../job.sh "$disp_dir/"
        cp KPOINTS "$disp_dir"
        
        for poscar_file in POSCAR-00*; do
            if [ -f "$poscar_file" ]; then
                disp_dir="disp-${poscar_file##*-}"
                mkdir -p "$disp_dir"
                echo "I made directory $disp_dir"
                cd $disp_dir
                cd ..

                cp "$poscar_file" "$disp_dir/POSCAR"

                # Copy INCAR, POTCAR, and job.sh from two levels up
                cp ../INCAR ../POTCAR ../job.sh "$disp_dir/"
                cp KPOINTS "$disp_dir"
            fi
        done
        echo "-------------------------------------"
        cd ..
    else
        echo "Directory $dir does not exist."
    fi
done
