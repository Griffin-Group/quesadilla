#!/bin/bash

# Loop over directories sc-1 to sc-5
for i in {1..8}; do
    dir="sc-$i"
    
    if [ -d "$dir" ]; then
        cd "$dir" || exit
        
        rm POSCAR-* SPOSCAR *.hdf5 FORCE_SETS &> /dev/null
        phonopy make_disp.conf
        
        # Find all POSCAR-00x files
        for poscar_file in POSCAR-00*; do
            if [ -f "$poscar_file" ]; then
                # Extract the 3-digit number from the file name
                disp_dir="disp-${poscar_file##*-}"

                # Create the corresponding disp-00x directory
                mkdir -p "$disp_dir"
                cd $disp_dir
                rm -f OSZICAR PCDAT OUTCAR IBZKPT EIGENVAL CONTCAR DOSCAR CHG CHGCAR XDATCAR WAVECAR REPORT job.out job.err vasprun.xml vaspout.h5 PROCAR WAVEDER SOURCEPOT XCPOT dg.md job_details.out &> /dev/null; 
                cd ..

                # Copy the POSCAR file to the new directory as POSCAR
                cp "$poscar_file" "$disp_dir/POSCAR"

                # Copy INCAR, POTCAR, and job.sh from two levels up
                cp ../INCAR ../POTCAR ../job.sh "$disp_dir/"
                cp KPOINTS "$disp_dir"
            fi
        done
        cd ..
    else
        echo "Directory $dir does not exist."
    fi
done
