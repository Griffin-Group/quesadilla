#!/bin/bash

# Loop over directories sc-1 to sc-5
for i in {1..5}; do
    dir="sc-00$i"
    
    cd "$dir"
    echo "Working on $dir"
        
    for disp_dir in disp-*; do
        cd $disp_dir || exit
        ./job.sh &
        sleep 3 # To avoid issues with multiple job steps running too quickly
        cd ..
    done
    wait
    echo "Directory $disp_dir is done"
    cd ..
done
