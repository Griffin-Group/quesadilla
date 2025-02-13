# Si tutorial

See the tutorial here <add link>. You will need to add the correct `POTCAR` file to this directory first, then you can follow the tutorial. Make sure you modify the `job.sh` file for your cluster/machine.

When you get to step 2, you can run
```
# Generates dips-* directories, correct KPOINTS, copies files, etc
./prepare_file_tree.sh 
# Runs ./job.sh in each dips-* directory
./run.sh
```