## Simulation code

Code used to perform simulations in SLiM v.4.0.

The scripts in [hpc_scripts](https://github.com/pinskylab/codPolyEvol/tree/main/slim_code/hpc_scripts) were used to perform simulations in parallel on Rutgers' Amarel HPCC. If you have a SLURM-type system, using 'sbatch' to run the .sb files should execute the corresponding .py file to generate 20 replicates for each migration scenario (though it may need to be tweaked for your particular system). process_slimouts.sh can then be used to calculate allle frequencies and Fsts. 

[codqtl_reduced_7k_3type_gammarec.slim](https://github.com/pinskylab/codPolyEvol/blob/main/slim_code/codqtl_reduced_7k_3type_gammarec.slim) uses the same demographic template but generates 3 different kinds of mutations (neutral, deleterious, and QT). 
