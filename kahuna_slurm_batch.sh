#!/bin/bash
# Usage: sbatch kahuna_slurm.batch "[matlab script]"

# --- Send the output to "out.txt" ---
#SBATCH -o out_%A.txt

# --- Use the "compute" queue ---
#SBATCH -p compute

# --- Set the max time to 10 hours ---
#SBATCH -t 5-10:00:00

# --- Use one single node, with 1 process, and 4 cores ---
# --- Plus exclusively claim the resources ---
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --exclusive
#SBATCH -N 1


module load matlab
echo Hostname:
hostname
echo MATLAB Command: $1
matlab -nosplash -nodesktop -nodisplay -r "addpath('/home/bwlarse/tensor_toolbox_sparse'); addpath('/home/bwlarse/matlab-tools'); $1; quit"
