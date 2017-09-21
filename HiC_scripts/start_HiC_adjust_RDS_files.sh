#!/bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --mem=46gb
#SBATCH --time=60:00:00

module load intel/15.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

base1=${1##*/}
base2=${base1%.*}

echo "Submitting job for files : $1 and $2"
Rscript --vanilla /projects/ag-papan/Heatmaps_HiC/HiC_adjust_RDS_files.R $1 $2
