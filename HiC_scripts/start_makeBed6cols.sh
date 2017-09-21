#!/bin/bash -l

#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=2:00:00


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Example command
#/projects/ag-papan/mapping/HiC_K562_public/makeBed6cols.sh /projects/ag-papan/mapping/HiC_K562_public/K562_SRR1658695_hic_hg19.bed


