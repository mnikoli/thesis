#!/bin/bash -l

#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mem=32gb
#SBATCH --time=96:00:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Example command
#/projects/ag-papan/mapping/HiC_K562_public/map_hic_sample.sh /projects/ag-papan/mapping/HiC_K562_public/K562_101_hic_1.fastq /projects/ag-papan/mapping/HiC_K562_public/K562_101_hic_2.fastq K562_101_hic_hg19 101 /projects/ag-papan/genomes/hg19/hs_ref_GRCh37_p5.fa




