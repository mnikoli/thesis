#!/bin/bash -x

# Example path
file=/projects/ag-papan/mapping/HiC_K562_public/K562_SRR1658702_hic_hg19.bed_6cols.bed

#echo "Submitting files for binning of intra- and interchromosomal interactions: "

sbatch /projects/ag-papan/Heatmaps_HiC/start_assignFragmentsToBins_intra.sh $file	
sbatch /projects/ag-papan/Heatmaps_HiC/start_assignFragmentsToBins_inter1.sh $file	

