#!/bin/bash

# This script is the first step (second is to use rs3.R and then rs4.R scripts) for post-processing of format output from scripts from publication 
# "Kelly et all 2015 - Splicing of many human genes involves sites embedded within introns."
# Paths need to be hard-codded in all three scripts


# ./gen_garlic_input_from_rsites.sh /home/milos/R/recursive_splicing/polyA_results
#./gen_garlic_input_from_rsites.sh /home/milos/R/recursive_splicing/RS_results

path=$1

file="$path/new_results_Adrenal_gland.xls"

#for file in $path/new_results_d*_merged.xls;
#for file in $path/*.xls;
#do
# use -c60- for polyA folder and -c57- for everything else (not corrected for polyA)
    name="$(echo $file | cut -c57- | rev | cut -c5- | rev)"
    echo $name
    threeprime="$path/three_prime_to_cis_intron_$name.txt"
    echo $threeprime
    
    ssconvert $file "new_results_"$name".csv"
    awk -F "," '{print $1"\t"$15"\t"$18}' "new_results_"$name".csv" > "additional_2_exon_cols_"$name".csv"
    awk -F "," '{print $1"\t"$5"\t"$12}' "new_results_"$name".csv" > "additional_2_exon_cols1_"$name".csv"

    Rscript ~/R/recursive_splicing/rs3_final3.R $threeprime "additional_2_exon_cols_"$name".csv" "additional_2_exon_cols1_"$name".csv"


    rm "additional_2_exon_cols_"$name".csv"
    rm "additional_2_exon_cols1_"$name".csv"
    rm "new_results_"$name".csv"
    
#done
