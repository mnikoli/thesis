#!/bin/bash -l

if [ $# -ne 1 ]
  then
    echo "Usage: You need to provide the name of the BAM file as an argument"
    exit 1
fi

# First strand
input_file1=$1
base1=${1##*/}		# Get the file name from the path
prefix1=${base1%.*}	# Get the file name without the extension

dir=$input_file1	# Both input files should be in the same directory
parentdir="$(dirname "$dir")"

awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $input_file1 > $parentdir"/"$base1"_6cols.bed"
