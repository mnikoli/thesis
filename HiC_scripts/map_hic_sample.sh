#!/bin/bash -l

### Put files from both strands in the same folder! ###

if [ $# -ne 5 ]
  then
    echo "Usage: You need to provide 5 arguments"
    echo "A note: Put files from both strands in the same folder!" 
    echo "Command:./b.sh <strand1_to_be_mapped> <strand2_to_be_mapped> <name> <read_length> <ref_genome>"
    echo "EXAMPLE: ./b.sh  TNF1.fastq.gz TNF2.fastq.gz TNF 100"
    exit 1
fi

# First strand
input_file1=$1
base1=${1##*/}		# Get the file name from the path
prefix1=${base1%.*}	# Get the file name without the extension

# Second strand
input_file2=$2
base2=${2##*/}		# Get the file name from the path
prefix2=${base2%.*}	# Get the file name without the extension

dir=$input_file1	# Both input files should be in the same directory
parentdir="$(dirname "$dir")"

# Initial paramaters
i=2
name=$3
rlen=19
trim_size=$(($4 - rlen))
ref_genome=$5

#####################################
# Map first with BWA MEM -- START --#
#####################################
/projects/ag-papan/install/bwa_0.7.12/bwa-0.7.12/bwa mem -t 12 -T 19 -M $ref_genome $input_file1 $input_file2 > $parentdir"/"$name"_it1.sam"

### Create a BAM file from both strands
/projects/ag-papan/install/samtools-1.2/samtools view -Sb $parentdir"/"$name"_it1.sam" > $parentdir"/"$name"_it1.bam"

### Sort BAM # Do I really need this???
/projects/ag-papan/install/samtools-1.2/samtools sort -@ 12 $parentdir"/"$name"_it1.bam" $parentdir"/"$name"_it1_sorted"

### Get unmapped reads from BAM
/projects/ag-papan/install/samtools-1.2/samtools view -b -f 4 $parentdir"/"$name"_it1_sorted.bam" > $parentdir"/"$name"_it1_unmapped.bam"

### Unmapped BAM to unmapped FASTQ
/projects/ag-papan/install/BEDTools_2.17.0/bedtools-2.17.0/bin/bamToFastq -i $parentdir"/"$name"_it1_unmapped.bam" -fq $parentdir"/"$name"_it1_unmapped1.fastq" -fq2 $parentdir"/"$name"_it1_unmapped2.fastq" 

### Get mapped reads from BAM (the final product from this part)
/projects/ag-papan/install/samtools-1.2/samtools view -b -F 4 $parentdir"/"$name"_it1_sorted.bam" > $parentdir"/"$name"_it1_mapped.bam"

input_file1=$parentdir"/"$name"_it1_unmapped1.fastq"
input_file2=$parentdir"/"$name"_it1_unmapped2.fastq"

###--END--###

##############################################
# Try to map the rest of the reads --START --# 
##############################################

### Trim both strands
/projects/ag-papan/install/homer/bin/homerTools trim -3 $trim_size $input_file1
/projects/ag-papan/install/homer/bin/homerTools trim -3 $trim_size $input_file2	

### Map the first strand
/projects/ag-papan/install/bwa_0.7.12/bwa-0.7.12/bwa aln -t 12 $ref_genome $input_file1".trimmed" > $parentdir"/"$prefix1"_it"$i".sai"

### Map the second strand
/projects/ag-papan/install/bwa_0.7.12/bwa-0.7.12/bwa aln -t 12 $ref_genome $input_file2".trimmed" > $parentdir"/"$prefix2"_it"$i".sai"

### Create a SAM file from both strands
/projects/ag-papan/install/bwa_0.7.12/bwa-0.7.12/bwa sampe $ref_genome $parentdir"/"$prefix1"_it"$i".sai" $parentdir"/"$prefix2"_it"$i".sai" $input_file1".trimmed" $input_file2".trimmed" -f $parentdir"/"$name"_it"$i".sam"
	
### SAM to BAM
/projects/ag-papan/install/samtools-1.2/samtools view -Sb $parentdir"/"$name"_it"$i".sam" > $parentdir"/"$name"_it"$i".bam"

### Sort BAM
/projects/ag-papan/install/samtools-1.2/samtools sort -@ 12 $parentdir"/"$name"_it"$i".bam" $parentdir"/"$name"_it"$i"_sorted"
 
### Get mapped reads from BAM (the final product from this part)
/projects/ag-papan/install/samtools-1.2/samtools view -b -F 4 $parentdir"/"$name"_it"$i"_sorted.bam" > $parentdir"/"$name"_it"$i"_mapped.bam"

###--END--###

################################################
# Merge 2 BAM files with mapped reads --START--#
################################################

/projects/ag-papan/install/samtools-1.2/samtools merge -n $parentdir"/"$name"_merged.bam" $parentdir"/"$name"_it1_mapped.bam" $parentdir"/"$name"_it2_mapped.bam"

/projects/ag-papan/install/samtools-1.2/samtools sort -@ 12 $parentdir"/"$name"_merged.bam" $parentdir"/"$name"_merged_sorted"

java -Xmx16g -Djava.io.tmpdir=/projects/ag-papan/install/picard-tools-1.119/tmp -jar /projects/ag-papan/install/picard-tools-1.119/MarkDuplicates.jar INPUT=$parentdir"/"$name"_merged_sorted.bam" REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT OUTPUT=$parentdir"/"$name"_merged_nodup.bam" M="report_"$name".txt"

#--END--#

/projects/ag-papan/install/samtools-1.2/samtools flagstat $parentdir"/"$name"_merged_nodup.bam"
/projects/ag-papan/install/samtools-1.2/samtools sort -n -@ 12 $parentdir"/"$name"_merged_nodup.bam" $parentdir"/"$name"_merged_nodup_sorted"
/projects/ag-papan/install/BEDTools_2.17.0/bedtools-2.17.0/bin/bamToBed -bedpe -i $parentdir"/"$name"_merged_nodup_sorted.bam" > $parentdir"/"$name".bed"

# Clean up a bit
rm *.sai
rm *_unmapped.bam
rm *_unmapped.fastq
rm *fastq.trimmed


