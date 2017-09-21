library(IRanges)
library(GenomicRanges)
library(Matrix)
library(HiTC)
library(MASS)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(tools) # get filename without extension

options(scipen=999) #exponential notation off

##################
# Read the input #
##################

# Input files
chr_len<-read.table("/projects/ag-papan/Heatmaps_HiC/inputR/hg19.chrom.sizes")

args=(commandArgs(trailingOnly = TRUE))
file=args[1]

filename=file_path_sans_ext(basename(file))

sample<-fread(file, data.table=FALSE)
chr_len<-chr_len[,1]

# Initial parameters
bin_size<-10e3 # Bin size
n<-0 # Number of bins
chr_id<-1

##############################################
# Bin the genome and make one GRanges object #
##############################################

bins.all<-NULL 
HTC_list<-NULL

for (chr in 1:length(chr_len)){
  # Break each chromosome into bins of a given size
  chunks<-breakInChunks(chr_len[chr_id], bin_size)
  n=n+length(chunks) 
  
  # Create chromosome names for each chunk
  seqnames<-rep(paste("chr", chr, sep=""), length(chunks))  
  
  # Merge chromosome names column with chunks  
  tmp.df<-cbind(data.frame(seqnames), chunks)
  
  # Add tmp.df data frame to the resulting bins variable bins  
  bins.all<-rbind(bins.all, tmp.df)
  
  # Move to the next chromosome
  chr_id=chr_id+1
}

bins.all<-unique(bins.all)
names<-paste('bin', 1:nrow(bins.all), sep="_")
bins.all<-GRanges(seqnames=bins.all[[1]], ranges=IRanges(start=bins.all[[2]], end=bins.all[[3]], names=names))

seqlevels(bins.all)<-sub('chr23','chrX', seqlevels(bins.all))
seqlevels(bins.all)<-sub('chr24','chrY', seqlevels(bins.all))

# Split GRanges object to GRangesList
bins.list<-split(bins.all, seqnames(bins.all))

###################
# Parse the input #
###################

# Get the start ("from") positions from BEDPE input   
from<-sample[,1:3]
names(from)<-c("chr","start","end")
from.gr<-makeGRangesFromDataFrame(from)

# Get the end ("to") positions from BEDPE input   
to<-sample[,4:6]
names(to)<-names(from)
to.gr<-makeGRangesFromDataFrame(to)

rm(sample)
gc()

###############################
# Assign interactions to bins #
###############################

for (from.chr in 1:length(chr_len)){
  for (to.chr in 1:length(chr_len)){
    if(from.chr != to.chr){
      M<-matrix(0, nrow=nrow(as.data.frame(bins.list[from.chr])), ncol=nrow(as.data.frame(bins.list[to.chr])))

      cat ("Calculating interaction matrix between chromosomes",from.chr," and ", to.chr, "\n")
      # Find interactions starting from bin i
      interactions.from<-findOverlaps(from.gr, unlist(bins.list[from.chr]), ignore.strand=TRUE)
      index.from.bed<-queryHits(interactions.from)  # get line indexes from the bed file that have interactions from the starting current chr
      index.from.bin<-subjectHits(interactions.from)  # bin ids of interacting start points from above
  
      # From interactions starting from bin i, find bins that they interact with (j coord. of the interacting matrix)
      sub.to.gr<-to.gr[index.from.bed,] # subset of intervals starting interactions only from the current chromosome
      interactions.to<-findOverlaps(sub.to.gr, unlist(bins.list[to.chr]), ignore.strand=TRUE)

      if (length(interactions.to)!=0){
        index.to.bin<-subjectHits(interactions.to)  
  
        cc<-as.data.frame(sub.to.gr) # all possible ("to") positions from the current chr
        cc1<-cbind(cc, index.from.bin) #starting bin id
        # cc1 now contains all contacts starting from the current start chr, but end chr is not yet fixed!
  
        cc1.sub<-cc1[queryHits(interactions.to),]
        cc2<-cbind(cc1.sub, subjectHits(interactions.to))

        # Populating matrix M as requested for HTCexp object (y intervals are put in rows, x intervals in columns)
        for (k in 1:nrow(cc2)){  
          M[cc2[[6]][k], cc2[[7]][k]]=M[cc2[[6]][k], cc2[[7]][k]]+1
        }

        a<-bins.list
        b<-bins.list

        seqlevels(a, force=TRUE)<-seqlevels(bins.list)[from.chr]
        seqlevels(b, force=TRUE)<-seqlevels(bins.list)[to.chr]
        
        o<-strsplit(names(a[[1]][1]), split="_")
        a_start<-as.integer(o[[1]][2])
        o<-strsplit(names(a[[1]][length(a[[1]])]), split="_")
        a_end<-as.integer(o[[1]][2])
        o<-strsplit(names(b[[1]][1]), split="_")
        b_start<-as.integer(o[[1]][2])
        o<-strsplit(names(b[[1]][length(b[[1]])]), split="_")
        b_end<-as.integer(o[[1]][2])

        a<-makeGRangesFromDataFrame(as.data.frame(a))
        names(a)<-names[a_start:a_end]
      
        b<-makeGRangesFromDataFrame(as.data.frame(b))
        names(b)<-names[b_start:b_end]

        colnames(M)<-names(b)
        rownames(M)<-names(a)
        MM<-Matrix(M)

        f<-HTCexp(MM, b, a, forceSymmetric=TRUE)
  
        if (length(HTC_list) != 0){
          HTC_list<-c(HTC_list, HTClist(f))
        } else {
          HTC_list<-HTClist(f)
        }
      }
     }
  }  
} 

################################################
# Make the heatmap symmetrical & plot raw data #
################################################
filename.out<-paste("/projects/ag-papan/mapping/HiC_K562_public/", filename, sep="")
saveRDS(HTC_list, file=paste(filename.out, "_inter.rds", sep=""))
