########################################################################################
# Script to get 2 consecutive CTCF loops and plot their interaction pilups genome wide # 
########################################################################################

library(GenomicRanges)
library(Matrix)
library(HiTC)
library(MASS)
library(data.table)
library(tools) #get the filename without the extension

options(scipen=999) #exponential notation off

# Read the list of HMGB2 peaks on 80kb extended borders
### HUVEC ###
peak_list<-"/media/milos/NAS/UCSC_data/GSE63525_HUVEC_HiCCUPS_looplist.bed"
# peak_list<-"/media/milos/Data/senescence_final/HMGB2_on_80kb_ext_TAD_borders_HUVy.csv"
# peak_list<-"/media/milos/Data/senescence_final/adjust_peaks/HMGB2_HUVEC_peaks_final.csv"

# peak_list<-"/media/milos/Data/senescence_final/CTCF_on_20kb_TAD_borders_40kb_ext.csv"
# peak_list<-"/media/milos/NAS/UCSC_data/ENCODE_ChIPseq_BED_files_hg18/hg19_marks/CTCF.HUVEC_fdr0.001_hg19.csv"

# peak_list<-"/home/milos/CTCF_on_HMGB2_peaks_HUV.csv"
# peak_list<-"/home/milos/HMGB2_and_CTCF_on_TAD_borders_HUV.csv"


### IMR90 ###
# peak_list<-"/media/milos/Data/senescence_final/HMGB2_on_80kb_ext_20kbTAD_borders_IMR90y.csv"
# peak_list<-"/media/milos/NAS/UCSC_data/GSE63525_IMR90_HiCCUPS_looplist.bed"
# peak_list<-"/media/milos/Data/senescence_final/adjust_peaks/IMR90_p10_peaks_final.csv"
# peak_list<-"/media/milos/NAS/AG_Papantonis/senescence_project/HiChIP/output_run1/tmp/66189.intra.loop_counts_g150k.bedpe" # young
# peak_list<-"/media/milos/NAS/AG_Papantonis/senescence_project/HiChIP/output_run1/tmp/66190.intra.loop_counts_g150k.bedpe"   # old

# p.all<-fread(peak_list, data.table=F)

### Adjustment for CTCF loops
p.all<-read.table(peak_list, header=FALSE, sep="\t", stringsAsFactors = FALSE )
indexes<-which((p.all[,5]-p.all[,3])>=150000)

tmp.1<-as.data.frame((p.all[indexes,2]+p.all[indexes,3])/2)
tmp.2<-as.data.frame((p.all[indexes,5]+p.all[indexes,6])/2)

p.from<-cbind(p.all[indexes,1], tmp.1,(tmp.1+1))
p.to<-cbind(p.all[indexes,4], tmp.2,(tmp.2+1))


# tmp.1<-as.data.frame((p.all[,2]+p.all[,3])/2)
# tmp.2<-as.data.frame((p.all[,5]+p.all[,6])/2)
# 
# p.from<-cbind(p.all[,1], tmp.1,(tmp.1+1))
# p.to<-cbind(p.all[,4], tmp.2,(tmp.2+1))

names(p.from)<-c("c", "s", "e")
names(p.to)<-c("c", "s", "e")

### Load the HiC HUVECy MERGED sample at 20kb resolution
# hic_file<-"/home/milos/R/Dy_MERGED_10Kb_cis.rds"
# hic_file<-"/home/milos/R/Do_MERGED_10Kb_cis.rds"

### Corrected HUVEC samples (merged)
# hic_file<-"/home/milos/R/HUVy_merged_corrected_normAnnot_10kb_cis.rds"
# hic_file<-"/home/milos/R/HUVo_merged_corrected_normAnnot_10kb_cis.rds"

# hic_file<-"~/R/Dy_MERGED_20Kb_cis.rds"
# hic_file<-"~/R/Do_MERGED_20Kb_cis.rds"
# hic_file<-"/home/milos/R/NTC_HiC_hg19_20Kb.rds"
# hic_file<-"/home/milos/R/SiHMGB2_HiC_hg19_20Kb.rds"

# hic_file<-"/media/milos/TOSHIBA EXT/HiC_10kb_separated_samples/Anne_Feb17/NTC_HiC_hg19.bed_6cols_intra_NormAnnotated_hg19_feb17.rds"
hic_file<-"/media/milos/TOSHIBA EXT/HiC_10kb_separated_samples/Anne_Feb17/SiHMGB2_HiC_hg19.bed_6cols_intra_NormAnnotated_hg19_feb17.rds"
### IMR90 ###
# hic_file<-"~/R/IMR90o_HiC_cis_20Kb.rds"
# hic_file<-"/home/milos/R/IMR90y_MERGED_10Kb_cis.rds"
# hic_file<-"/home/milos/R/IMR90o_MERGED_10Kb_cis.rds"

hic<-readRDS(hic_file)
filename<-file_path_sans_ext(basename(hic_file))

eps<-5 # Epsilon neighbourhoud (in bins)
n<-2*eps+1

### Find peak coordinates in hic matrix
M.res<-matrix(rep(0,n^2), nrow=n, ncol=n)

for (curr_chr in 1:length(hic)){
  chr<-paste("chr",curr_chr,sep="")
  bins<-as.data.frame(x_intervals(hic[[curr_chr]]))
  
  A<-intdata(hic[[curr_chr]])
  
  if (curr_chr==23){ chr<-"chrX"}
  if (curr_chr==24){ chr<-"chrY"}
  
  # get the list index for the current chr
  p.f<-p.from[p.from[,1]==chr,]
  p.t<-p.to[p.to[,1]==chr,]
  # get the number of peaks on the current chr
  (num_peaks<-nrow(p.f))
  for (i in 1:num_peaks){
    (bin.index.f<-which(row.names(bins) == row.names(bins[(bins[,2]<=p.f[i,2] & bins[,3]>=p.f[i,3]),])))
    (row.names(bins[(bins[,2]<=p.f[i,2] & bins[,3]>=p.f[i,3]),]))
    (bin.index.t<-which(row.names(bins) == row.names(bins[(bins[,2]<=p.t[i,2] & bins[,3]>=p.t[i,3]),])))
    (row.names(bins[(bins[,2]<=p.t[i,2] & bins[,3]>=p.t[i,3]),]))
    # if (length(bin.index)>0 && bin.index-eps>0 && (bin.index+eps)<=nrow(bins) && bin.index.f>0 && bin.index.t>0){
    if (length(bin.index.f)>0 && length(bin.index.t)>0 && bin.index.f-eps>0 && bin.index.t-eps>0 && (bin.index.f+eps)<=nrow(bins) && (bin.index.t+eps)<=nrow(bins)){
      if (A[bin.index.f,bin.index.t]>0){
        # out<-A[bin.index.f,bin.index.t]
        cat(i," bin index on ",chr,"f: ", bin.index.f,"t: ", bin.index.t, A[bin.index.f,bin.index.t], p.f[i,2], p.t[i,2], row.names(bins[(bins[,2]<=p.f[i,2] & bins[,3]>=p.f[i,3]),]), row.names(bins[(bins[,2]<=p.t[i,2] & bins[,3]>=p.t[i,3]),]), "\n",sep=" ")
        T<-A[((bin.index.f-eps):(bin.index.f+eps)), ((bin.index.t-eps):(bin.index.t+eps))]
        M.res<-M.res+T
      }
    }
   } 
}

# write.table(as.data.frame(as.matrix(M.res)), paste(filename, "_HiChIP_loops_IMR90_6689_pileup.csv", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
# M.res.mean<-M.res/mean(M.res)
# write.table(as.data.frame(as.matrix(M.res.mean)), paste(filename, "_HMGB2_eps_meanNorm.csv", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
M.res.median<-M.res/median(M.res)
# write.table(as.data.frame(as.matrix(M.res.median)), paste("HiChIP_loops_IMR90y_6689_g150k_10kb_medianNorm.csv", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
write.table(as.data.frame(as.matrix(M.res.median)), paste("SiHMGB2_HUV_loops_10kb_medianNorm_g150k.csv", sep=""), quote=F, row.names=F, col.names=F, sep="\t")