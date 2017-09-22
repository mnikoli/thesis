### Postprocessing step after using rs3_III.R script -- the purpose of the 
### script is to remove novel exons from the RS list and generate input for GARLIC

library(GenomicRanges)
library(data.table)
library(tools)
                              #####################################################
correct<-1                    ### ADJUST THIS PARAMETER (AND THE RC_PATH) ! ! ! ### 
                              #####################################################
### Use this line if there is only one file in a folder
# file1<-fread("/home/milos/R/Adrenal_gland_filtered_table_tpi4_icic1_final.csv", data.table=FALSE)
# name<-"Adrenal_gland"
# rc_path<-"/home/milos/R/recursive_splicing/testing/newest"                              #### Path to totalRNAseq samples
rc_path<-"/home/milos/R/recursive_splicing/testing/samples_for_polyA_correction"      #### Path to RNAseq samples that should be corrected using polyA RNAseq samples
# rc_path<-"/home/milos/R/recursive_splicing/testing/missing_samples"                   #### Tmp path when some samples are missing or are repeated
filenames<-list.files(rc_path, pattern="*_final_TEST.csv", full.name=TRUE) 
lrc<-lapply(filenames,sep='\t', header=TRUE, read.csv) # read each file as a list of elements

for (i in 1:length(lrc)){
  name<-file_path_sans_ext(basename(filenames[i]))
  name<-strsplit(name, "_filtered_table_tpi4_icic1_final_TEST")[[1]]
  file1<-lrc[i]
  cat(name, "\n")
  if (!correct){
    result<-as.data.frame(file1)
  }else{
    libsizes<-fread("/home/milos/R/recursive_splicing/manuscript/lib_sizes.csv", data.table=FALSE)  # TGF and TNF samples (2 in total) are excluded
    
    # polyA_name<-paste(name, "_filtered_table_tpi4_icic1_III_polyA.csv", sep='')
    polyA_name<-paste(name, "_filtered_table_tpi4_icic1_final_polyA.csv", sep='')
    file2_polyA<-fread(paste("/home/milos/R/recursive_splicing/polyA_results/",polyA_name, sep=''), data.table=FALSE)

    file1<-as.data.frame(file1)
    f1<-file1[order(file1$rs_pos),]
    f2<-file2_polyA[order(file2_polyA$rs_pos),]

    ### Arbitrary cutoff --- take only those RS sites having at least 40 reads (10x more than we require for RS sites) --- we are now looking for novel exons
    f2<-f2[f2$Three_prime_to_intron>=40,]
    head(f1)
    head(f2)

    ### Get indexes from same rows
    same_rows1<-which(f1$ccds_id %in% f2$ccds_id & f1$rs_pos %in% f2$rs_pos)
    same_rows2<-which(f2$ccds_id %in% f1$ccds_id & f2$rs_pos %in% f1$rs_pos)
    
    sub1<-f1[same_rows1,]
    sub2<-f2[same_rows2,]

    n1<-libsizes[libsizes[,1]==name,2]        # 274034110  # lib from f1 adrenal gland
    n2<-libsizes[libsizes[,1]==name,3]        # 112742654  # lib from f2 adrenal gland

    fact<-n1/n2
    fact

    # small fix (merging of same rs sites on differently reported (consecutive) exons) - before the final filtering
    sub1.index<-which(diff(sub1$offset)==0)
    sub2.index<-which(diff(sub2$offset)==0)
    
    if (length(sub1.index)>0){
      sub1[sub1.index,4]<-sub1[sub1.index,4]+sub1[sub1.index+1,4]
      sub1<-sub1[-(sub1.index+1),]
    }
    if (length(sub2.index)>0){
      sub2[sub2.index,4]<-sub2[sub2.index,4]+sub2[sub2.index+1,4]
      sub2<-sub2[-(sub2.index+1),]
    }
    
    table(sub1$ccds_id==sub2$ccds_id & sub2$Three_prime_to_intron*fact-sub1$Three_prime_to_intron>=0.1*sub1$Three_prime_to_intron & sub1$rs_pos==sub2$rs_pos)
    for_removal<-which(sub1$ccds_id==sub2$ccds_id & sub2$Three_prime_to_intron*fact-sub1$Three_prime_to_intron>=1.1*sub1$Three_prime_to_intron & sub1$rs_pos==sub2$rs_pos)
    # RS sites without novel exons
    # result<-f1[-for_removal,]

    ### PolyA removed --- output list ###
    write.table(f1[for_removal,], paste(name, "_polyA_removed.csv", sep=""), quote=F, row.names=F, sep="\t")
    result<-f1[-for_removal,]
  }
  
  ### LOGO3 --- output list ###
  out<-cbind(paste("chr", result[[6]], sep=""),result[,c(25,26)])#,result[,c(1:5,7:24,27:28)])
  colnames(out)[1]<-"chr"
  write.table(unique(out), paste(name, "_logo2.csv", sep=""), quote=F, col.names=F, row.names=F, sep="\t")
  
  ### Merge overlapping intervals not to have a bias when tested with GARLIC ###
  for_merging<-which(abs(diff(result$rs_pos))<=10)
  nrow(result)
  
  if(length(for_merging)>0){
    for (i in 1:length(for_merging)){
    vmax<-max(result[for_merging[i],26], result[(for_merging+1)[i],26])
    vmin<-min(result[for_merging[i],25], result[(for_merging+1)[i],25])
    # cat("(vmin,vmax)=",vmin,vmax, sep="")
    result[for_merging[i],25]<-vmin
    result[for_merging[i],26]<-vmax
    }
  
    result<-result[-(for_merging+1),]
  }
  nrow(result)
  out.garlic<-cbind(paste("chr", result[[6]], sep=""), result[,c(25,26)])
  out.garlic<-out.garlic[complete.cases(out.garlic),]                         #### COOL FUNCTION TO SKIP NAs IN ANY COLUMN ####
  write.table(out.garlic, paste(name, "_RS2.csv", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
}

