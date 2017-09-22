library(data.table)
library(tools)
library(GenomicRanges)
library(stringr)

# File with exons hg19
ccds<-fread("/home/milos/R/recursive_splicing/CCDS.current.txt", data.table=FALSE)#, colClasses='character')

### Used when script is called from the command line ###
args<-commandArgs(TRUE)
file_f<-args[1]

# file_f<-"/home/milos/R/recursive_splicing/RS_results/three_prime_to_cis_intron_Adrenal_gland.txt"

f<-fread(file_f, data.table=FALSE)
f<-f[,c(1,3)]     # f<-f[,c(1,2)]
names(f)<-c("ccds_id","offset")
f.sorted<-unique(f[order(f[,1]),])

### Correct the offset from RS pipeline for 200nts (so that gene/exon coordinates comply with UCSC genome coords ***)
f.sorted$offset<-as.numeric(f.sorted$offset)
#####################################f.sorted$offset<-f.sorted$offset-200

f.sorted<-f.sorted[which(f.sorted[,2]>0),]

filename<-strsplit(file_path_sans_ext(basename(file_f)),"three_prime_to_cis_intron_")[[1]][[2]]
# filename<-"test"
### Remove "Withdrawn" ids and keep the "Public" ones and keep the "relevant columns" only

ccds<-ccds[,c(5,1,8,9,7,3,10)]
names(ccds)[[2]]<-"chr"
ccds.sorted<-ccds[order(ccds[,1]),]


### Read the second file needed for the final table ###
# file_f2<-"/home/milos/R/recursive_splicing/testing/additional_2_exon_cols_Adrenal_gland.csv"

### Used when script is called from the command line ###
file_f2<-args[2]

f2<-fread(file_f2, data.table=FALSE)
f2.sorted<-f2[order(f2[,1]),]
names(f2.sorted)[1]<-c("ccds_id")
m<-merge(f.sorted, f2.sorted, by="ccds_id", all.x=TRUE)

### Adjust the "Three prime to cis intron" ###

# a) Split the first column
f1<-data.frame(matrix(unlist(strsplit(m[,1], '_')),ncol=3, byrow=T))[,c(1,3)]  ### Cool option to convert a list to df!!!
f1<-cbind(f1,m[,2:4])
names(f1)[1:2]<-c("ccds_id","exon")

# b) Order by first column and keep unique rows
f1.sorted<-unique(f1[order(f1[,1]),])
f1.sorted[,2]<-as.integer(as.character(f1.sorted[[2]]))
merged<-merge(f1.sorted, ccds.sorted, by="ccds_id", all.x=TRUE)
merged<-merged[-which(is.na(merged[,7])),]              
merged<-merged[-which(merged[,7]=='-'),]

merged[,7]<-as.integer(as.matrix(merged[,7]))
merged[,8]<-as.integer(as.matrix(merged[,8]))

merged$offset<-merged$offset-200

#########################################################
### Correct exon values for genes on the minus strand ###
#########################################################
tmp<-merged[merged$cds_strand=='-',]
tmp_merged<-merged

### Generate columns with exon, next intron and exon coordinates ###

ex1<-NULL
ex2<-NULL
next.intr.x1<-NULL
next.intr.x2<-NULL
next.ex1<-NULL
next.ex2<-NULL

for (i in 1:nrow(merged)){
  # Find the corresponding exon and get its start and end
  if(merged[i,9]=='+'){ ### In case of the plus strand ###
    (t<-unlist(strsplit(merged[i,11],";"))[merged[i,2]])
    (t<-(strsplit(t, "-")))

    ### Get coords from next exon
    (t_next<-unlist(strsplit(merged[i,11],";"))[(merged[i,2]+1)])
    (t_next<-(strsplit(t_next,"-")))
    (next.ex1<-rbind(next.ex1, as.integer(t_next[[1]][[1]])))
    if(is.na(next.ex1[i])){#,rbind(next.ex2,"NA"),rbind(next.ex2, as.integer(t_next[[1]][[2]])))
      next.ex2<-rbind(next.ex2,"NA")
    }else{
      next.ex2<-rbind(next.ex2, as.integer(t_next[[1]][[2]]))
    }

    # cat(merged[i,2], t, "\n", sep=" ")
    ex1<-rbind(ex1, as.integer(t[[1]][1]))
    ex2<-rbind(ex2, as.integer(t[[1]][2]))
    next.intr.x1<-rbind(next.intr.x1,(as.integer(t[[1]][2])))
    t1<-unlist(strsplit(merged[i,11],";"))[(merged[i,2])+1]
    if(is.na(t1)){ next.intr.x2<-rbind(next.intr.x2,-1)} #as.integer(t[[1]][2]))}
    else{
        t1<-(strsplit(t1, "-"))
        next.intr.x2<-rbind(next.intr.x2,as.integer(t1[[1]][1]))
        # cat(merged[(i+1),2], t[[1]][1],t[[1]][2], "\n", sep=" ")
      }
   }else{    ### Minus strand ###
             ####################
    # k is tmp var
    (k<-str_count(merged[i,11], ';') - merged[i,2] + 2 )
    (t<-unlist(strsplit(merged[i,11],";"))[k])
    (t<-(strsplit(t, "-")))

    ### Get coords from next exon
    if (k>1){
      (t_next<-unlist(strsplit(merged[i,11],";"))[k-1])
      (t_next<-(strsplit(t_next,"-")))
      (next.ex1<-rbind(next.ex1, as.integer(t_next[[1]][[1]])))
      if(is.na(next.ex1[i]) || length(t_next)==0){#,rbind(next.ex2,"NA"),rbind(next.ex2, as.integer(t_next[[1]][[2]])))
        next.ex2<-rbind(next.ex2,-1)
      }else{ (next.ex2<-rbind(next.ex2, as.integer(t_next[[1]][[2]]))) }
      
      if(!is.na(t_next)){
        next.intr.x1<-rbind(next.intr.x1,(as.integer(t_next[[1]][2])))
      }else{next.intr.x1<-rbind(next.intr.x1, -1) }
    }else{
      (next.ex2<-rbind(next.ex2, -1))
      (next.ex1<-rbind(next.ex1, -1))
      next.intr.x1<-rbind(next.intr.x1, -1)
      }

    ex1<-rbind(ex1, as.integer(t[[1]][1]))
    ex2<-rbind(ex2, as.integer(t[[1]][2]))

    if (i<nrow(merged)){
      t1<-unlist(strsplit(merged[i,11],";"))[k-1]
      if(length(t1)==0){#is.na(t1)){
        next.intr.x2<-rbind(next.intr.x2, "NA")#x2,as.integer(t[[1]][2]))
      }
      else{
        t1<-(strsplit(t1, "-"))
        next.intr.x2<-rbind(next.intr.x2,as.integer(t[[1]][1]))
      }
    }else{
      next.intr.x2<-rbind(next.intr.x2,as.integer(t[[1]][1]))
    }
  }
}

next.ex1<-as.data.frame(as.matrix(next.ex1))
next.ex2<-as.data.frame(as.matrix(next.ex2))

ex1<-as.data.frame(as.matrix(ex1))
ex2<-as.data.frame(as.matrix(ex2))

next.intr.x1<-as.data.frame(as.matrix(next.intr.x1))
next.intr.x2<-as.data.frame(as.matrix(next.intr.x2))
merged1<-cbind(merged[,1:10],ex1,ex2)

merged1$cds_strand<-as.character(merged1$cds_strand)
merged1<-cbind(merged1,rep(0,nrow(merged1)))
names(merged1)[11:13]<-c("e_x1","e_x2","rs_pos")
merged1[,13]<-as.integer(merged1[,13])
merged1$offset<-as.numeric(merged1$offset)

merged1[merged1$cds_strand=="+",13]<-merged1[merged1$cds_strand=="+",7] + as.integer(merged1[merged1$cds_strand=="+",3]) #7
merged1[merged1$cds_strand=="-",13]<-merged1[merged1$cds_strand=="-",8] - as.integer(merged1[merged1$cds_strand=="-",3]) + 1 #8

merged1<-cbind(merged1, between(as.integer(merged1$rs_pos), merged1$cds_from,merged1$cds_to))
names(merged1)[14]<-c("within_gene")

merged1<-cbind(merged1,next.intr.x1,next.intr.x2)
names(merged1)[15:16]<-c("consec_intron_x1", "consec_intron_x2")
put.to.zero.ind<-which(merged1[,15]==merged1[,16])

merged1[[15]]<-suppressWarnings(as.numeric(as.character(merged1[[15]])))
merged1[[16]]<-suppressWarnings(as.numeric(as.character(merged1[[16]])))

merged1[merged1$cds_strand=='+',15]<-merged1[merged1$cds_strand=='+',15]+1
merged1[merged1$cds_strand=='+',16]<-merged1[merged1$cds_strand=='+',16]-1

merged1[merged1$cds_strand=='-',15]<-merged1[merged1$cds_strand=='-',15]+1
merged1[merged1$cds_strand=='-',16]<-merged1[merged1$cds_strand=='-',16]-1

merged1<-cbind(merged1, between(merged1$rs_pos, merged1$consec_intron_x1, merged1$consec_intron_x2))
names(merged1)[17]<-"within_consec_intron"

### Add intron sizes and coordinates
intron_x1<-NULL
intron_x2<-NULL
next_e_x1<-NULL
next_e_x2<-NULL
exons<-merged[,11]
exons1<-strsplit(exons, ';') # splitted exons
intron.lengths<-NULL
for(i in 1:nrow(merged1)){
  e<-exons1[[i]] ### take one row with exons (current)
  if(length(e)>1){
    e1<-data.frame(matrix(unlist(strsplit(e,"-")), ncol=2, byrow=T)) ### split even more so that intron regions can be constructed ---                ***VERY COOL OPTION***
    introns<-cbind(as.data.frame(as.integer(as.matrix(e1[1:(nrow(e1)-1),2])))+1, as.data.frame(as.integer(as.matrix(e1[2:nrow(e1),1])))-1)
    names(introns)<-c("intron_x1", "intron_x2")
    intron.chk<-between(as.integer(merged1$rs_pos[i]), introns[,1],introns[,2])
    if(any(intron.chk)){
      intron.ind<-which(intron.chk==TRUE)
      intron.len<-introns[intron.ind,2]-introns[intron.ind,1]
      intron.lengths<-rbind(intron.lengths, intron.len)
      # Get intron coordinates
      intron_x1<-rbind(intron_x1, introns[intron.ind,1])
      intron_x2<-rbind(intron_x2, introns[intron.ind,2])
      next_e_x1<-rbind(next_e_x1, introns[intron.ind,2]+1)
      next_e_x2<-rbind(next_e_x2, introns[intron.ind+1,1]-1)
    }else{
      intron.lengths<-rbind(intron.lengths, 0)
      intron_x1<-rbind(intron_x1, 0)
      intron_x2<-rbind(intron_x2, 0)
      next_e_x1<-rbind(next_e_x1, 0)
      next_e_x2<-rbind(next_e_x2, 0)
    }
  }else{
    intron.lengths<-rbind(intron.lengths, 0)
    intron_x1<-rbind(intron_x1, 0)
    intron_x2<-rbind(intron_x2, 0)
    next_e_x1<-rbind(next_e_x1, 0)
    next_e_x2<-rbind(next_e_x2, 0)
  }
}

merged1<-cbind(merged1[,1:17], intron.lengths)
names(merged1)[18]<-"intron_len"

merged1<-cbind(merged1, intron_x1, intron_x2)
names(merged1)[19:20]<-c("intron_x1","intron_x2")

merged1<-cbind(merged1, next.ex1,next.ex2)
names(merged1)[21:22]<-c("next_exon_x1","next_exon_x2")

merged2<-merged1[-which(is.na(merged1[,4])),]
merged2<-merged2[-which(merged2$intron_len==0),]

# Add columns representing rs-exon distances 
merged2<-cbind(merged2,abs(merged2$rs_pos-merged2$e_x2), abs(merged2$next_exon_x1-merged2$rs_pos))
names(merged2)[(ncol(merged2)-1):ncol(merged2)]<-c("dist_5prime_exon","dist_3prime_exon")

#write.table(merged2, paste(filename, "_full_table.csv", sep=""), quote=F, row.names=F, col.names=T, sep="\t")

##############################################
### ADJUST THE FILTERING CRITERIA HERE !!! ###
##############################################

tpi_filter<-4     # Three prime to intron column
icic_filter<-1    # Inferred cis intron connections column

# Arbitrary filter criteria
indexes<-which(merged2[,4]>=tpi_filter & merged2[,5]>=icic_filter & merged2$within_gene==TRUE)
result<-merged2[indexes,]
nrow(result)
# If the downstream exon's coord x2 is NA
result<-result[!is.na(result$next_exon_x2),]
nrow(result)
# Sum of distances from RS site to upstream and downstream exon >=2500
result<-result[abs(result$dist_3prime_exon)+abs(result$dist_5prime_exon)>=2500,]
nrow(result)
# Distance >=500 from both upstream and downstream exons
result<-result[abs(result$dist_3prime_exon)>=500 & abs(result$dist_5prime_exon)>=500,]
nrow(result)
# Remove circles (everything that is bound from 3'exon to anything that is upstream)
ind.plus<-which(result$rs_pos>result$e_x2 & result$cds_strand=="+")
ind.minus<-which(result$rs_pos<result$e_x2 & result$cds_strand=="-")
ind.rows<-unique(sort(c(ind.plus, ind.minus)))
result<-result[ind.rows,]
nrow(result)

### Used when script is called from the command line ###
file_f3<-args[3]

### Filter even more using 2 more columns from new_results_*.xls files ###
# file_f3<-"/home/milos/R/recursive_splicing/testing/additional_2_exon_cols1_Adrenal_gland.csv"

f3<-fread(file_f3, data.table=FALSE)
f3.sorted<-f3[order(f3[,1]),]
names(f3.sorted)[1]<-c("ccds_id")
m1<-merge(f.sorted, f3.sorted, by=c("ccds_id"), all.x=TRUE)
# m1<-merge(f.sorted, f3.sorted, by=c("ccds_id","exon"), all.x=TRUE) ## no

f3<-data.frame(matrix(unlist(strsplit(m1[,1], '_')),ncol=3, byrow=T))[,c(1,3)]  ### Cool option to convert a list to df!!!
f3<-cbind(f3,m1[,2:4])
names(f3)[1:2]<-c("ccds_id","exon")
f3<-f3[,-3]

f3.sorted<-unique(f3[order(f3[,1]),])
f3.sorted[,2]<-as.integer(as.character(f3.sorted[[2]]))

### Extend RS sites ###
res.garlic<-cbind(result,result$rs_pos,result$rs_pos)
res.garlic[res.garlic$cds_strand=='+',25]<-res.garlic[res.garlic$cds_strand=='+',13]-45
res.garlic[res.garlic$cds_strand=='+',26]<-res.garlic[res.garlic$cds_strand=='+',13]+5
res.garlic[res.garlic$cds_strand=='-',25]<-res.garlic[res.garlic$cds_strand=='-',13]-5
res.garlic[res.garlic$cds_strand=='-',26]<-res.garlic[res.garlic$cds_strand=='-',13]+45

colnames(res.garlic)[25:26]<-c("rs_pos_x1", "rs_pos_x2")

### Merge f3 df with the rest of the data ###
merged_final<-merge(res.garlic, f3.sorted, by=c("ccds_id","exon"), all.x=TRUE)
nrow(merged_final)
# write.table(merged_final, paste(filename, "_filtered_table_","tpi",tpi_filter,"_icic",icic_filter,"_final_DELME.csv", sep=""), quote=F, row.names=F, col.names=T, sep="\t")
merged_final<-merged_final[-which(merged_final$Observed_three_prime_cis_connections<100),]
nrow(merged_final)
merged_final<-merged_final[-which(merged_final$Five_prime_read_through<1),]
nrow(merged_final)

write.table(merged_final, paste(filename, "_filtered_table_","tpi",tpi_filter,"_icic",icic_filter,"_final_TEST.csv", sep=""), quote=F, row.names=F, col.names=T, sep="\t")

### This line is moved to rs4.R
write.table(unique(cbind(paste("chr", merged_final[[6]], sep=""), merged_final[,c(25,26)])), paste(filename, "_filtered_table_","tpi",tpi_filter,"_icic",icic_filter,"_unique_TEST.csv", sep=""), quote=F, row.names=F, col.names=T, sep="\t")

table(merged_final[merged_final$cds_strand=='-',17])
table(merged_final[merged_final$cds_strand=='+',17])