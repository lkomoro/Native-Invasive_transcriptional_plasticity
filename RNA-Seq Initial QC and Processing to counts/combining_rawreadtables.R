# combining tables
#setwd("~/Desktop/test.DDIG_counts")#comment this out when running on cluster, navigate with command line
library(plyr)
#read in data
file.list <- list.files(pattern = "counttable.txt")
sample_ID<-strsplit(file.list, '_')
f <- function(x) x[[1]]
t1 <- unlist(lapply(sample_ID, f))
#t <- unlist(sample_ID)
#maxval<-length(file.list)
#t1 <- t[seq(1, maxval, by=2)]
#t1 <- t[c(1, 3, 5, 7)]  # was the clunky way to do same thing
t2 <- unique(t1)   # unique names-NB if there are samples that have common #s like HC-1 and HC-12, they won't
#be considered unique, so may need to reformat file names accordingly (like I just had one HC-12, so I replaced it with HC-twelve for this step)

for (i in 1:length(t2)){
  idx <- grep(t2[i], file.list)
  file1 <- file.list[idx[1]]
  file2 <- file.list[idx[2]]
print(file1)
print(file2)#this goes to nohup.out (if specified in command line),
  #just as a double check that the right files are being read in together
  data_PE<-read.table(file.list[idx[1]], header=FALSE)
  data_S<-read.table(file.list[idx[2]], header=FALSE)
#rename columns
data_PE<-rename(data_PE,c("V1"="GeneID_TranscriptID","V2"="Length_PE","V3"="Mapped_Counts_PE","V4"="Unmapped_Counts_PE"))
data_S<-rename(data_S,c("V1"="GeneID_TranscriptID","V2"="Length_S","V3"="Mapped_Counts_S","V4"="Unmapped_Counts_S"))
#Merge PE with singleton files, get rid of extra columns
data_comb<-merge(data_PE,data_S, all=TRUE)
data_comb_short<-data_comb[,c(1,3,6)]
#Add together,get rid of individual columns:
data_comb_short$Mapped_Counts_PE_S<-data_comb_short$Mapped_Counts_PE+data_comb_short$Mapped_Counts_S
data_comb_short1<-data_comb_short[,c(1,4)]
#write new file
write.csv(data_comb_short1, file = paste(strsplit(file.list[idx[1]], "_")[[1]][1],"PEandS_rawcounts", ".csv", sep = ""))
#can change this to be a write .txt if prefer, send it to its own new directory, etc.

#generate some summary stats of contribution of PE vs singletons
PE_mapped<-sum(data_comb_short$Mapped_Counts_PE)
S_mapped<-sum(data_comb_short$Mapped_Counts_S)
Total_mapped<-sum(data_comb_short$Mapped_Counts_PE_S)
#str(data_comb_short)
Proportion_S_mapped_of_total<-S_mapped/Total_mapped
Sample.name<-strsplit(file.list[idx[1]], "_")[[1]][1]
sum_stats<-cbind(Sample.name,PE_mapped, S_mapped, Total_mapped, Proportion_S_mapped_of_total)
write.table(sum_stats,file="mapping_sum_stats.txt",append=TRUE,sep="\t",col.names = FALSE,quote=FALSE)#hacky way to just get a table summarizing all the coverage of PE and S across samples
}
