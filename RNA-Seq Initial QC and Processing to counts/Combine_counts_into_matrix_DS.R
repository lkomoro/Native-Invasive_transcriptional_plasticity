rm(list=ls())
#setwd("~/Desktop/test.DDIG_counts")#comment this out when running on cluster, navigate with command line
#note-need to run seperately for DS and Menidia since different transcriptomes
library(plyr)
file.list <- list.files(pattern = "_rawcounts.csv")

base<-read.csv(file.list[1])
base1<-base[,2]
base1<-data.frame(base1)
base1<-rename(base1,c("base1"="GeneID_TranscriptID"))#this nonsense is just bc it changes it to a list when you drop to one variable, and then you have to change the column name back

for (i in file.list) {
  sample<-read.csv(i)
  new.name <- paste(strsplit(i, "PE")[[1]][1])
sample<-sample[,2:3]
sample<-rename(sample,c("GeneID_TranscriptID"="GeneID_TranscriptID","Mapped_Counts_PE_S"=new.name))
base1<-merge(base1,sample,by="GeneID_TranscriptID")
}

write.table(base1,file="DS_combined_rawcounts.txt",sep="\t",quote=FALSE)#hacky way to just get a table summarizing all the coverage of PE and S across samples

