# Roll up to counts-draft
rm(list=ls())
#setwd("~/Desktop/DDIG_count_tables")#comment this out when running on cluster, navigate with command line
getwd()
library(reshape2)
#Part1: Menidia:####
Mb<-read.table("menidia_combined_rawcounts.txt")
Mb<-Mb[2:61280,]#just removing the top row with *=counts of unmapped reads
vars<-colsplit(Mb$GeneID_TranscriptID, "_s", c("GeneID","TranscriptID"))
Mb1<-cbind(Mb,vars)
Mb1<-Mb1[,c(1,55:56,2:54)]#just moving crap around so I can see it, if running as script can comment out
length(unique(Mb1$GeneID))#how many genes are there?

samples<-colnames(Mb1, do.NULL = TRUE)
samples<-samples[c(4:56)]
print(samples)
#aggregate and write to file####
Menidia_RNASeq_by_genes<-aggregate(cbind(Mb.CTMax.2.t0.1, Mb.CTMax.2.t0.2,  Mb.CTMax.2.t0.3,  
                      Mb.CTMax.2.t0.4,  Mb.CTMax.2.t0.5, Mb.CTMax.2.t0.6,  
                      Mb.CTMax.2.t60.1, Mb.CTMax.2.t60.2, Mb.CTMax.2.t60.3, 
                      Mb.CTMax.2.t60.4,Mb.CTMax.2.t60.5, Mb.HC.t0.1,       
                      Mb.HC.t0.2,       Mb.HC.t0.5,       Mb.HC.t0.6,      
                      Mb.HC.t0.twelve,  Mb.HC.t60.1,      Mb.HC.t60.2,      
                      Mb.HC.t60.3,      Mb.HC.t60.4,     Mb.HC.t60.5,      
                      Mb.Tenv.t0.1,     Mb.Tenv.t0.2,     Mb.Tenv.t0.3,     
                      Mb.Tenv.t0.4,    Mb.Tenv.t0.5,     Mb.Tenv.t0.6,     
                      Mb.Tenv.t60.1,    Mb.Tenv.t60.2,    Mb.Tenv.t60.4,   
                      Mb.Tenv.t60.5,    Mb.Tenv.t60.6,    Mb.Thigh.t0.1,    
                      Mb.Thigh.t0.2,    Mb.Thigh.t0.3,   Mb.Thigh.t0.4,    
                      Mb.Thigh.t0.6,    Mb.Thigh.t0.8,    Mb.Thigh.t60.1,   
                      Mb.Thigh.t60.2,  Mb.Thigh.t60.3,   Mb.Thigh.t60.4,   
                      Mb.Thigh.t60.8,   Mb.Tlow.t0.1,     Mb.Tlow.t0.2,    
                      Mb.Tlow.t0.3,     Mb.Tlow.t0.4,     
                      Mb.Tlow.t0.5,     Mb.Tlow.t60.2,    Mb.Tlow.t60.3,   
                      Mb.Tlow.t60.4,    Mb.Tlow.t60.6,    Mb.Tlow.t60.7) ~ GeneID, Mb1, sum)
####
#to do later-fix for loop or look into switching to matrices/apply functions-
#just did by brute force for now bc really pressed for time
#for (i in samples) {
  #new.name <- paste(i)
  #print(new.name)
  #test<-aggregate(new.name~GeneID, Mb1,sum)
  #sample<-rename(test,c("GeneID"="GeneID","x"=new.name))
  #base1<-merge(base1,sample,by="GeneID")
#}
write.table(Menidia_RNASeq_by_genes,file="Menidia_RNASeq_rolledupto_genes.txt",sep="\t",quote=FALSE)

#Part2: Delta smelt:####
rm(list=ls())#just in case, start fresh
#setwd("~/Desktop/DDIG_count_tables")#comment this out when running on cluster, navigate with command line
getwd()
library(reshape2)
DS<-read.table("DS_combined_rawcounts.txt")
DS<-DS[2:60738,]#just removing the top row with *=counts of unmapped reads
DS_key_table<-read.table("delta_isotigs_plus_trinity_assembly_correspondence_table.txt", header=TRUE)
DS_comb<-merge(DS,DS_key_table, by="GeneID_TranscriptID", all.x=TRUE)
DS_comb<-DS_comb[,c(1,45,2:44)]

length(unique(DS_comb$GeneID))#how many genes are there?
samples<-colnames(DS_comb, do.NULL = TRUE)
samples<-samples[c(3:45)]
print(samples)

#aggregate and write to file####
Deltasmelt_RNASeq_by_genes<-aggregate(cbind(Ds.HC.t0.1,  Ds.HC.t0.2,  Ds.HC.t0.3,  Ds.HC.t0.4,  Ds.HC.t0.5,  Ds.HC.t60.1,
                                         Ds.HC.t60.2, Ds.HC.t60.3, Ds.HC.t60.4, Ds.HC.t60.5, Ds.T2.t0.1,  Ds.T2.t0.2, 
                                         Ds.T2.t0.3,  Ds.T2.t0.4,  Ds.T2.t0.5,  Ds.T2.t0.6,  Ds.T2.t60.1, Ds.T2.t60.2,
                                         Ds.T2.t60.3, Ds.T2.t60.4, Ds.T2.t60.5, Ds.T3.t0.1,  Ds.T3.t0.2,  Ds.T3.t0.3, 
                                         Ds.T3.t0.4,  Ds.T3.t0.5,  Ds.T3.t0.6,  Ds.T3.t60.1, Ds.T3.t60.2, Ds.T3.t60.3,
                                         Ds.T3.t60.4, Ds.T3.t60.5, Ds.T4.t0.1,  Ds.T4.t0.2,  Ds.T4.t0.3,  Ds.T4.t0.4, 
                                         Ds.T4.t0.5,  Ds.T4.t0.6,  Ds.T4.t60.1, Ds.T4.t60.2, Ds.T4.t60.3, Ds.T4.t60.5,
                                         Ds.T4.t60.6) ~ GeneID, DS_comb, sum)
write.table(Deltasmelt_RNASeq_by_genes,file="Deltasmelt_RNASeq_rolledupto_genes.txt",sep="\t",quote=FALSE)

