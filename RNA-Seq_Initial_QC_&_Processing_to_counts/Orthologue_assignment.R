#Intro:####

library(tidyr)
library(dplyr)
library(reshape2)
library(here)
here()

#list.files()

#Read in Data:####
orthologues<-read.table(here("./DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthofinder_analyses/analysis_2020/Orthologues_LMKworked.txt"), header=T)
Mb_ortho<-orthologues[,c(1:3)]
Ds_ortho<-orthologues[,c(1:2,4)]
#seperate list of transcripts

#Menidia Part 1: ####
#use the summary info from the Statistics_PerSpecies.tsv output file to determine how many columns are needed
#see that Menidia has only one orthogroup with 21-50 genes (and sorting in sheet actually seems like the max is 20?)
Mb_ortho1<-Mb_ortho  %>%  separate(Menidia_pep, paste("col1",1:21,sep="."),", ") 
Mb_ortho_long <- gather(Mb_ortho1, column, ORF_transcriptID, col1.1:col1.21) #convert from wide to long format
CR<-sum(is.na (Mb_ortho_long$ORF_transcriptID))#how many rows should be removed below
Mb_ortho_long1<-Mb_ortho_long[complete.cases(Mb_ortho_long[ , 4]),]
Keep<-length(Mb_ortho_long1$ORF_transcriptID)
CR+Keep #(this is just a check: should add up to equal the original number of cases in the long dataframe above)

#Transdecoder adds a ".pX" (e.g., .p1 or .p2, etc.) to the end of the transcript IDs, so need to rectify before merging with counts table
vars1<-colsplit(Mb_ortho_long1$ORF_transcriptID, "_p", c("GeneID_TranscriptID","transdecoder_ORF")) #NB I changed the ".p" to be "_p" in the import file because R is wildcarding "."
Mb_ortho_long1<-cbind(Mb_ortho_long1,vars1)

length(unique(Mb_ortho_long1$GeneID_TranscriptID))#how many unique transcript IDs end with (which should match ref transcriptome)?
length(unique(Mb_ortho_long1$ORF_transcriptID))#how many unique transcript IDs start with, including transdecoder ORF suffix?
length(Mb_ortho_long1$ORF_transcriptID) #how many cases in column total?
#Good, these are now all the same, so the single_best_orf function in transcdecoder has fixed the problem of multiple results for the same transcript!
Mb_dups<-Mb_ortho_long1[duplicated(Mb_ortho_long1$GeneID_TranscriptID)|duplicated(Mb_ortho_long1$GeneID_TranscriptID, fromLast=TRUE),] #(this should now be 0)

write.csv(Mb_ortho_long1,here("./DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthofinder_analyses/analysis_2020/Mb_ortho_key_by_transcript.csv"),row.names = FALSE)

#Menidia Part II: pull in counts per transcript table ####
Mb<-read.table(here("DDIG_data&analyses/DDIG_count_tables and differential expression analysis/Filtered_transcriptomes_rm2reads/menidia_combined_rawcounts.txt"))#using the count tables with secondary reads removed
Mb<-Mb[2:61280,]#just removing the top row with *=counts of unmapped reads
vars<-colsplit(Mb$GeneID_TranscriptID, "_s", c("GeneID","TranscriptID"))
Mb1<-cbind(Mb,vars)
Mb1<-Mb1[,c(1,55:56,2:54)]#just moving crap around so I can see it, if running as script can comment out
length(unique(Mb1$GeneID))#how many genes are there? #26833

#So there are 61279 transcripts in the Mb1 df, which is the raw counts for all the transcripts in the Menidia reference
#There are 22104 cases in the Mb_ortho_long1 df, which is the transcripts that have been assigned an orthogroup/orthologue that is shared between the two species
#merge:
Menidia_sharedortholog_rawcounts<-merge(Mb_ortho_long1,Mb1)#Keeps only entries that are shared between the datasets. The number of cases should match the original number in Mb_ortho_long1 since should just be a subset of Mb1
#Menidia_sharedortholog_rawcounts_ally<-merge(Mb_ortho_long1,Mb1,all.y=TRUE) #keep all the entries that have transcripts in the reference

#check that have all unique transcripts from transdecoder input
length(unique(Menidia_sharedortholog_rawcounts$GeneID_TranscriptID))#how many unique transcript IDs end with (which should match ref transcriptome)?
length(unique(Menidia_sharedortholog_rawcounts$ORF_transcriptID))#how many unique transcript IDs start with, including transdecoder ORF suffix?
Mb_dups<-Menidia_sharedortholog_rawcounts[duplicated(Menidia_sharedortholog_rawcounts$GeneID_TranscriptID)|duplicated(Menidia_sharedortholog_rawcounts$GeneID_TranscriptID, fromLast=TRUE),] #should be zero

SO.Mb<-Menidia_sharedortholog_rawcounts[,c(1:3,9:61)]
write.csv(SO.Mb,here("./DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthofinder_analyses/analysis_2020/Mb_orthogroup_rawcounts.csv"),row.names = FALSE)

#Menidia Part III: Bring in Annotation ####
Mb_annot<-read.csv(here("DDIG_data&analyses//DDIG_annotation&Functional_analyses/menidia_annotation_files/M_beryllina_fish_B2G_top_blasts.csv")) #note this is just for sequence description, not GO IDs, see section below for GO IDs
Mb_annot<-rename(Mb_annot, GeneID_TranscriptID=Sequence.name)
Mb_SO_annot<-merge(Mb_annot, SO.Mb, all.y=TRUE)
Mb_SO_annot1<-Mb_SO_annot[,c(1:2,11:65)]

#Delta smelt Part I: ####
Ds_ortho1<-Ds_ortho  %>%  separate(Deltasmelt_pep, paste("col1",1:31,sep="."),", ") 
sum(is.na(Ds_ortho1$col1.31))#just checking where the cutoff is here

Ds_ortho_long <- gather(Ds_ortho1, column, ORF_transcriptID, col1.1:col1.31) #convert from wide to long format
CR<-sum(is.na (Ds_ortho_long$ORF_transcriptID))#how many rows should be removed below
Ds_ortho_long1<-Ds_ortho_long[complete.cases(Ds_ortho_long[ , 4]),]
Keep<-length(Ds_ortho_long1$ORF_transcriptID)
CR+Keep #(this is just a check: should add up to equal the original number of cases in the long dataframe above)

#Transdecoder adds a ".pX" (e.g., .p1 or .p2, etc.) to the end of the transcript IDs, so need to rectify before merging with counts table
vars1<-colsplit(Ds_ortho_long1$ORF_transcriptID, "_p", c("GeneID_TranscriptID","transdecoder_ORF"))#NB I changed the ".p" to be "_p" in the import file because R is wildcarding "."
Ds_ortho_long1<-cbind(Ds_ortho_long1,vars1)

length(unique(Ds_ortho_long1$GeneID_TranscriptID))#how many unique transcript IDs end with (which should match ref transcriptome)?
length(unique(Ds_ortho_long1$ORF_transcriptID))#how many unique transcript IDs start with, including transdecoder ORF suffix?
length(Ds_ortho_long1$ORF_transcriptID) #how many cases in column total?

#Good, these are now all the same, so the single_best_orf function in transcdecoder has fixed the problem of multiple results for the same transcript!

Ds_dups<-Ds_ortho_long1[duplicated(Ds_ortho_long1$GeneID_TranscriptID)|duplicated(Ds_ortho_long1$GeneID_TranscriptID, fromLast=TRUE),] #(this should now be 0)

write.csv(Ds_ortho_long1,here("./DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthofinder_analyses/analysis_2020/Ds_ortho_key_by_transcript.csv"),row.names = FALSE)

#Delta smelt Part II: pull in counts per transcript table ####
Ds<-read.table("DDIG_data&analyses/DDIG_count_tables and differential expression analysis/Filtered_transcriptomes_rm2reads/DS_combined_rawcounts.txt")
Ds<-Ds[2:60738,]#just removing the top row with *=counts of unmapped reads
length(unique(Ds$GeneID_TranscriptID))#There are 60737 transcripts in the Ds df, which is the raw counts for all the transcripts in the Delta smelt reference
length(unique(Ds_ortho_long1$GeneID_TranscriptID))#There are 16406 cases in the Ds_ortho_long1 df, which is the transcripts that have been assigned an orthogroup/orthologue that is shared between the two species

#merge:
Ds_sharedortholog_rawcounts<-merge(Ds_ortho_long1,Ds)#Keeps only entries that are shared between the datasets. The number of cases should match the original number in Ds_ortho_long1 since should just be a subset of Ds
#Ds_sharedortholog_rawcounts_ally<-merge(Ds_ortho_long1,Ds,all.y=TRUE) #keep all the entries that have transcripts in the reference

#check that have all unique transcripts from transdecoder input
length(unique(Ds_sharedortholog_rawcounts$GeneID_TranscriptID))#how many unique transcript IDs end with (which should match ref transcriptome)?
length(unique(Ds_sharedortholog_rawcounts$ORF_transcriptID))#how many unique transcript IDs start with, including transdecoder ORF suffix?
Ds_dups<-Ds_sharedortholog_rawcounts[duplicated(Ds_sharedortholog_rawcounts$GeneID_TranscriptID)|duplicated(Ds_sharedortholog_rawcounts$GeneID_TranscriptID, fromLast=TRUE),] #should be zero

#Merge Menidia and DS Shared Orthogroup files

SO.Ds<-Ds_sharedortholog_rawcounts[,c(1:3,7:49)]
write.csv(SO.Ds,here("/DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthofinder_analyses/analysis_2020/Ds_orthogroup_rawcounts.csv"),row.names = FALSE)

#Bring in Annotation:
Ds_annot<-read.table(here("DDIG_data&analyses/DDIG_annotation&Functional_analyses/delta_smelt_annotation_files/delta_merged_final.annot.txt"),header=FALSE, sep = "\t",fill = TRUE,na.strings=c("","NA"))#need fill=TRUE here because unequal columns in input file
CR1<-sum(is.na (Ds_annot$V3))#how many rows should be removed below (note: this makes it so GO terms are truncated, use this only for sequence descriptions! see below section for GO terms)
Ds_annot1<-Ds_annot[complete.cases(Ds_annot[ , 3]),]
Keep1<-length(Ds_annot1$V3)
CR1+Keep1#should be original observation number
Ds_annot1<-rename(Ds_annot1, GeneID_TranscriptID=V1)
Ds_annot1<-rename(Ds_annot1, Sequence.desc.=V3)
Ds_SO_annot<-merge(Ds_annot1, SO.Ds, all.y=TRUE)
Ds_SO_annot1<-Ds_SO_annot[,c(1,3,4:48)]

#-----------------
#Merge Ds and Mb Orthologues to compare annotations (gene descriptions)####
X<-Ds_SO_annot1[,1:4]
X$Species<-"Ds"
Y<-Mb_SO_annot1[,1:4]
Y$Species<-"Mb"
XY<-rbind(X,Y)
write.csv(XY, here("DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthofinder_analyses/analysis_2020/Annotations of Shared Orthogroups by Transcript.csv"), row.names = F)

#Roll up to Orthogroup & combining across species:####
#note this is the brute force code I used before to roll up to gene, need to change to be loop so transferrable
SO.Mb.rolledup<-aggregate(cbind(Mb.CTMax.2.t0.1, Mb.CTMax.2.t0.2,  Mb.CTMax.2.t0.3,  
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
                                         Mb.Tlow.t60.4,    Mb.Tlow.t60.6,    Mb.Tlow.t60.7) ~ Orthogroup, SO.Mb, sum)#change to orthologue here if want to use that instead
SO.Ds.rolledup<-aggregate(cbind(Ds.HC.t0.1,  Ds.HC.t0.2,  Ds.HC.t0.3,  Ds.HC.t0.4,  Ds.HC.t0.5,  Ds.HC.t60.1,
                                            Ds.HC.t60.2, Ds.HC.t60.3, Ds.HC.t60.4, Ds.HC.t60.5, Ds.T2.t0.1,  Ds.T2.t0.2, 
                                            Ds.T2.t0.3,  Ds.T2.t0.4,  Ds.T2.t0.5,  Ds.T2.t0.6,  Ds.T2.t60.1, Ds.T2.t60.2,
                                            Ds.T2.t60.3, Ds.T2.t60.4, Ds.T2.t60.5, Ds.T3.t0.1,  Ds.T3.t0.2,  Ds.T3.t0.3, 
                                            Ds.T3.t0.4,  Ds.T3.t0.5,  Ds.T3.t0.6,  Ds.T3.t60.1, Ds.T3.t60.2, Ds.T3.t60.3,
                                            Ds.T3.t60.4, Ds.T3.t60.5, Ds.T4.t0.1,  Ds.T4.t0.2,  Ds.T4.t0.3,  Ds.T4.t0.4, 
                                            Ds.T4.t0.5,  Ds.T4.t0.6,  Ds.T4.t60.1, Ds.T4.t60.2, Ds.T4.t60.3, Ds.T4.t60.5,
                                            Ds.T4.t60.6) ~ Orthogroup, SO.Ds, sum)#change to orthologue here if want to use that instead
#Combine 
Rolledup_comb<-merge(SO.Ds.rolledup,SO.Mb.rolledup) #has same number of observations as input (which are also the same)
write.csv(Rolledup_comb, "Orthogroup_count_speciescomb.csv", row.names = FALSE)
#Filter for including only orthologues or orthogroups that are shared between species
Rolledup_comb$Dscount<-rowSums(Rolledup_comb [2:44])
Rolledup_comb$Mbcount<-rowSums(Rolledup_comb [45:97])
Rolledup_combfilt<-subset(Rolledup_comb,Dscount>0 )
Rolledup_combfilt1<-subset(Rolledup_comb,Mbcount>0 ) #only reduces by 9, and justification for this unclear anyway, so use original for now


#Clean up and refresh: ####
#rm(list=ls()) #(uncomment to use)
library(tidyr)
library(dplyr)
library(reshape2)
annot.comb<-read.csv(here("./DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthofinder_analyses/analysis_2020/Annotations of Shared Orthogroups by Transcript.csv"))
counts.comb<-read.csv(here("./DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthofinder_analyses/analysis_2020/Orthogroup_count_speciescomb.csv"))

#double checking that there are not species specific orthogroups present in shared combined file:####
check<-annot.comb[,c(4,5)]
check$m<-paste0(check$Orthogroup,"_",check$Species)
t<-(unique(check$m))
t1<-as.data.frame(t)
t1<-colsplit(t1$t, "_", c("Orthogroup","Species"))
write.csv(t1, "shared.groupcheck.csv")#further checked visually/with count formulas in excel that there are idential orthogroups in each species. so these are in fact the shared orthogroups where there is at least one, maybe multiple transcripts belonging to each species.
countD<-subset(t1,Species=="Ds")#11782
countM<-subset(t1,Species=="Mb")#11782

#Making merged annotation files with GO info for functional analyses:####
Mb_GO<-read.delim(here("DDIG_data&analyses/DDIG_annotation&Functional_analyses/menidia_annotation_files/M_beryllina_B2G_final.annot"), header=F)

Ds_GO<-read.table(here("DDIG_data&analyses/DDIG_annotation&Functional_analyses/delta_smelt_annotation_files/delta_merged_final.annot.txt"),header=FALSE, sep = "\t",fill = TRUE,na.strings=c("","NA"))#need fill=TRUE here because unequal columns in input file

Mb_key<-read.csv(here("./DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthofinder_analyses/analysis_2020/Mb_ortho_key_by_transcript.csv"))
Mb_key<-Mb_key[,c(2,5)]
Ds_key<-read.csv(here("./DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthofinder_analyses/analysis_2020/Ds_ortho_key_by_transcript.csv"))
Ds_key<-Ds_key[,c(2,5)]

#note the GO ID's include evidence codes, so may need to filter out first, come back if determine it's needed
Mb_GO1<-Mb_GO[,1:2]
Ds_GO1<-Ds_GO[,1:2] 
Ds_GO1$Species<-"Ds"
Mb_GO1$Species<-"Mb"
Ds_GO1<-rename(Ds_GO1, GeneID_TranscriptID=V1,GO_terms=V2)
Ds_GO2<-merge(Ds_GO1,Ds_key)
Mb_GO1<-rename(Mb_GO1, GeneID_TranscriptID=V1,GO_terms=V2)
Mb_GO2<-merge(Mb_GO1,Mb_key)
ZZ<-rbind(Ds_GO2,Mb_GO2)
write.csv(ZZ, here("DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthofinder_analyses/analysis_2020/GO terms of Shared Orthogroups by Transcript_raw.csv"), row.names = F)

ZZ$GO_OG<-paste0(ZZ$GO_terms,"_",ZZ$Orthogroup)
ZZZ<-distinct(ZZ,GO_OG, .keep_all= TRUE)#get rid of duplicate GOID's within orthogroups
write.csv(ZZZ, here("DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthofinder_analyses/analysis_2020/GO terms of Shared Orthogroups by Transcript_nodups.csv"), row.names = F)

#Lisa start here: 
#1. look at TOPGO script to look at file input-how need GO terms? change input format accordingly
#2. Figure out how to pull in significant orthogroup lists for comparisons

#Next step (see following scripts): Perform DE on orthologues, make figures####
