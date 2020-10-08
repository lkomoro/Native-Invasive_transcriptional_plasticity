##Enriched GO Comparisons/Analyses####
#Reading in files that are output from TopGO analyses
#rm(list=ls())
library(here)
here()
setwd(here("./DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthogroup_TopGO_analyses"))

#KS Test Results:
list.files(pattern="KSGO.txt")

CTMax2.0.HC.Ds<-read.delim("CTMax2.0.HC.Ds_KSGO.txt", header=T)
Key<-CTMax2.0.HC.Ds[,c(1:3,5)]

colnames(CTMax2.0.HC.Ds)[4] <- "CTMax2.0.HC.Dspval"
CTMax2.0.HC.Ds<-CTMax2.0.HC.Ds[,c(1,4)]

CTMax2.60.HC.Ds<-read.delim("CTMax2.60.HC.Ds_KSGO.txt", header=T)
colnames(CTMax2.60.HC.Ds)[4] <- "CTMax2.60.HC.Dspval"
CTMax2.60.HC.Ds<-CTMax2.60.HC.Ds[,c(1,4)]

CTMax4.0.HC.Ds<-read.delim("CTMax4.0.HC.Ds_KSGO.txt", header=T)
colnames(CTMax4.0.HC.Ds)[4] <- "CTMax4.0.HC.Dspval"
CTMax4.0.HC.Ds<-CTMax4.0.HC.Ds[,c(1,4)]

CTMax4.60.HC.Ds<-read.delim("CTMax4.60.HC.Ds_KSGO.txt", header=T)
colnames(CTMax4.60.HC.Ds)[4] <- "CTMax4.60.HC.Dspval"
CTMax4.60.HC.Ds<-CTMax4.60.HC.Ds[,c(1,4)]

CTMax6.0.HC.Ds<-read.delim("CTMax6.0.HC.Ds_KSGO.txt", header=T)
colnames(CTMax6.0.HC.Ds)[4] <- "CTMax6.0.HC.Dspval"
CTMax6.0.HC.Ds<-CTMax6.0.HC.Ds[,c(1,4)]

CTMax6.60.HC.Ds<-read.delim("CTMax6.60.HC.Ds_KSGO.txt", header=T)
colnames(CTMax6.60.HC.Ds)[4] <- "CTMax6.60.HC.Dspval"
CTMax6.60.HC.Ds<-CTMax6.60.HC.Ds[,c(1,4)]

CTMax2.0.HC.Mb<-read.delim("CTMax2.0.HC.Mb_KSGO.txt", header=T)
colnames(CTMax2.0.HC.Mb)[4] <- "CTMax2.0.HC.Mbpval"
CTMax2.0.HC.Mb<-CTMax2.0.HC.Mb[,c(1,4)]

CTMax2.60.HC.Mb<-read.delim("CTMax2.60.HC.Mb_KSGO.txt", header=T)
colnames(CTMax2.60.HC.Mb)[4] <- "CTMax2.60.HC.Mbpval"
CTMax2.60.HC.Mb<-CTMax2.60.HC.Mb[,c(1,4)]

CTMax4.0.HC.Mb<-read.delim("CTMax4.0.HC.Mb_KSGO.txt", header=T)
colnames(CTMax4.0.HC.Mb)[4] <- "CTMax4.0.HC.Mbpval"
CTMax4.0.HC.Mb<-CTMax4.0.HC.Mb[,c(1,4)]

CTMax4.60.HC.Mb<-read.delim("CTMax4.60.HC.Mb_KSGO.txt", header=T)
colnames(CTMax4.60.HC.Mb)[4] <- "CTMax4.60.HC.Mbpval"
CTMax4.60.HC.Mb<-CTMax4.60.HC.Mb[,c(1,4)]

CTMax8.0.HC.Mb<-read.delim("CTMax8.0.HC.Mb_KSGO.txt", header=T)
colnames(CTMax8.0.HC.Mb)[4] <- "CTMax8.0.HC.Mbpval"
CTMax8.0.HC.Mb<-CTMax8.0.HC.Mb[,c(1,4)]

CTMax8.60.HC.Mb<-read.delim("CTMax8.60.HC.Mb_KSGO.txt", header=T)
colnames(CTMax8.60.HC.Mb)[4] <- "CTMax8.60.HC.Mbpval"
CTMax8.60.HC.Mb<-CTMax8.60.HC.Mb[,c(1,4)]

CTMax12.0.HC.Mb<-read.delim("CTMax12.0.HC.Mb_KSGO.txt", header=T)
colnames(CTMax12.0.HC.Mb)[4] <- "CTMax12.0.HC.Mbpval"
CTMax12.0.HC.Mb<-CTMax12.0.HC.Mb[,c(1,4)]

CTMax12.60.HC.Mb<-read.delim("CTMax12.60.HC.Mb_KSGO.txt", header=T)
colnames(CTMax12.60.HC.Mb)[4] <- "CTMax12.60.HC.Mbpval"
CTMax12.60.HC.Mb<-CTMax12.60.HC.Mb[,c(1,4)]

ls(pattern="CT")#just to get object list to check/paste names to combine
comb<-merge(Key,CTMax2.0.HC.Ds)
comb<-merge(comb,CTMax4.0.HC.Ds)
comb<-merge(comb,CTMax6.0.HC.Ds)
comb<-merge(comb,CTMax2.60.HC.Ds)
comb<-merge(comb,CTMax4.60.HC.Ds)
comb<-merge(comb,CTMax6.60.HC.Ds)
comb<-merge(comb,CTMax2.0.HC.Mb)
comb<-merge(comb,CTMax4.0.HC.Mb)
comb<-merge(comb,CTMax8.0.HC.Mb)
comb<-merge(comb,CTMax12.0.HC.Mb)
comb<-merge(comb,CTMax2.60.HC.Mb)
comb<-merge(comb,CTMax4.60.HC.Mb)
comb<-merge(comb,CTMax8.60.HC.Mb)
comb<-merge(comb,CTMax12.60.HC.Mb)

#subset & venn diagram
pval<-0.05 
#pval<-0.01 #for stricter threshold
CTMax2.0.HC.Dssig<-subset(CTMax2.0.HC.Ds, CTMax2.0.HC.Dspval<=pval)
CTMax2.0.HC.Dssig<-CTMax2.0.HC.Dssig[,1]

CTMax4.0.HC.Dssig<-subset(CTMax4.0.HC.Ds, CTMax4.0.HC.Dspval<=pval)
CTMax4.0.HC.Dssig<-CTMax4.0.HC.Dssig[,1]

CTMax6.0.HC.Dssig<-subset(CTMax6.0.HC.Ds, CTMax6.0.HC.Dspval<=pval)
CTMax6.0.HC.Dssig<-CTMax6.0.HC.Dssig[,1]

CTMax2.60.HC.Dssig<-subset(CTMax2.60.HC.Ds, CTMax2.60.HC.Dspval<=pval)
CTMax2.60.HC.Dssig<-CTMax2.60.HC.Dssig[,1]

CTMax4.60.HC.Dssig<-subset(CTMax4.60.HC.Ds, CTMax4.60.HC.Dspval<=pval)
CTMax4.60.HC.Dssig<-CTMax4.60.HC.Dssig[,1]

CTMax6.60.HC.Dssig<-subset(CTMax6.60.HC.Ds, CTMax6.60.HC.Dspval<=pval)
CTMax6.60.HC.Dssig<-CTMax6.60.HC.Dssig[,1]

CTMax2.0.HC.Mbsig<-subset(CTMax2.0.HC.Mb, CTMax2.0.HC.Mbpval<=pval)
CTMax2.0.HC.Mbsig<-CTMax2.0.HC.Mbsig[,1]

CTMax2.60.HC.Mbsig<-subset(CTMax2.60.HC.Mb, CTMax2.60.HC.Mbpval<=pval)
CTMax2.60.HC.Mbsig<-CTMax2.60.HC.Mbsig[,1]

CTMax4.0.HC.Mbsig<-subset(CTMax4.0.HC.Mb, CTMax4.0.HC.Mbpval<=pval)
CTMax4.0.HC.Mbsig<-CTMax4.0.HC.Mbsig[,1]

CTMax4.60.HC.Mbsig<-subset(CTMax4.60.HC.Mb, CTMax4.60.HC.Mbpval<=pval)
CTMax4.60.HC.Mbsig<-CTMax4.60.HC.Mbsig[,1]

CTMax8.0.HC.Mbsig<-subset(CTMax8.0.HC.Mb, CTMax8.0.HC.Mbpval<=pval)
CTMax8.0.HC.Mbsig<-CTMax8.0.HC.Mbsig[,1]

CTMax8.60.HC.Mbsig<-subset(CTMax8.60.HC.Mb, CTMax8.60.HC.Mbpval<=pval)
CTMax8.60.HC.Mbsig<-CTMax8.60.HC.Mbsig[,1]

CTMax12.0.HC.Mbsig<-subset(CTMax12.0.HC.Mb, CTMax12.0.HC.Mbpval<=pval)
CTMax12.0.HC.Mbsig<-CTMax12.0.HC.Mbsig[,1]

CTMax12.60.HC.Mbsig<-subset(CTMax12.60.HC.Mb, CTMax12.60.HC.Mbpval<=pval)
CTMax12.60.HC.Mbsig<-CTMax12.60.HC.Mbsig[,1]

#use/put in apply function if want #s for table:
#length(intersect(CTMax2.HC.60sig, CTMax2.HC.0sig))
#length(setdiff(CTMax2.HC.60sig, CTMax2.HC.0sig))
#length(setdiff(CTMax2.HC.0sig, CTMax2.HC.60sig))


require("gplots")
ls(pattern="sig")#just to get object list to check/paste names to combine
venn(list(CTMax2.0.HC.Ds = CTMax2.0.HC.Dssig, CTMax4.0.HC.Ds = CTMax4.0.HC.Dssig, CTMax6.0.HC.Ds = CTMax6.0.HC.Dssig))
venn(list(CTMax2.60.HC.Ds = CTMax2.60.HC.Dssig, CTMax4.60.HC.Ds = CTMax4.60.HC.Dssig, CTMax6.60.HC.Ds = CTMax6.60.HC.Dssig))

venn(list(CTMax2.0.HC.Mb = CTMax2.0.HC.Mbsig, CTMax4.0.HC.Mb = CTMax4.0.HC.Mbsig, CTMax8.0.HC.Mb = CTMax8.0.HC.Mbsig,CTMax12.0.HC.Mb = CTMax12.0.HC.Mbsig))
venn(list(CTMax2.60.HC.Mb = CTMax2.60.HC.Mbsig, CTMax4.60.HC.Mb = CTMax4.60.HC.Mbsig, CTMax8.60.HC.Mb = CTMax8.60.HC.Mbsig,CTMax12.60.HC.Mb = CTMax12.60.HC.Mbsig))

#Combine by temp treatment
CTMax2.0 <- unique(c(CTMax2.0.HC.Dssig,CTMax2.0.HC.Mbsig))
length(CTMax2.0)
venn(list(CTMax2.0.HC.Ds = CTMax2.0.HC.Dssig,CTMax2.0.HC.Mb = CTMax2.0.HC.Mbsig))

CTMax4.0 <- unique(c(CTMax4.0.HC.Dssig,CTMax4.0.HC.Mbsig))
length(CTMax4.0)
venn(list(CTMax4.0.HC.Ds = CTMax4.0.HC.Dssig,CTMax4.0.HC.Mb = CTMax4.0.HC.Mbsig))

CTMax6.0 <- unique(c(CTMax6.0.HC.Dssig,CTMax8.0.HC.Mbsig))
length(CTMax6.0)
venn(list(CTMax6.0.HC.Ds = CTMax6.0.HC.Dssig,CTMax8.0.HC.Mb = CTMax8.0.HC.Mbsig))

#add in others to compare as desired

###############################
t1<-comb[,c(1,5:18)]
t<-t1[apply(t1[, -1], MARGIN = 1, function(x) any(x <= pval)), ]#save only with at least one treatment at or lower than set threshold
comb1<-comb[,1:4]
comb2<-merge(comb1,t)
write.csv(comb2,"Combined_orthogroup_sigGOenrichKSp0.05.csv")

################################################
#Note that some treatments had low numbers of significant genes based on 0.05 adj FDR, and perhaps it doesn't make sense to include results (especially in heatmap)?
#CTMax4.0.Ds: 14 sig genes
#CTMax6.0.Ds: 7 sig genes
#CTMax6.60.Ds: 0 sig genes
#CTMax8.0.HC.Mb: 9 sig genes
#CTMax12.0.HC.Mb: 11 sig genes
#(see saved table with added column: orthogroups_summary_sig_genes_LK.xlsx)
###############################
