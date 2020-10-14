#Differential expression analysis with limma 
#Resources: see Limma manual, and https://www.bioconductor.org/help/workflows/RNAseq123/

#References for code/ideas: (embed when morph this to cleaned RMarkdown doc)
#https://www.bioconductor.org/help/workflows/RNAseq123/#useful-graphical-representations-of-differential-expression-results

#I. Analyzing species + temp_time treatment####

# 1. Setup and Load the packages####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("edgeR")
#BiocManager::install("Glimma")
library(edgeR)
library(Glimma)
library(erer)
library(limma)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(reshape2)
library(here)
here()

#2. Read in data & metadata####
# 2a. Read in counts and sample information
counts<-read.csv(here("Differential_expression_analyses/Orthogroup_count_speciescomb_GHcp.csv"), row.names = 1)

#2b. Add in metadata:
#Can include as much companionate info as desired here to test different predictors, co-factors, confounding variables. But just (as always) remember that for inclusion in model, gets more challenging with more factors. Creating combined treatments here since we are focused on main effects between HC and temp challenges, mostly within each species.
metadata<-read.table(here("Differential_expression_analyses/Combined.DsMb.metadata_GHcp.txt"), header = T, stringsAsFactors = F)
metadata$Time_point<-factor(metadata$Time_point)
metadata$Temp_treat_level<-factor(metadata$Temp_treat_level)
metadata$Temp_time<-paste0(metadata$Temp_treat,"_",metadata$Time_point)

#3. Filtering & Initial Sanity Checks####
#3a. Filtering out lowly expressed and non-informative genes 
#Option 1-Filter out genes that are only present in 0 or 1 samples (softer filtering)
# filt <- which(apply(counts, 1, function(x)length(x[which(x > 0)])) <= 1)
# counts2 <- counts[-filt,]
# dim(counts) ## the number of features you started with
# dim(counts2) ## count counts of features you have left over after initial filter
#but when run this for Mb (through code below), see dip in low end of mean-variance relationship, indicating should do more filtering (see manual and quick explanation: https://stats.stackexchange.com/questions/160255/voom-mean-variance-trend-plot-how-to-interpret-the-plot)

#Option 2-stronger filtering-# In edgeR, it is recommended to remove features without at least 1 read per million in n of the samples, where n is the size of the smallest group of replicates'; I tested a suite of values to see when dip disappears; compared to the species specific DE analyses, don't have strong dip to start; tested across values 1-10, consistently gets smaller but probably just want to strike a balance of reducing but not getting rid of too much data. for now use 5 since that is the smallest group size (and what I used on the species specific analyses)

cpms <- cpm(counts)  ## counts per million
keep <- rowSums(cpms >1) >=5
counts2 <- counts[keep,]
dim(counts) ## the number of features you started with
dim(counts2) ## count counts of features you have left over after initial filter

#3b. Check filtering effects:
lcpm <- cpm(counts, log=TRUE)
nsamples <- ncol(counts)
col <- brewer.pal(nsamples, "Paired") #if have small # of samples is meaningful, but larger #s run out of unique colors in palette...helpful if see different pattern from one or a handful of samples and want to go back and figure out which ones are causing issues; could also come back and add loop to color by treatment to visualize different distributions by treatment here
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", samplenames, text.col=col, bty="n") #see note above about color use
lcpm <- cpm(counts2, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

#4. Normalize####
#From limma manual: 'It is usual to apply scale normalization to RNA-seq read counts, and the TMM normalization method in particular has been found to perform well in comparative studies. To apply TMM normalization, it is convenient to create a DGEList object using the edgeR package:'

ed <- DGEList(counts2)# Create a DGEList object (edgeR's container for RNA-seq count data)
ed <- calcNormFactors(ed)# Calculate normalization factors
ed$samples #print list of the calculated library sizes and normalization factors for each of the samples

#5. Set up design matrix####
metadata$Condition<-factor(paste(metadata$Temp_time,metadata$Species,sep="."))#combining temp time because want to treat as one condition (only looking at main effects, no interactions for our questions)
summary(metadata$Condition)
mod.mat<-model.matrix(~0+Condition, data=metadata) 
colnames(mod.mat)<-levels(metadata$Condition)
colnames(mod.mat)

#Set up contrasts:
my.contrasts <- makeContrasts(
  CTMax2.0.HC.Ds = CTMax2_0.Ds-HC_0.Ds,
  CTMax4.0.HC.Ds = CTMax4_0.Ds-HC_0.Ds,
  CTMax6.0.HC.Ds = CTMax6_0.Ds-HC_0.Ds,
  CTMax2.60.HC.Ds = CTMax2_60.Ds-HC_60.Ds,
  CTMax4.60.HC.Ds = CTMax4_60.Ds-HC_60.Ds,
  CTMax6.60.HC.Ds = CTMax6_60.Ds-HC_60.Ds,
  CTMax2.0.HC.Mb = CTMax2_0.Mb-HC_0.Mb,
  CTMax4.0.HC.Mb = CTMax4_0.Mb-HC_0.Mb,
  CTMax8.0.HC.Mb = CTMax8_0.Mb-HC_0.Mb,
  CTMax12.0.HC.Mb = CTMax12_0.Mb-HC_0.Mb,
  CTMax2.60.HC.Mb = CTMax2_60.Mb-HC_60.Mb,
  CTMax4.60.HC.Mb = CTMax4_60.Mb-HC_60.Mb,
  CTMax8.60.HC.Mb = CTMax8_60.Mb-HC_60.Mb,
  CTMax12.60.HC.Mb = CTMax12_60.Mb-HC_60.Mb,
  HC.0.Ds.Mb = HC_0.Ds-HC_0.Mb,
  HC.60.Ds.Mb = HC_60.Ds-HC_60.Mb,
  CTMax2.0.Ds.Mb = CTMax2_0.Mb-CTMax2_0.Ds,
  CTMax4.0.Ds.Mb = CTMax4_0.Mb-CTMax4_0.Ds,
  CTMax6_8.0.Ds.Mb = CTMax8_0.Mb-CTMax6_0.Ds,
  CTMax6_12.0.Ds.Mb = CTMax12_0.Mb-CTMax6_0.Ds,
  CTMax2.60.Ds.Mb = CTMax2_60.Mb-CTMax2_60.Ds,
  CTMax4.60.Ds.Mb = CTMax4_60.Mb-CTMax4_60.Ds,
  CTMax6_8.60.Ds.Mb = CTMax8_60.Mb-CTMax6_60.Ds,
  CTMax6_12.60.Ds.Mb = CTMax12_60.Mb-CTMax6_60.Ds,
  levels=colnames(mod.mat))

my.contrasts

#6. Apply the voom transformation, which “converts the counts to log-counts per million with associated precision weights”####
par(mfrow=c(2,1))
v <- voom(ed,mod.mat,plot=TRUE)
vfit <- lmFit(v,mod.mat)
vfit <- contrasts.fit(vfit, contrasts=my.contrasts)
efit <- eBayes(vfit)

par(mfrow=c(2,1))
v <- voom(ed,mod.mat,plot=TRUE)
plotSA(efit, main="Final model: Mean−variance trend")
par(mfrow=c(1,1))#just clearing format

#7. Assess differentially expressed genes:####
summary(decideTests(efit))#this gives a nice quick summary table of how many genes are significant for each contrast (default sig set at adjusted 5% p-value threshold)
ds<-decideTests(efit)#store sig genes in an object for later use
#ds1<-decideTests(efit, adjust.method = "BH",p.value = 0.05)#checking the defaults of this function
#summary(ds1)#-yep same summary results
#ds2<-decideTests(efit, adjust.method = "fdr",p.value = 0.05)#manual says FDR should be the same as BH method
#summary(ds2)#yep same results

#save for later inspection:
sig<-as.data.frame(summary(decideTests(efit)))
#write.table(sig,file=here("DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthogroup_DE_analyses/Temptime_splitanalysisorthogroups_summary_sig_genes.txt"), sep = "\t", quote = F, row.names = F)

#alternatively may want to set significance to be above a certain log fold change cutoff
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

#Examining individual DE genes from top to bottom
voom_results_CTMax2.0.HC.Ds <- topTable(efit, coef=1, n=Inf)
head(voom_results_CTMax2.0.HC.Ds)
colnames(efit)#to determine which contrasts are which coef to call them (in order below)
row<-rownames(efit)#to make orthogroup rowname for all below

#set for graphs:
voom_results_CTMax2.0.HC.Ds <- topTable(efit, coef=1, n=Inf)
voom_results_CTMax4.0.HC.Ds <- topTable(efit, coef=2, n=Inf)
voom_results_CTMax6.0.HC.Ds <- topTable(efit, coef=3, n=Inf) 
voom_results_CTMax2.60.HC.Ds <- topTable(efit, coef=4, n=Inf) 
voom_results_CTMax4.60.HC.Ds<- topTable(efit, coef=5, n=Inf)
voom_results_CTMax6.60.HC.Ds <- topTable(efit, coef=6, n=Inf)
voom_results_CTMax2.0.HC.Mb <- topTable(efit, coef=7, n=Inf)
voom_results_CTMax4.0.HC.Mb <- topTable(efit, coef=8, n=Inf)
voom_results_CTMax8.0.HC.Mb <- topTable(efit, coef=9, n=Inf)
voom_results_CTMax12.0.HC.Mb <- topTable(efit, coef=10, n=Inf)
voom_results_CTMax2.60.HC.Mb <- topTable(efit, coef=11, n=Inf)
voom_results_CTMax4.60.HC.Mb <- topTable(efit, coef=12, n=Inf)
voom_results_CTMax8.60.HC.Mb <- topTable(efit, coef=13, n=Inf)
voom_results_CTMax12.60.HC.Mb <- topTable(efit, coef=14, n=Inf)
voom_results_HC.0.Ds.Mb <- topTable(efit, coef=15, n=Inf)
voom_results_HC.60.Ds.Mb <- topTable(efit, coef=16, n=Inf)
voom_results_CTMax2.0.Ds.Mb <- topTable(efit, coef=17, n=Inf)
voom_results_CTMax4.0.Ds.Mb <- topTable(efit, coef=18, n=Inf)
voom_results_CTMax6_8.0.Ds.Mb <- topTable(efit, coef=19, n=Inf)
voom_results_CTMax6_12.0.Ds.Mb <- topTable(efit, coef=20, n=Inf)
voom_results_CTMax2.60.Ds.Mb <- topTable(efit, coef=21, n=Inf)
voom_results_CTMax4.60.Ds.Mb <- topTable(efit, coef=22, n=Inf)
voom_results_CTMax6_8.60.Ds.Mb <- topTable(efit, coef=23, n=Inf)
voom_results_CTMax6_12.60.Ds.Mb <- topTable(efit, coef=24, n=Inf)

#write.table(voom_results_CTMax2.0.HC.Ds, paste(here("DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthogroup_DE_analyses/Temptime_splitanalysis/voom_results_CTMax2.0.HC.Ds.txt"), sep = ""), col.names = NA, row.names = TRUE, sep = "\t", quote= FALSE)#repeat for others if want seperate for each contrast

#to write results to file all together:
#write.fit(efit, ds, file=here("DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthogroup_DE_analyses/Temptime_splitanalysis/orthogroup.voom.test.results_rm2rds.txt"),adjust="BH")
#Help for this function says: This function writes a tab-delimited text file containing for each gene (1) the average log-intensity, (2) the log-ratios, (3) moderated t-statistics, (4) t-statistic P-values, (5) F-statistic if available, (6) F-statistic P-values if available, (7) classification if available and (8) gene names and annotation.'

#what genes are common across the contrasts?
head(ds)#to checks which columns are which contrasts
de.common <- which(ds[,1]!=0 & ds[,2]!=0)#which genes aren't equal to 0 (=not sig) that are common to column 1 & 2
length(de.common)#indicating 12 common sig genes between CTMax2.0.HC.Ds and CTMax4.0.HC.Ds lists

#8. Visualize####
#8a. Venn diagram showing the number of genes DE in the comparison of specified contrasts, and the number of genes that are DE in both comparisons (center). The number of genes that are not DE in either comparison are marked in the bottom-right.

#within species:
pdf(here("DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthogroup_DE_analyses/Temptime_splitanalysis/venndiagrams_withinspecies.pdf"), width=12, height=6)
par(mfrow=c(2,2))
vennDiagram(ds[,1:3], circle.col=c("salmon","yellow","turquoise"))
vennDiagram(ds[,4:6], circle.col=c("salmon","yellow","turquoise"))
vennDiagram(ds[,7:10], circle.col=c("salmon","yellow","turquoise","green" ))
vennDiagram(ds[,11:14], circle.col=c("salmon","yellow","turquoise","green" ))
dev.off()

pdf(here("DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthogroup_DE_analyses/Temptime_splitanalysis/venndiagrams_acrosstimepoints.pdf"), width=12, height=6)
par(mfrow=c(2,4))
vennDiagram(ds[,c(1,4)], circle.col=c("salmon","yellow"))
vennDiagram(ds[,c(2,5)], circle.col=c("salmon","yellow"))
vennDiagram(ds[,c(3,6)], circle.col=c("salmon","yellow"))
vennDiagram(ds[,c(3,6)], circle.col=c("salmon","yellow"))

vennDiagram(ds[,c(7,11)], circle.col=c("salmon","yellow"))
vennDiagram(ds[,c(8,12)], circle.col=c("salmon","yellow"))
vennDiagram(ds[,c(9,13)], circle.col=c("salmon","yellow"))
vennDiagram(ds[,c(10,14)], circle.col=c("salmon","yellow"))
dev.off()

pdf(here("DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthogroup_DE_analyses/Temptime_splitanalysis/venndiagrams_acrossspecies_sametreatment.pdf"), width=12, height=6)
par(mfrow=c(2,3))
vennDiagram(ds[,c(1,7)], circle.col=c("salmon","yellow"))
vennDiagram(ds[,c(2,8)], circle.col=c("salmon","yellow"))
vennDiagram(ds[,c(3,9)], circle.col=c("salmon","yellow"))

vennDiagram(ds[,c(4,11)], circle.col=c("salmon","yellow"))
vennDiagram(ds[,c(5,12)], circle.col=c("salmon","yellow"))
vennDiagram(ds[,c(6,13)], circle.col=c("salmon","yellow"))
dev.off()

#8b. MA Plots-Useful graphical representations of differential expression results:####
par(mfrow=c(1,1))
plotMD(efit, column=1,  #column of object to be plotted
       status=ds[,1],#sig status (-1,0,1) for sig down, not sig, or sig up 
       main=colnames(efit)[1])+#contrast for title of graph
  abline(h=c(-1, 1), col="blue", lty="dashed")
#repeat for any others want to explore here

#Figure:
#within Ds & HC between species:
pdf(here("DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthogroup_DE_analyses/Temptime_splitanalysis/temptreat_matrix_withinDs&HCDsMb_MAplotspanel.pdf"), width=12, height=6)
par(mfrow=c(2,4))
plotMD(efit, column=3, status=ds[,3],col=c("red","blue1"),ylim=c(-8,8),legend=FALSE,  main=colnames(efit)[3])+ abline(h=c(-1, 1), lty="dashed",col="gray40")
plotMD(efit, column=2, status=ds[,2],col=c("red","blue1"),ylim=c(-8,8),legend=FALSE,  main=colnames(efit)[2])+ abline(h=c(-1, 1), lty="dashed",col="gray40")
plotMD(efit, column=1, status=ds[,1],col=c("red","blue1"),ylim=c(-8,8),legend=FALSE, main=colnames(efit)[1])+ abline(h=c(-1, 1),  lty="dashed",col="gray40")
plotMD(efit, column=15, status=ds[,15],col=c("red","blue1"),ylim=c(-20,20),legend=FALSE,  main=colnames(efit)[15])+ abline(h=c(-1, 1), lty="dashed",col="gray40")#note changed legend here, to pull out these figures seperately in post-processing
plotMD(efit, column=6, status=ds[,6],col=c("red","blue1"),ylim=c(-8,8),legend=FALSE,  main=colnames(efit)[6])+ abline(h=c(-1, 1), lty="dashed",col="gray40")
plotMD(efit, column=5, status=ds[,5],col=c("red","blue1"),ylim=c(-8,8),legend=FALSE,  main=colnames(efit)[5])+ abline(h=c(-1, 1),  lty="dashed",col="gray40")
plotMD(efit, column=4, status=ds[,4],col=c("red","blue1"),ylim=c(-8,8),legend=FALSE,  main=colnames(efit)[4])+ abline(h=c(-1, 1),lty="dashed",col="gray40")
plotMD(efit, column=16, status=ds[,16],col=c("red","blue1"),ylim=c(-20,20),legend=FALSE,  main=colnames(efit)[16])+ abline(h=c(-1, 1), lty="dashed",col="gray40")#note changed legend here, to pull out these figures seperately in post-processing
dev.off()

#within Mb
pdf(here("DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthogroup_DE_analyses/Temptime_splitanalysis/temptreat_matrix_withinMb_MAplotspanel.pdf"), width=12, height=6)
par(mfrow=c(2,4))
plotMD(efit, column=10, status=ds[,10],col=c("red","blue1"),ylim=c(-8,8),legend=FALSE,  main=colnames(efit)[10])+ abline(h=c(-1, 1), lty="dashed",col="gray40")
plotMD(efit, column=9, status=ds[,9],col=c("red","blue1"),ylim=c(-8,8),legend=FALSE,  main=colnames(efit)[9])+ abline(h=c(-1, 1), lty="dashed",col="gray40")
plotMD(efit, column=8, status=ds[,8],col=c("red","blue1"),ylim=c(-8,8),legend=FALSE,  main=colnames(efit)[8])+ abline(h=c(-1, 1),  lty="dashed",col="gray40")
plotMD(efit, column=7, status=ds[,7],col=c("red","blue1"),ylim=c(-8,8),legend=FALSE,  main=colnames(efit)[7])+ abline(h=c(-1, 1),lty="dashed",col="gray40")
plotMD(efit, column=14, status=ds[,14],col=c("red","blue1"),ylim=c(-8,8),legend=FALSE,  main=colnames(efit)[14])+ abline(h=c(-1, 1),lty="dashed",col="gray40")
plotMD(efit, column=13, status=ds[,13],col=c("red","blue1"),ylim=c(-8,8),legend=FALSE,  main=colnames(efit)[13])+ abline(h=c(-1, 1),lty="dashed",col="gray40")
plotMD(efit, column=12, status=ds[,12],col=c("red","blue1"),ylim=c(-8,8),legend=FALSE,  main=colnames(efit)[12])+ abline(h=c(-1, 1),lty="dashed",col="gray40")
plotMD(efit, column=11, status=ds[,11],col=c("red","blue1"),ylim=c(-8,8),legend=FALSE,  main=colnames(efit)[11])+ abline(h=c(-1, 1),lty="dashed",col="gray40")
dev.off()

#8c. Heatmaps: see other script####

#9.Pulling in annotation:####
#note: the here function was throwing all sorts of weird errors below so I just scrapped it and checked the WD to make the appropriate relative file paths
library(plyr)
library(tidyr)
getwd()
annot.comb<-read.csv("DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthofinder_analyses/analysis_2020/Annotations of Shared Orthogroups by Transcript.csv")
annot.comb1<-annot.comb[,c(2,4)]
annot.unique<-(unique(annot.comb1))#removing redudancies
annot.unique1<-annot.unique[complete.cases(annot.unique[ , 1]),]
m1<-aggregate(Sequence.desc. ~ Orthogroup, data = annot.unique1, FUN = paste, collapse = ",")
results<-read.table("DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthofinder_analyses/analysis_2020/orthogroup.voom.test.results_rm2rds.txt",header=TRUE)
Orthogroup<-rownames(efit)
c<-cbind(Orthogroup, results)
annot_compiled_contrasts<-merge(m1,c,by="Orthogroup",all.y=TRUE)#merging with annotation but keeping data even for those without annotation
write.table(annot_compiled_contrasts, file = "DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthogroup_DE_analyses/Temptime_splitanalysis/Orthogroup_annot_compiled_contrastslimma_rm2rds.txt", sep = "\t", quote = F, row.names = F)
