############Distributions & Heatmaps of DE Genes###########
#Set Conditions####
#rm(list=ls())
library(gplots)
library(RColorBrewer)
library(plyr)#need this for the rename function-note that if dplyr is loaded (and/or maybe tidyr?) it will override and cause an error bc different syntax
library(ggplot2)
library(tidyr)

#data readin & setup:####
voomresults<-read.delim("Differential_expression_analyses/Orthogroup_annot_compiled_contrasts_GHcp.txt")
grep("adj",names(voomresults),value=T)
adj<-grep("adj",names(voomresults))#extract #s of which columns I want
voomresults.adjp<-voomresults[,c(1,76:89)]#keeping just within species contrasts
grep("adj",names(voomresults.adjp),value=T)#double check correct columns
voomresults.key<-voomresults[,1:2]
siggenes<-voomresults.adjp[apply(voomresults.adjp[, -1], MARGIN = 1, function(x) any(x <= 0.05)), ]#keep only genes where at least one treatment has FDR<=0.05
siggeneskey<-merge(voomresults.key,siggenes)
siggeneskey<-siggeneskey[,1:2]

#Ds:
CTMax2.0.HC.Ds<-read.delim("Differential_expression_analyses/voom_results_CTMax2.0.HC.Ds.txt")
CTMax2.0.HC.Ds<-CTMax2.0.HC.Ds[,1:2]
CTMax2.0.HC.Ds<-rename(CTMax2.0.HC.Ds, c("X"="Orthogroup", "logFC"="CTMax2.0.HC.DslogFC"))

CTMax4.0.HC.Ds<-read.delim("Differential_expression_analyses/voom_results_CTMax4.0.HC.Ds.txt")
CTMax4.0.HC.Ds<-CTMax4.0.HC.Ds[,1:2]
CTMax4.0.HC.Ds<-rename(CTMax4.0.HC.Ds, c("X"="Orthogroup", "logFC"="CTMax4.0.HC.DslogFC"))

CTMax6.0.HC.Ds<-read.delim("Differential_expression_analyses/voom_results_CTMax6.0.HC.Ds.txt")
CTMax6.0.HC.Ds<-CTMax6.0.HC.Ds[,1:2]
CTMax6.0.HC.Ds<-rename(CTMax6.0.HC.Ds, c("X"="Orthogroup", "logFC"="CTMax6.0.HC.DslogFC"))

CTMax2.60.HC.Ds<-read.delim("Differential_expression_analyses/voom_results_CTMax2.60.HC.Ds.txt")
CTMax2.60.HC.Ds<-CTMax2.60.HC.Ds[,1:2]
CTMax2.60.HC.Ds<-rename(CTMax2.60.HC.Ds, c("X"="Orthogroup", "logFC"="CTMax2.60.HC.DslogFC"))

CTMax4.60.HC.Ds<-read.delim("Differential_expression_analyses/voom_results_CTMax4.60.HC.Ds.txt")
CTMax4.60.HC.Ds<-CTMax4.60.HC.Ds[,1:2]
CTMax4.60.HC.Ds<-rename(CTMax4.60.HC.Ds, c("X"="Orthogroup", "logFC"="CTMax4.60.HC.DslogFC"))

CTMax6.60.HC.Ds<-read.delim("Differential_expression_analyses/voom_results_CTMax6.60.HC.Ds.txt")
CTMax6.60.HC.Ds<-CTMax6.60.HC.Ds[,1:2]
CTMax6.60.HC.Ds<-rename(CTMax6.60.HC.Ds, c("X"="Orthogroup", "logFC"="CTMax6.60.HC.DslogFC"))

CTMax6.60.HC.Ds<-read.delim("Differential_expression_analyses/voom_results_CTMax6.60.HC.Ds.txt")
CTMax6.60.HC.Ds<-CTMax6.60.HC.Ds[,1:2]
CTMax6.60.HC.Ds<-rename(CTMax6.60.HC.Ds, c("X"="Orthogroup", "logFC"="CTMax6.60.HC.DslogFC"))

#Mb:
CTMax2.0.HC.Mb<-read.delim("Differential_expression_analyses/voom_results_CTMax2.0.HC.Mb.txt")
CTMax2.0.HC.Mb<-CTMax2.0.HC.Mb[,1:2]
CTMax2.0.HC.Mb<-rename(CTMax2.0.HC.Mb, c("X"="Orthogroup", "logFC"="CTMax2.0.HC.MblogFC"))

CTMax4.0.HC.Mb<-read.delim("Differential_expression_analyses/voom_results_CTMax4.0.HC.Mb.txt")
CTMax4.0.HC.Mb<-CTMax4.0.HC.Mb[,1:2]
CTMax4.0.HC.Mb<-rename(CTMax4.0.HC.Mb, c("X"="Orthogroup", "logFC"="CTMax4.0.HC.MblogFC"))

CTMax8.0.HC.Mb<-read.delim("Differential_expression_analyses/voom_results_CTMax8.0.HC.Mb.txt")
CTMax8.0.HC.Mb<-CTMax8.0.HC.Mb[,1:2]
CTMax8.0.HC.Mb<-rename(CTMax8.0.HC.Mb, c("X"="Orthogroup", "logFC"="CTMax8.0.HC.MblogFC"))

CTMax12.0.HC.Mb<-read.delim("Differential_expression_analyses/voom_results_CTMax12.0.HC.Mb.txt")
CTMax12.0.HC.Mb<-CTMax12.0.HC.Mb[,1:2]
CTMax12.0.HC.Mb<-rename(CTMax12.0.HC.Mb, c("X"="Orthogroup", "logFC"="CTMax12.0.HC.MblogFC"))

CTMax2.60.HC.Mb<-read.delim("Differential_expression_analyses/voom_results_CTMax2.60.HC.Mb.txt")
CTMax2.60.HC.Mb<-CTMax2.60.HC.Mb[,1:2]
CTMax2.60.HC.Mb<-rename(CTMax2.60.HC.Mb, c("X"="Orthogroup", "logFC"="CTMax2.60.HC.MblogFC"))

CTMax4.60.HC.Mb<-read.delim("Differential_expression_analyses/voom_results_CTMax4.60.HC.Mb.txt")
CTMax4.60.HC.Mb<-CTMax4.60.HC.Mb[,1:2]
CTMax4.60.HC.Mb<-rename(CTMax4.60.HC.Mb, c("X"="Orthogroup", "logFC"="CTMax4.60.HC.MblogFC"))

CTMax8.60.HC.Mb<-read.delim("Differential_expression_analyses/voom_results_CTMax8.60.HC.Mb.txt")
CTMax8.60.HC.Mb<-CTMax8.60.HC.Mb[,1:2]
CTMax8.60.HC.Mb<-rename(CTMax8.60.HC.Mb, c("X"="Orthogroup", "logFC"="CTMax8.60.HC.MblogFC"))

CTMax12.60.HC.Mb<-read.delim("Differential_expression_analyses/voom_results_CTMax12.60.HC.Mb.txt")
CTMax12.60.HC.Mb<-CTMax12.60.HC.Mb[,1:2]
CTMax12.60.HC.Mb<-rename(CTMax12.60.HC.Mb, c("X"="Orthogroup", "logFC"="CTMax12.60.HC.MblogFC"))

#################################################################
#Figure 4, Distributions: (post processing in Adobe)####
#1. Delta smelt
Ds.voomresults.adjp<-voomresults.adjp[,1:7]
Ds.siggenes<-Ds.voomresults.adjp[apply(Ds.voomresults.adjp[, -1], MARGIN = 1, function(x) any(x <= 0.05)), ]
Ds.siggeneskey<-merge(voomresults.key,Ds.siggenes)
Ds.siggeneskey<-Ds.siggeneskey[,1:2]
c<-merge(CTMax2.0.HC.Ds,CTMax4.0.HC.Ds)
c<-merge(c,CTMax6.0.HC.Ds)
c<-merge(c,CTMax2.60.HC.Ds)
c<-merge(c,CTMax4.60.HC.Ds)
c<-merge(c,CTMax6.60.HC.Ds)
Ds.siggeneFC<-merge(Ds.siggeneskey,c)
#write.csv(Ds.siggeneskey,"orthogroup.adjp0.05.Ds.sig.genes.key.csv",row.names = F)
Ds.siggeneFC.L<-gather(Ds.siggeneFC, key = "Treatment", value = "logFC", CTMax6.0.HC.DslogFC,CTMax4.0.HC.DslogFC,CTMax2.0.HC.DslogFC,CTMax6.60.HC.DslogFC,CTMax4.60.HC.DslogFC,CTMax2.60.HC.DslogFC)

Ds.siggeneFC.L$Treatment<-factor(Ds.siggeneFC.L$Treatment)
levels(Ds.siggeneFC.L$Treatment)
Ds.siggeneFC.L$Treatment<-ordered(Ds.siggeneFC.L$Treatment,c("CTMax6.0.HC.DslogFC"  , "CTMax6.60.HC.DslogFC",  "CTMax4.0.HC.DslogFC",  "CTMax4.60.HC.DslogFC","CTMax2.0.HC.DslogFC" ,"CTMax2.60.HC.DslogFC"))

#pdf("orthogroup.dist.logFC.Ds.0.05.pdf", width=6, height=12) 
ggplot(Ds.siggeneFC.L, aes(x=logFC)) +geom_density()+
  theme_bw()+xlim(-6,6)+
  facet_grid(Treatment~.)
#dev.off() 
#2. Menidia
Mb.voomresults.adjp<-voomresults.adjp[,c(1,8:15)]
Mb.siggenes<-Mb.voomresults.adjp[apply(Mb.voomresults.adjp[, -1], MARGIN = 1, function(x) any(x <= 0.05)), ]
Mb.siggeneskey<-merge(voomresults.key,Mb.siggenes)
Mb.siggeneskey<-Mb.siggeneskey[,1:2]
d<-merge(CTMax2.0.HC.Mb,CTMax4.0.HC.Mb)
d<-merge(d,CTMax8.0.HC.Mb)
d<-merge(d,CTMax12.0.HC.Mb)
d<-merge(d,CTMax2.60.HC.Mb)
d<-merge(d,CTMax4.60.HC.Mb)
d<-merge(d,CTMax8.60.HC.Mb)
d<-merge(d,CTMax12.60.HC.Mb)
Mb.siggeneFC<-merge(Mb.siggeneskey,d)
#write.csv(Mb.siggeneskey,"orthogroup.adjp0.Mb.05.sig.genes.key.csv",row.names = F)

Mb.siggeneFC.L<-gather(Mb.siggeneFC, key = "Treatment", value = "logFC", CTMax2.0.HC.MblogFC,CTMax4.0.HC.MblogFC,CTMax8.0.HC.MblogFC,CTMax12.0.HC.MblogFC,CTMax2.60.HC.MblogFC,CTMax4.60.HC.MblogFC,CTMax8.60.HC.MblogFC,CTMax12.60.HC.MblogFC)

Mb.siggeneFC.L$Treatment<-factor(Mb.siggeneFC.L$Treatment)
levels(Mb.siggeneFC.L$Treatment)
Mb.siggeneFC.L$Treatment<-ordered(Mb.siggeneFC.L$Treatment,c("CTMax12.0.HC.MblogFC" , "CTMax12.60.HC.MblogFC" ,  "CTMax8.0.HC.MblogFC"  , "CTMax8.60.HC.MblogFC" ,"CTMax4.0.HC.MblogFC","CTMax4.60.HC.MblogFC" ,  "CTMax2.0.HC.MblogFC" , "CTMax2.60.HC.MblogFC"))

#pdf("orthogroup.dist.logFC.Mb.0.05.pdf", width=6, height=12) 
ggplot(Mb.siggeneFC.L, aes(x=logFC)) +geom_density()+
  theme_bw()+xlim(-6,6)+
  facet_grid(Treatment~.)
#dev.off() 
#################################################################
#Figure 5: Heatmap :####
#setup:
#reorder for heatmap visualization
siggeneFC1<-siggeneFC[,c(1,2,5,8,4,7,3,6,12,16,11,15,10,14,9,13)]
real.data2<-siggeneFC1
rnames <- real.data2[,1] 
mat_real.data2 <- data.matrix(real.data2[,3:ncol(real.data2)])
rownames(mat_real.data2) <- rnames

#to determine range and breaks, look over the distributions:
summary(siggeneFC[,3:16])
ggplot(siggeneFC.L, aes(x=Treatment, y=logFC)) + 
  geom_boxplot() + theme(axis.text.x=element_text(angle = -45, hjust = 0))
#min of all treatments is -6.6(Mb), -3.11(Ds);max is 7.5(Mb), 5(Ds)
#Looks like most within -4 to 4, some outliers especially  for Mb, unsuprisingly

#if have mostly normal distribution with one or a few outliers:
# #breaks for the core of the distribution:
# breaks=seq(-4.5,4.5, by=0.1)
# #add in for the outliers:
# breaks=append(breaks,7.6)
# breaks=append(breaks, -6.6,0) # the zero just puts this on the negative end where it is supposed to be
# breaks #check

#alternately can adjust further to shift, etc if you have skewed (and/vs. outlier) issues
breaks = c(seq(-6.6,-4,by=0.5),
           seq(-3.9999,4,by=0.1),
           seq(4.001,7.5,by=0.5))
breaks
#create colour panel with length(breaks)-1 colours
my_palette <- colorRampPalette(c("dodgerblue4","dodgerblue", "gray90","firebrick","firebrick4"))(n=length(breaks)-1)
#my_palette<-(redgreen(n=length(breaks)-1)) #can also just use preset color palettes for heatmaps

#create heatmap:####
#pdf("orthogroups.logFC.heatmap.padj0.05.pdf", width=6, height=8)
heatmap.2(mat_real.data2, 
  main = "Temp 0.05 Sig logFC", # heat map title
  notecol="black",# change font color of cell labels to black
  density.info="histogram",# turns on/off density plot in color legend
  trace="none",# turns off trace lines inside the heat map
  margins =c(8,8),# widens margins around plot
  col=my_palette,# use on color palette defined earlier 
  breaks=breaks,
  key=TRUE,#turn on/off color key
  symkey=FALSE,#don't know wtf this does, look up later
  scale="none",#normalize for each gene (turn off if you've done it previously):
  cexRow=1, cexCol=1,#text size for labels
  dendrogram="row",# only draw a row dendrogram
  Colv="NA"#,# turn off column clustering
 # lmat = rbind(c(0,3),c(2,1),c(0,4)), #layout: mat = lmat
#  lwid = c(0.5,7), #layout: widths = lwid
 # lhei = c(0.5,25,0
)#)
#dev.off()