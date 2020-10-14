##Enriched GO Comparisons/Analyses####
#Reading in files that are output from TopGO analyses
#rm(list=ls())
library(here)
here()
setwd(here("./DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthogroup_TopGO_analyses"))

#Read in Fisher Exact Test Results:####
list.files(pattern="FisherGO.txt")

CTMax2.0.HC.Ds<-read.delim("CTMax2.0.HC.Ds_FisherGO.txt", header=T)
Key<-CTMax2.0.HC.Ds[,c(1:3,5)]

colnames(CTMax2.0.HC.Ds)[4] <- "CTMax2.0.HC.Dspval"
CTMax2.0.HC.Ds<-CTMax2.0.HC.Ds[,c(1,4)]

CTMax2.60.HC.Ds<-read.delim("CTMax2.60.HC.Ds_FisherGO.txt", header=T)
colnames(CTMax2.60.HC.Ds)[4] <- "CTMax2.60.HC.Dspval"
CTMax2.60.HC.Ds<-CTMax2.60.HC.Ds[,c(1,4)]

CTMax4.0.HC.Ds<-read.delim("CTMax4.0.HC.Ds_FisherGO.txt", header=T)
colnames(CTMax4.0.HC.Ds)[4] <- "CTMax4.0.HC.Dspval"
CTMax4.0.HC.Ds<-CTMax4.0.HC.Ds[,c(1,4)]

CTMax4.60.HC.Ds<-read.delim("CTMax4.60.HC.Ds_FisherGO.txt", header=T)
colnames(CTMax4.60.HC.Ds)[4] <- "CTMax4.60.HC.Dspval"
CTMax4.60.HC.Ds<-CTMax4.60.HC.Ds[,c(1,4)]

CTMax6.0.HC.Ds<-read.delim("CTMax6.0.HC.Ds_FisherGO.txt", header=T)
colnames(CTMax6.0.HC.Ds)[4] <- "CTMax6.0.HC.Dspval"
CTMax6.0.HC.Ds<-CTMax6.0.HC.Ds[,c(1,4)]

CTMax6.60.HC.Ds<-read.delim("CTMax6.60.HC.Ds_FisherGO.txt", header=T)
colnames(CTMax6.60.HC.Ds)[4] <- "CTMax6.60.HC.Dspval"
CTMax6.60.HC.Ds<-CTMax6.60.HC.Ds[,c(1,4)]

CTMax2.0.HC.Mb<-read.delim("CTMax2.0.HC.Mb_FisherGO.txt", header=T)
colnames(CTMax2.0.HC.Mb)[4] <- "CTMax2.0.HC.Mbpval"
CTMax2.0.HC.Mb<-CTMax2.0.HC.Mb[,c(1,4)]

CTMax2.60.HC.Mb<-read.delim("CTMax2.60.HC.Mb_FisherGO.txt", header=T)
colnames(CTMax2.60.HC.Mb)[4] <- "CTMax2.60.HC.Mbpval"
CTMax2.60.HC.Mb<-CTMax2.60.HC.Mb[,c(1,4)]

CTMax4.0.HC.Mb<-read.delim("CTMax4.0.HC.Mb_FisherGO.txt", header=T)
colnames(CTMax4.0.HC.Mb)[4] <- "CTMax4.0.HC.Mbpval"
CTMax4.0.HC.Mb<-CTMax4.0.HC.Mb[,c(1,4)]

CTMax4.60.HC.Mb<-read.delim("CTMax4.60.HC.Mb_FisherGO.txt", header=T)
colnames(CTMax4.60.HC.Mb)[4] <- "CTMax4.60.HC.Mbpval"
CTMax4.60.HC.Mb<-CTMax4.60.HC.Mb[,c(1,4)]

CTMax8.0.HC.Mb<-read.delim("CTMax8.0.HC.Mb_FisherGO.txt", header=T)
colnames(CTMax8.0.HC.Mb)[4] <- "CTMax8.0.HC.Mbpval"
CTMax8.0.HC.Mb<-CTMax8.0.HC.Mb[,c(1,4)]

CTMax8.60.HC.Mb<-read.delim("CTMax8.60.HC.Mb_FisherGO.txt", header=T)
colnames(CTMax8.60.HC.Mb)[4] <- "CTMax8.60.HC.Mbpval"
CTMax8.60.HC.Mb<-CTMax8.60.HC.Mb[,c(1,4)]

CTMax12.0.HC.Mb<-read.delim("CTMax12.0.HC.Mb_FisherGO.txt", header=T)
colnames(CTMax12.0.HC.Mb)[4] <- "CTMax12.0.HC.Mbpval"
CTMax12.0.HC.Mb<-CTMax12.0.HC.Mb[,c(1,4)]

CTMax12.60.HC.Mb<-read.delim("CTMax12.60.HC.Mb_FisherGO.txt", header=T)
colnames(CTMax12.60.HC.Mb)[4] <- "CTMax12.60.HC.Mbpval"
CTMax12.60.HC.Mb<-CTMax12.60.HC.Mb[,c(1,4)]

#Combine results: ####
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

#Subset & venn diagram####
#pval<-0.05 
pval<-0.01 #for stricter threshold
CTMax2.0.HC.Dssig<-subset(CTMax2.0.HC.Ds, CTMax2.0.HC.Dspval<=pval)
CTMax2.0.HC.Dssig<-CTMax2.0.HC.Dssig[,1]
#37 p0.01
#179 p0.05

CTMax4.0.HC.Dssig<-subset(CTMax4.0.HC.Ds, CTMax4.0.HC.Dspval<=pval)
CTMax4.0.HC.Dssig<-CTMax4.0.HC.Dssig[,1]
#7 p0.01
#50 p0.05

CTMax6.0.HC.Dssig<-subset(CTMax6.0.HC.Ds, CTMax6.0.HC.Dspval<=pval)
CTMax6.0.HC.Dssig<-CTMax6.0.HC.Dssig[,1]
#34 p0.01
#83 p0.05

CTMax2.60.HC.Dssig<-subset(CTMax2.60.HC.Ds, CTMax2.60.HC.Dspval<=pval)
CTMax2.60.HC.Dssig<-CTMax2.60.HC.Dssig[,1]
#70 p0.01
#187 p0.05

CTMax4.60.HC.Dssig<-subset(CTMax4.60.HC.Ds, CTMax4.60.HC.Dspval<=pval)
CTMax4.60.HC.Dssig<-CTMax4.60.HC.Dssig[,1]
#23 p0.01
#90 p0.05

CTMax6.60.HC.Dssig<-subset(CTMax6.60.HC.Ds, CTMax6.60.HC.Dspval<=pval)
CTMax6.60.HC.Dssig<-CTMax6.60.HC.Dssig[,1]
#0 p0.01
#0 p0.05

CTMax2.0.HC.Mbsig<-subset(CTMax2.0.HC.Mb, CTMax2.0.HC.Mbpval<=pval)
CTMax2.0.HC.Mbsig<-CTMax2.0.HC.Mbsig[,1]
#102 p0.01
#284 p0.05

CTMax2.60.HC.Mbsig<-subset(CTMax2.60.HC.Mb, CTMax2.60.HC.Mbpval<=pval)
CTMax2.60.HC.Mbsig<-CTMax2.60.HC.Mbsig[,1]
#79 p0.01
#304 p0.05

CTMax4.0.HC.Mbsig<-subset(CTMax4.0.HC.Mb, CTMax4.0.HC.Mbpval<=pval)
CTMax4.0.HC.Mbsig<-CTMax4.0.HC.Mbsig[,1]
#74 p0.01
#277 p0.05

CTMax4.60.HC.Mbsig<-subset(CTMax4.60.HC.Mb, CTMax4.60.HC.Mbpval<=pval)
CTMax4.60.HC.Mbsig<-CTMax4.60.HC.Mbsig[,1]
#80 p0.01
#269 p0.05

CTMax8.0.HC.Mbsig<-subset(CTMax8.0.HC.Mb, CTMax8.0.HC.Mbpval<=pval)
CTMax8.0.HC.Mbsig<-CTMax8.0.HC.Mbsig[,1]
#21 p0.01
#102 p0.05

CTMax8.60.HC.Mbsig<-subset(CTMax8.60.HC.Mb, CTMax8.60.HC.Mbpval<=pval)
CTMax8.60.HC.Mbsig<-CTMax8.60.HC.Mbsig[,1]
#56 p0.01
#263 p0.05

CTMax12.0.HC.Mbsig<-subset(CTMax12.0.HC.Mb, CTMax12.0.HC.Mbpval<=pval)
CTMax12.0.HC.Mbsig<-CTMax12.0.HC.Mbsig[,1]
#27 p0.01
#68 p0.05

CTMax12.60.HC.Mbsig<-subset(CTMax12.60.HC.Mb, CTMax12.60.HC.Mbpval<=pval)
CTMax12.60.HC.Mbsig<-CTMax12.60.HC.Mbsig[,1]
#7 p0.01
#79 p0.05


#use/put in apply function if want #s for table:####
length(intersect(CTMax2.0.HC.Dssig, CTMax2.60.HC.Dssig))#example
#length(setdiff(CTMax2.0.HC.Dssig, CTMax2.60.HC.Dssig))
#note this can only handle two at a time, and then have overlap need to account for (so see code below to generate common lists by temp between species)

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


###############################
t1<-comb[,c(1,5:18)]
t<-t1[apply(t1[, -1], MARGIN = 1, function(x) any(x <= pval)), ]#save only with at least one treatment at or lower than set threshold
comb1<-comb[,1:4]
comb2<-merge(comb1,t)
#write.csv(comb2,"Combined_orthogroup_sigGOenrichFisherp0.01.csv") #list including those only meeting threshold in one treatment
write.csv(comb,"Combined_orthogroup_enrichFisher_all.csv")
################################################
#Note that some treatments had low numbers of significant genes based on 0.05 adj FDR, and perhaps it doesn't make sense to include results (especially in heatmap)?
#CTMax4.0.Ds: 14 sig genes
#CTMax6.0.Ds: 7 sig genes
#CTMax6.60.Ds: 0 sig genes
#CTMax8.0.HC.Mb: 9 sig genes
#CTMax12.0.HC.Mb: 11 sig genes
#(see saved table with added column: orthogroups_summary_sig_genes_LK.xlsx)
###############################

#Pulling out specific shared GO terms between treatments:
#bc just exploring data, for now just viewing results and putting into written results in MS, not writing to file
pval<-0.01
pval<-0.05
#Import GO heatmap grouping IDs and names ####
library(dplyr)
g<-read.csv(here("./DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthogroup_GOheatmap_analyses/GO_SemSim_groups_0.92clustering.csv"))
h<-read.csv(here("./DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthogroup_GOheatmap_analyses/GO_SemSim_groupnames_0.92clustering.csv"))
i<-inner_join(g, h, by = "group")
j<-merge(i,comb2, all.y=T) 

write.table(j,"Combined_OG_sigGOenrichFisherp0.01_wheatmap_groupnames.txt",sep = "\t") #list including those only meeting threshold in one treatment, with heatmap IDs labeled


#Common relative temps:####
#CTMax-2####
t2<-t1[,c(1,2,5,8,12)]#all together
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) any(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-merge(i,t2.1)
#LK start here to look if there are shared groupings between species

t2<-t1[,c(1,2,8)]#0-0
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) all(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-t2.2[,1:2]
#shared between species in at least one time point (each) at 0.01 level:
#GO:0000122	negative regulation of transcription by RNA polymerase II
#GO:0007095	mitotic G2 DNA damage checkpoint
#GO:0009910	negative regulation of flower development

#shared between species in at least one time point (each) at 0.05 level(in addition to ones above listed for 0.01 level):
# GO:0001892	embryonic placenta development
# GO:0006379	mRNA cleavage
# GO:0009223	pyrimidine deoxyribonucleotide catabolic process
# GO:0016567	protein ubiquitination
# GO:0022008	neurogenesis
# GO:0033523	histone H2B ubiquitination
# GO:0033687	osteoblast proliferation
# GO:0035136	forelimb morphogenesis
# GO:0035666	TRIF-dependent toll-like receptor signaling pathway
# GO:0042254	ribosome biogenesis
# GO:0043984	histone H4-K16 acetylation
# GO:0048133	male germ-line stem cell asymmetric division
# GO:0048469	cell maturation
# GO:0090398	cellular senescence
# GO:1901216	positive regulation of neuron death

t2<-t1[,c(1,5,12)]#60-60
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) all(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-t2.2[,1:2]
#shared between species in at least one time point (each) at 0.01 level:
#GO:0006986	response to unfolded protein
#GO:0051726	regulation of cell cycle

#shared between species in at least one time point (each) at 0.05 level(in addition to ones above listed for 0.01 level):
# GO:0007595	lactation
# GO:0009719	response to endogenous stimulus
# GO:0010700	negative regulation of norepinephrine secretion
# GO:0034097	response to cytokine
# GO:0035264	multicellular organism growth
# GO:0035767	endothelial cell chemotaxis
# GO:0035914	skeletal muscle cell differentiation
# GO:0040026	positive regulation of vulval development
# GO:0045746	negative regulation of Notch signaling pathway
# GO:0051085	chaperone cofactor-dependent protein refolding

t2<-t1[,c(1,2,12)]#DS0-Mb60
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) all(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-t2.2[,1:2]
#shared between species in at least one time point (each) at 0.01 level:
#GO:0000122	negative regulation of transcription by RNA polymerase II

#shared between species in at least one time point (each) at 0.05 level(in addition to ones above listed for 0.01 level):
#GO:0006325	chromatin organization
#GO:0006513	protein monoubiquitination
#GO:0016567	protein ubiquitination
#GO:0030910	olfactory placode formation
#GO:0034138	toll-like receptor 3 signaling pathway
#GO:0035264	multicellular organism growth
#GO:0035666	TRIF-dependent toll-like receptor signaling pathway
#GO:0043281	regulation of cysteine-type endopeptidase activity involved in apoptotic process
#GO:0045666	positive regulation of neuron differentiation
#GO:0046133	pyrimidine ribonucleoside catabolic process
#GO:0048208	COPII vesicle coating
#GO:0048469	cell maturation
#GO:0048852	diencephalon morphogenesis
#GO:0090398	cellular senescence
#GO:1900429	negative regulation of filamentous growth of a population of unicellular organisms

t2<-t1[,c(1,5,8)]#DS60-Mb0
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) all(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-t2.2[,1:2]
#shared between species in at least one time point (each) at 0.01 level:
#GO:0045600	positive regulation of fat cell differentiation
#GO:0051726	regulation of cell cycle

#shared between species in at least one time point (each) at 0.05 level(in addition to ones above listed for 0.01 level):
#GO:0007623	circadian rhythm
#GO:0019371	cyclooxygenase pathway
#GO:0021612	facial nerve structural organization
#GO:0021658	rhombomere 3 morphogenesis
#GO:0035767	endothelial cell chemotaxis
#GO:0035914	skeletal muscle cell differentiation
#GO:0035924	cellular response to vascular endothelial growth factor stimulus
#GO:0040026	positive regulation of vulval development
#GO:0045648	positive regulation of erythrocyte differentiation
#GO:0046697	decidualization
#GO:0048025	negative regulation of mRNA splicing, via spliceosome
#GO:0050872	white fat cell differentiation
#GO:0051085	chaperone cofactor-dependent protein refolding
#GO:0070935	3'-UTR-mediated mRNA stabilization
#GO:1900119	positive regulation of execution phase of apoptosis

#CTMax-4####
t2<-t1[,c(1,3,6,9,13)]#all together
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) any(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-merge(i,t2.1)

t2<-t1[,c(1,3,9)]#0-0
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) all(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-t2.2[,1:2]
#shared between species in at least one time point (each) at 0.01 level:NONE
#shared between species in at least one time point (each) at 0.05 level(in addition to ones above listed for 0.01 level):
#GO:0051085	chaperone cofactor-dependent protein refolding

t2<-t1[,c(1,6,13)]#60-60
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) all(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-t2.2[,1:2]
#shared between species in at least one time point (each) at 0.01 level:
#GO:0006986	response to unfolded protein

#shared between species in at least one time point (each) at 0.05 level(in addition to ones above listed for 0.01 level):
# GO:0001892	embryonic placenta development
# GO:0001938	positive regulation of endothelial cell proliferation
# GO:0006686	sphingomyelin biosynthetic process
# GO:0008380	RNA splicing
# GO:0010114	response to red light
# GO:0021575	hindbrain morphogenesis
# GO:0043401	steroid hormone mediated signaling pathway
# GO:0043434	response to peptide hormone
# GO:0045740	positive regulation of DNA replication
# GO:0051726	regulation of cell cycle
# GO:1900119	positive regulation of execution phase of apoptosis
# GO:2001240	negative regulation of extrinsic apoptotic signaling pathway in absence of ligand

t2<-t1[,c(1,3,13)]#DS0-Mb60
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) all(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-t2.2[,1:2]
#shared between species in at least one time point (each) at 0.01 level:
#GO:0043401	steroid hormone mediated signaling pathway

#shared between species in at least one time point (each) at 0.05 level(in addition to ones above listed for 0.01 level):
# GO:0042997	negative regulation of Golgi to plasma membrane protein transport
# GO:0070370	cellular heat acclimation


t2<-t1[,c(1,6,9)]#DS60-Mb0
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) all(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-t2.2[,1:2]
#shared between species in at least one time point (each) at 0.01 level:
#GO:0006986	response to unfolded protein
#GO:0031048	chromatin silencing by small RNA

#shared between species in at least one time point (each) at 0.05 level(in addition to ones above listed for 0.01 level):
#GO:0001892	embryonic placenta development
#GO:0001938	positive regulation of endothelial cell proliferation
#GO:0008380	RNA splicing
#GO:0009162	deoxyribonucleoside monophosphate metabolic process
#GO:0009223	pyrimidine deoxyribonucleotide catabolic process
#GO:0017143	insecticide metabolic process
#GO:0045648	positive regulation of erythrocyte differentiation
#GO:0051085	chaperone cofactor-dependent protein refolding
#GO:0051726	regulation of cell cycle

#CTMax-6/8####
t2<-t1[,c(1,4,7,10,14)]#all together
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) any(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-merge(i,t2.1)

t2<-t1[,c(1,4,10)]#0-0
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) all(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-t2.2[,1:2]
#shared between species in at least one time point (each) at 0.01 level:NONE
#shared between species in at least one time point (each) at 0.05 level(in addition to ones above listed for 0.01 level):
#GO:0040040	thermosensory behavior

t2<-t1[,c(1,7,14)]#60-60
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) all(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-t2.2[,1:2]
#shared between species in at least one time point (each) at 0.01 level:NONE
#shared between species in at least one time point (each) at 0.05 level(in addition to ones above listed for 0.01 level):NONE

t2<-t1[,c(1,4,14)]#DS0-Mb60
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) all(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-t2.2[,1:2]
#shared between species in at least one time point (each) at 0.01 level:NONE
#shared between species in at least one time point (each) at 0.05 level(in addition to ones above listed for 0.01 level):
#GO:0000122	negative regulation of transcription by RNA polymerase II
#GO:0043374	CD8-positive, alpha-beta T cell differentiation
#GO:0060017	parathyroid gland development

t2<-t1[,c(1,7,10)]#DS60-Mb0
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) all(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-t2.2[,1:2]
#shared between species in at least one time point (each) at 0.01 level:NONE
#shared between species in at least one time point (each) at 0.05 level(in addition to ones above listed for 0.01 level):NONE


#Common absolute temps:####
#22C:
t2<-t1[,c(1,4,7,11,15)]#all together
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) any(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-merge(i,t2.1)

#shared between species in at least one time point (each) at 0.01 level:
#NONE

#shared between species in at least one time point (each) at 0.05 level:
#GO:0070370 cellular heat acclimation 
#GO:0006825 copper ion transport 	

#26C:
t2<-t1[,c(1,2,5,10,14)]#all together
t2.1<-t2[apply(t2[, -1], MARGIN = 1, function(x) any(x <= pval)), ]
t2.2<-merge(comb1,t2.1)
t2.3<-merge(i,t2.1)
#shared between species in at least one time point (each) at 0.01 level:
#GO:0007623	circadian rhythm

#shared between species in at least one time point (each) at 0.05 level:
#GO:0000122	negative regulation of transcription by RNA polymerase II
#GO:0021781 glial cell fate commitment 	
#GO:0002052 positive regulation of neuroblast proliferation 
#GO:0048852	diencephalon morphogenesis 
#GO:0030910	olfactory placode formation
#GO:0043984	histone H4-K16 acetylation
#GO:0001892	embryonic placenta development
#GO:0016567	protein ubiquitination
#GO:0035666	TRIF-dependent toll-like receptor signaling pathway
#GO:0033687	osteoblast proliferation
#GO:0034138	toll-like receptor 3 signaling pathway
#GO:0051412	response to corticosterone
#GO:0046697	decidualization
#GO:0032227	negative regulation of synaptic transmission, dopaminergic
#GO:0032811	negative regulation of epinephrine secretion
#GO:0046426	negative regulation of JAK-STAT cascade
#GO:0035914	skeletal muscle cell differentiation
#GO:0072525	pyridine-containing compound biosynthetic process
#GO:0034097	response to cytokine

#SCRATCH TO QUICKLY LOOK AT GROUPS FOR EACH TREATMENT
t2<-t1[,c(1,15)]#change column as desired here
names(t2)[2]
names(t2)[2] <- "pval"
t2.1<-subset(t2, pval<=.01)
t2.3<-merge(i,t2.1)
t2.3$group<-factor(t2.3$group)
unique(t2.3$group)
