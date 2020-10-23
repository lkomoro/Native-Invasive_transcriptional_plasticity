#This is using GoSemSim and supporting packages, adapted based on script from Kyle Wellband-see his email/'Goby_GO_heatmap.R' script and Evol Apps 2017 paper for details
#Additional resources:
#https://bioconductor.org/packages/devel/bioc/vignettes/GOSemSim/inst/doc/GOSemSim.html
#http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
library(here)
setwd("./Functional_Enrichment_analyses/")

# Setup environment
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("GOSemSim")
#BiocManager::install("ComplexHeatmap")
#BiocManager::install("circlize")
library(GOSemSim)
library(dendextend)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(tidyr)

#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)#human database


# Read the list of GO terms to include in the heatmap
#results from Fisher's exact test, stricter threshold of p<=0.01
xx<-read.csv("Combined_orthogroup_sigGOenrichFisherp0.01.csv")

xx.1<-as.data.frame(xx$GO.ID)
# Expand all pairwise comparisons of unique GO terms; in this case using biological process 
xx.1 <- unique(as.character(xx.1[,1]))

# Extract human GO info for information content (IC)
goHs <- godata(OrgDb = org.Hs.eg.db, ont = "BP")

#NB- our GO terms were from fish blasts at the order level (see Monica's emails), and many aren't showing up in human or zebrafish. But, there are limited species that have the GO annotation databases

# Calculate similarity among GO
sim <- goSim(xx.1, xx.1, semData = goHs, measure = "Rel")
sum(is.na(sim))

# convert to a distance measure and coerce to a square matrix
mat <- matrix(c(1-sim), ncol = length(xx.1), nrow = length(xx.1)) 
colnames(mat) <- rownames(mat) <- xx.1

mat[!is.finite(mat)] <- 0 #replace NAs with 0s
#NB this is not ideal, but limited currently by databases and goal here is to look comparatively across the species so will achieve that.

# convert square matrix to a distance matrix and cluster
dist <- as.dist(mat)
clust <- hclust(dist)

# From Kyle: "you will want to play with the value for cutting the dendrogram based on the depth of the groupings you want on the y-axis. I used the plotting functions to see which clusters would result with different values of "h" below"
#pdf(file = "GO_cut_threshold_plot.pdf", width = 6, height = 16)
plot(clust, hang = -1, cex = 0.2)+
abline(h = 0.92, lty = 2)+
abline(h = 0.9, lty = 2)+
abline(h = 0.95, lty = 2)+
abline(h = 0.85, lty = 2)
#dev.off()

# after looking at plot above, convert to a dendrogram and cut branches at differece = 0.92 to create representative but managable groupings
dend <- as.dendrogram(clust)
group <- cutree(clust, h = 0.92)
group_ordered <- group[order.dendrogram(dend)]

# pull in GO enrichment results
yy<-xx[,c(2,6:19)]
str(yy)
yy.1<-gather(yy,key="Treatment",value="pval",CTMax6.0.HC.Dspval,CTMax6.60.HC.Dspval,CTMax4.0.HC.Dspval,CTMax4.60.HC.Dspval,CTMax2.0.HC.Dspval,CTMax2.60.HC.Dspval,CTMax12.0.HC.Mbpval,CTMax12.60.HC.Mbpval,CTMax8.0.HC.Mbpval,CTMax8.60.HC.Mbpval,CTMax4.0.HC.Mbpval,CTMax4.60.HC.Mbpval,CTMax2.0.HC.Mbpval,CTMax2.60.HC.Mbpval)

yy.1[is.na(yy.1)] <- 1 #get rid of any NAs (just in case)
yy.1$Treatment<-gsub("pval","",yy.1$Treatment)#just getting rid of extra crap in text
yy.1$Treatment<-gsub(".HC","",yy.1$Treatment)
yy.1$Treatment<-factor(yy.1$Treatment, levels=c("CTMax6.0.Ds","CTMax6.60.Ds","CTMax4.0.Ds","CTMax4.60.Ds", "CTMax2.0.Ds","CTMax2.60.Ds","CTMax12.0.Mb" ,"CTMax12.60.Mb","CTMax8.0.Mb","CTMax8.60.Mb","CTMax4.0.Mb","CTMax4.60.Mb", "CTMax2.0.Mb" ,"CTMax2.60.Mb"))
 
# log transform the pvalue for plotting if desired
yy.1$value <- log(yy.1$pval)
str(yy.1)
hist(yy.1$value)
ggplot(yy.1, aes(x=(value))) +geom_density()
hist(yy.1$pval)
ggplot(yy.1, aes(x=(pval))) +geom_density()
# change the sign of the GO terms from up-regulated genes to facilitate using a divergent color palette in the heatmap
#yy[which(yy$Direction == "Up"), "value"] <- yy[which(yy$Direction == "Up"), "value"] * -1

# build a dataframe by treatment for heatmap plotting
yy_l <- split(yy.1, f = yy.1$Treatment)
zz <- data.frame(row.names = xx.1,
                 CTMax6.0.Ds = numeric(length = length(xx.1)),
                 CTMax6.60.Ds = numeric(length = length(xx.1)),
                 CTMax4.0.Ds = numeric(length = length(xx.1)),
                 CTMax4.60.Ds = numeric(length = length(xx.1)),
                 CTMax2.0.Ds = numeric(length = length(xx.1)),
                 CTMax2.60.Ds = numeric(length = length(xx.1)),
                 CTMax12.0.Mb = numeric(length = length(xx.1)),
                 CTMax12.60.Mb = numeric(length = length(xx.1)),
                 CTMax8.0.Mb = numeric(length = length(xx.1)),
                 CTMax8.60.Mb = numeric(length = length(xx.1)),
                 CTMax4.0.Mb = numeric(length = length(xx.1)),
                 CTMax4.60.Mb = numeric(length = length(xx.1)),
                 CTMax2.0.Mb = numeric(length = length(xx.1)),
                 CTMax2.60.Mb = numeric(length = length(xx.1))
                 )
zz[as.character(yy_l$CTMax2.0.Ds$GO.ID),"CTMax2.0.Ds"] <- yy_l$CTMax2.0.Ds$pval
zz[as.character(yy_l$CTMax4.0.Ds$GO.ID),"CTMax4.0.Ds"] <- yy_l$CTMax4.0.Ds$pval
zz[as.character(yy_l$CTMax6.0.Ds$GO.ID),"CTMax6.0.Ds"] <- yy_l$CTMax6.0.Ds$pval
zz[as.character(yy_l$CTMax2.60.Ds$GO.ID),"CTMax2.60.Ds"] <- yy_l$CTMax2.60.Ds$pval
zz[as.character(yy_l$CTMax4.60.Ds$GO.ID),"CTMax4.60.Ds"] <- yy_l$CTMax4.60.Ds$pval
zz[as.character(yy_l$CTMax6.60.Ds$GO.ID),"CTMax6.60.Ds"] <- yy_l$CTMax6.60.Ds$pval
zz[as.character(yy_l$CTMax2.0.Mb$GO.ID),"CTMax2.0.Mb"] <- yy_l$CTMax2.0.Mb$pval
zz[as.character(yy_l$CTMax4.0.Mb$GO.ID),"CTMax4.0.Mb"] <- yy_l$CTMax4.0.Mb$pval
zz[as.character(yy_l$CTMax8.0.Mb$GO.ID),"CTMax8.0.Mb"] <- yy_l$CTMax8.0.Mb$pval
zz[as.character(yy_l$CTMax12.0.Mb$GO.ID),"CTMax12.0.Mb"] <- yy_l$CTMax12.0.Mb$pval
zz[as.character(yy_l$CTMax2.60.Mb$GO.ID),"CTMax2.60.Mb"] <- yy_l$CTMax2.60.Mb$pval
zz[as.character(yy_l$CTMax4.60.Mb$GO.ID),"CTMax4.60.Mb"] <- yy_l$CTMax4.60.Mb$pval
zz[as.character(yy_l$CTMax8.60.Mb$GO.ID),"CTMax8.60.Mb"] <- yy_l$CTMax8.60.Mb$pval
zz[as.character(yy_l$CTMax12.60.Mb$GO.ID),"CTMax12.60.Mb"] <- yy_l$CTMax12.60.Mb$pval

# create heatmap column annotation bars
df = data.frame(Species = c(rep("Delta smelt", each = 6),rep("Menidia", each = 8)))
ha = HeatmapAnnotation(df = df, col = list(Species = c("Menidia" = "blue", "Delta smelt" = "firebrick")))

# create heatmap row annotation bars
df_r = data.frame(Group = as.character(group))
even = seq(from = 2, to = length(unique(group)), by = 2)
odd = seq(from = 1, to = length(unique(group)), by = 2)
l = list(Group = unique(group))
l[["Group"]][unique(group_ordered)[odd]] = "black"
l[["Group"]][unique(group_ordered)[even]] = "darkgrey"
names(l[["Group"]]) = unique(group)
ha_r = HeatmapAnnotation(df = df_r, col = l, which = "row", show_legend = F, width = unit(2, "mm"))

my_palette <- colorRamp2(c(0,0.1), c("forestgreen", "grey80"))
# create heatmap
hm = Heatmap(as.matrix(zz),
             col = my_palette,
             # passing the dendrogram from above to group the y-axis here
             cluster_rows = as.dendrogram(clust), 
             cluster_columns = F,
             show_row_names = F,
            # row_names_gp = gpar(fontsize = 2),
            # row_names_side = "left",
             show_column_names = T,
             top_annotation = ha,
            # split = length(unique(group)),
             name = "pvalue",
             show_heatmap_legend = T
)

# export plot as a PDF
pdf(file = "GO_heatmap_0.92clust.pdf", width = 6, height = 16)
hm + ha_r
dev.off()

#cross check/extract GO term groupings
hm1 = Heatmap(as.matrix(zz),
             col = my_palette,
             # passing the dendrogram from above to group the y-axis here
             cluster_rows = as.dendrogram(clust), 
             cluster_columns = F,
             show_row_names = T,
              row_names_gp = gpar(fontsize = 2),
              row_names_side = "right",
             show_column_names = T,
             top_annotation = ha,
             # split = length(unique(group)),
             name = "pvalue",
             show_heatmap_legend = T
)


# export plot as a PDF
pdf(file = "GO_ordercheck_0.92clust.pdf", width = 6, height = 16)
hm1 # + ha_r  for some reason now if I combine, the GO terms don't print anymore so just leaving off since just for checking for group assignments etc see Adobe Photoshop merged check doc
dev.off()

#check GO groupings
g1<-as.data.frame(group)
g1$GO.ID<-row.names(g1)
a<-xx[,c(2,3)]
b<-merge(g1,a)
write.csv(b, "GO_groupings_0.92clustering.csv",row.names = F)
