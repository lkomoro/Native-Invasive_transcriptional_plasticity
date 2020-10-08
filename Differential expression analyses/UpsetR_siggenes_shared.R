#upsetR

library(tidyr)
require("UpSetR")
library(here)
here()
setwd(here("DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthogroup_DE_analyses/Temptime_splitanalysis"))
#################################################################################################

#Orthgroups
#Read in data & manipulate:####
voomresults<-read.delim("Orthogroup_annot_compiled_contrastslimma_rm2rds.txt")
adj<-grep("adj",names(voomresults))#extract #s of which columns I want
voomresults.adjp<-voomresults[,c(1,76:89)]#keeping just within species contrasts
grep("adj",names(voomresults.adjp),value=T)#double check correct columns
colnames(voomresults.adjp) <- c("Orthogroup","CTMax2.0.HC.Ds" ,  "CTMax4.0.HC.Ds" ,  "CTMax6.0.HC.Ds","CTMax2.60.HC.Ds" , "CTMax4.60.HC.Ds" , "CTMax6.60.HC.Ds" , "CTMax2.0.HC.Mb", "CTMax4.0.HC.Mb"  , "CTMax8.0.HC.Mb"  , "CTMax12.0.HC.Mb" , "CTMax2.60.HC.Mb", "CTMax4.60.HC.Mb",  "CTMax8.60.HC.Mb" , "CTMax12.60.HC.Mb")

s<-voomresults.adjp#just to be able to check the replacement worked
for (i in c(2:15))
{
  voomresults.adjp[,i] <- ifelse(voomresults.adjp[,i] <= 0.05, "1", "0")
  voomresults.adjp[,i] <- as.numeric(voomresults.adjp[,i])
}
pdf("../../Orthgroup_UpsetR/UpsetR_shared_sigorthogroups.pdf", width=30, height=12)#making super wide bc going to ax right side to only show sets with >=5
upset(voomresults.adjp, nsets = 14, nintersects = NA, mb.ratio = c(0.69, 0.31),sets=c("CTMax2.60.HC.Mb","CTMax2.0.HC.Mb",   "CTMax4.60.HC.Mb","CTMax4.0.HC.Mb", "CTMax8.60.HC.Mb","CTMax8.0.HC.Mb"  , "CTMax12.60.HC.Mb","CTMax12.0.HC.Mb" , "CTMax2.60.HC.Ds" , "CTMax2.0.HC.Ds" ,    "CTMax4.60.HC.Ds" , "CTMax4.0.HC.Ds" , "CTMax6.60.HC.Ds" ,  "CTMax6.0.HC.Ds"),
      order.by = "freq", decreasing = TRUE, number.angles = 30, point.size = 3.5, line.size = 1, keep.order=TRUE,
      mainbar.y.label = "No. Orthogroups Shared", sets.x.label = "Significant Orthogroups per Treatment")
dev.off() #can customize more attributes of plot with text size, order etc, see tutorials.
