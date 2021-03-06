################################## GO Enrichment Analysis with TopGO #########################################
# Using TopGO package for functional analyses, see Alexa et al. 2006

#Additional references for this package:
#http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html
#https://bioconductor.org/packages/release/bioc/manuals/topGO/man/topGO.pdf
#https://www.biostars.org/p/247636/
#https://stackoverflow.com/questions/12054816/r-finding-duplicates-in-one-column-and-collapsing-in-a-second-column (for ragged array conversion code options)
#browseVignettes("topGO")
############################################################################################################
#I. Setup####
#rm(list=ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("topGO")
#BiocManager::install("Rgraphviz")
#BiocManager::install("GO.db")
library(Rgraphviz)
library(topGO)
library(plyr)
library(here)
here()

#N.B. Before starting also want to check that have most recent version of Bioconductor, topGO and GO.db, since this is where the package pulls the info for the GO terms from
#check if needed:
sessionInfo()

if(!(file.exists("Orthogroup_Gene2GO.map"))){
  anno<-read.csv("./Functional_Enrichment_analyses/GO terms of Shared Orthogroups by Transcript_nodups.csv", header=T)
	anno2 <- anno[,c("Orthogroup", "GO_terms")]
	anno3<-ddply(anno2, "Orthogroup", summarize, GO_terms = paste(GO_terms, collapse=","))
	anno3<-rename(anno3,c("Orthogroup"="Gene","GO_terms"="GO_ID"))#just to make rest script match for now
	write.table(anno3,  file ="./Functional_Enrichment_analyses/Orthogroup_Gene2GO.map", sep = "\t", quote = F,
			row.names = F, col.names = F)
}

geneID2GO <- readMappings(file = ("./Functional_Enrichment_analyses/Orthogroup_Gene2GO.map"))
infiles <- list.files("./Functional_Enrichment_analyses/voom_results",pattern = "voom_results", full.names = T)#these are the output from writing the topTable results to file, need to include all results so have complete gene ('geneuniverse') set to draw from below

##########################################################################
#Adjust file location paths so write to correct subfolders as needed
#II. Fisher's Exact Testing####
for (file in infiles){
  tmp <- read.delim(file, stringsAsFactors = F)
  colnames(tmp)[1] <- "Gene"
  DE <- tmp
  pcutoff <- 0.05 # cutoff for defining significant genes
  tmp <- ifelse(DE$adj.P.Val < pcutoff, 1, 0)# Define gene list as 1's if adjP < cutoff, 0, otherwise
  geneList <- tmp
  names(geneList) <- unlist(DE$Gene, function(x)x[1])
 # names(geneList) <- unlist(lapply(strsplit(DE$Gene, split = ".", fixed = T), function(x)x[1]))   # geneList needs names that match those for GO terms, so stripping off part after dot again if have that format
  head(geneList) 
  # Create topGOData object
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,#this is your 'gene universe'
                geneSelectionFun = function(x)(x == 1), #genes of interest to use are the ones defined above 
                nodeSize=5,#optional argument-this is used to prune the GO hierarchy, eg. nodesize=5 prunes the GO hierarchy, to remove terms which have less than 5 annotated genes.
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
#Also can add 'description="My project",' argument: The 'description' argument has a short description of your project.
#get a summary of the data in the object:
GOdata

#The number of GO terms in the TopGO subset of the GO hierarchy (the GO terms annotated to genes in the gene 'universe' input file, plus the ancestors of those GO terms in the GO hierarchy) can be found using:
length(usedGO(GOdata))
#This should agree with the 'number of nodes' given for the GO graph when you type 'GOdata', as it is the number of nodes in TopGO's internal version of the GO hierarchy (which is filtered just to contain your terms annotated to your gene universe, plus their ancestors). 

#Note that if you use the weight01 (default), elim or parentchild algorithm, TopGO may report a GO term as significantly enriched which wasn't actually in your input annotations, but rather is an ancestor (in the GO hierarchy) of a GO term in your input annotations.

#	The list of genes of interest can be accessed using the method sigGenes():
sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)
#this should match what you know you put in from the DE results-run as a cross check that your code above worked to ID correctly at the FDR threshold

# Fisher testing
resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
#you can choose different algorithms and statistics here. See the manual and Alexa et al. paper for options and explanations

# get a summary of the results:
resultFisher #(prints to screen-set nohup.out file to capture/run from commandline)

#List the results found:
tab <- GenTable(GOdata, Fisher = resultFisher, topNodes = length(resultFisher@score),numChar = 120)#change length of topNodes if just want top 100 sig results, etc

t<-head(tab)#to see what results format look like
  tab <- tab[,c(1,2,3,6)]
  names(tab)[4] <- "Raw.P.Value"
  #tab$adj.P.Value <- p.adjust(tab$Raw.P.Value) #Note-The topGO authors recommend against using adjusted p-values for GO enrichment due to the high amount of correlation/nesting among terms.
  genes.in.term <- genesInTerm(GOdata)
  genes.in.term2 <- unlist(lapply(genes.in.term, function(x)paste(x, collapse = ", ")))
  genes.in.term2 <- cbind(names(genes.in.term), genes.in.term2)#these lines retrieve and convert the genes that have each GO term to then combine with results below 
  tmp <- match(tab$GO.ID, genes.in.term2[,1])
  tab$Genes.In.Term <- genes.in.term2[tmp,2]
  file2<-gsub(".txt","",file)
  file3<-gsub("./Functional_Enrichment_analyses/voom_results_","",file2)
  outfilename <- paste(file3,"_FisherGO.txt", sep = "")
  write.table(tab, file = outfilename, row.names = F, quote = F, sep = "\t")

# # Visualize the position of the statistically significant GO terms in the GO hierarchy using the following function:
#   showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, 
#                  useInfo ='all',#additional info to be plotted to each node
#                  reverse=T, #graph it upside down or not; default is T
#                  wantedNodes = NULL, #can specify nodes you want to be able to find and they will be colored so you can see them easily
#                  type=NULL,#can change to plot pie charts supposedly by no options are listed to specify this
#                  showEdges=T,#default is T; if F removed arrow showing relationships
#                  sigForAll=T #score/p-value of all nodes in the DAG is shown, otherwise only the score for the sigNodes; default is T
#                   )
# # The significant GO terms are shown as rectangles in the picture. The most significant terms are coloured red and least significant in yellow:
#I'm creating 5 & 10 just compare, adjust and reduce to smallest node cutoff that accurately depicts main relationships once get feel for data in different datasets
}
#################################################################
