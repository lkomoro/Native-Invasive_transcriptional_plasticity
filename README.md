# Native-Invasive transcriptional plasticity
Repository for scripts and workflow notes for project investigating transcriptional plasticity in native vs invasive fishes

## Bioinformatics Analyses General Workflow 
*(sections align with subfolders in repository)*
#### 1. RNA-Seq Initial QC and Processing to counts
  * De-multiplexed and combined (concatenated) sequences across lanes for each sample
  * FASTQC to assess initial data quality
  * Raw reads per sample filtered via adaptor trimming with Scythe and adaptive quality trimming with Sickle
  * Trimmed sequences aligned to previously generated reference transcriptomes for each species using BWA 
      * Additional filtering to remove secondary and supplementary alignments to avoid issues of paralog and chimeric sequences. 
  * Resulting filtered alignments processed to counts per transcript contig using SAMtools idxstats 
  * Transdecoder to predit the translated peptide sequences within each transcriptome 
  * Identified orthologues between the two transcriptomes using Orthofinder
  * Generated a combined count matrix including all samples for both species, transcripts summed within the same orthogroup. 

#### 2. Differential Expression Analyses
  * Differential expression analyses using the limma-voom transformation and a general linear model 
  * Graphical visualization of differentially expressed genes 
    * Reaction Norms 
    * Heatmap 
    * UpSetR

#### 3. Functional Analyses
  * Pre-processing: combined and removed duplicate GO terms for annotated transcripts within orthogroups
  * Use GO information associated with each orthogroup to assess functional enrichment for a priori contrasts specified above using Fisher’s exact tests as implemented in topGO
  * Created similarity matrix using GoSemSim followed by hierarchical clustering for visualization of results
