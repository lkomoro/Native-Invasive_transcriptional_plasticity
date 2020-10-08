# Native-Invasive transcriptional plasticity
Repository for scripts and workflow notes for project investigating transcriptional plasticity in native vs invasive fishes

## Bioinformatics Analyses General Workflow 
*(sections align with subfolders in repository)*
#### 1. RNA-Seq Initial QC and Processing to counts
  * De-multiplexed and combined (concatenated) sequences across lanes for each sample
    * *NGS.nested.for.loop.ahi.sh*
  * FASTQC to assess initial data quality
  * Raw reads per sample filtered via adaptor trimming with Scythe and adaptive quality trimming with Sickle
    * *sickleLK.sh & scytheLK.sh*
  * Trimmed sequences aligned to previously generated reference transcriptomes for each species using BWA 
      * *BWA.mem.DS.ONLY-LK.sh, BWA.mem.Mb.ONLY-LK.sh,BWA.mem.singletons.DS.ONLY-LK.sh, BWA.mem.singletons.Mb.ONLY-LK.sh*
  * Additional filtering to remove secondary and supplementary alignments to avoid issues of paralog and chimeric sequences. 
      * *samtools_filter_rm2reads.sh*
  * Generate mapping stats
      * *samtools_flagstats.sh, flagstat_tabbed.sh* 
  * Resulting filtered alignments processed to counts per transcript contig using SAMtools
    * *bam_index_all.sh, samtools_idxstats.sh*
  * 
  
**Transdecoder and Orthofinder to link Ds and Mb together (see detailed HackMD, pull in key steps here to upload)**
  * Transdecoder to predit the translated peptide sequences within each transcriptome
    * *transdecoder_ORF_DS.sh, transdecoder_ORF_Mb.sh, transdecoder_predict_DS.sh, transdecoder_predict_Mb.sh*
  * Identified orthologues between the two transcriptomes using Orthofinder
    * *Orthofinder.sh*

  * Linked orthogroup output with the transcript count files for each species (using the filtered/rm2reads versions), roll up to orthogroup counts, combine to generaet one combined count matrix for DGE (including all samples for both species, transcripts summed within the same orthogroup)
    * *Orthologue_assignment.R*
    

#### 2. Differential Expression Analyses
  * Differential expression analyses using the limma-voom transformation and a general linear model 
   Conduct DGE, using group assignments that combine species, temp and time. need to state this in methods and explain rationale, caveats. output is in subfolder “Temptime_splitanalysis”
   
Limma_temptimesplit_orthogroup_rm2rds.R
Use output of MD figures and Venn diagrams from this for MS figures
Created graphs of distribution shifts and heatmap to use in MS

Sig_genes_heatmaps_orthogroups.R
Visualizations filtered to include orthogroups that were significantly different in at least one treatment
UpsetR to visualize shared and unique orthogroups between treatments

UpsetR_siggenes_shared.R
Ran functional analyses in TopGO; ran Fisher’s Exact Test and KS tests with both 0.05 and 0.01 thresholds (though the way KS works, setting a threshold shouldn’t matter), compared results and read about pros/cons of each and why appropriate for different datasets

  * Graphical visualization of differentially expressed genes 
    * Reaction Norms 
    * Heatmap 
    * UpSetR

#### 3. Functional Analyses
  * Pre-processing: combined and removed duplicate GO terms for annotated transcripts within orthogroups
  * Use GO information associated with each orthogroup to assess functional enrichment for a priori contrasts specified above using Fisher’s exact tests as implemented in topGO
  * Created similarity matrix using GoSemSim followed by hierarchical clustering for visualization of results
  
TopGO_Orthogroup.R
Decided to use results from Fisher’s 0.01 threshold to be conservative; though KS has advantages, in our dataset some treatments had few (or no) sig DE genes so seemed odd to have enriched biological categories. But, if reviewers have suggestions otherwise, can change, have other results saved.
Resources:
KS vs Fisher’s Exact
TopGO manual
Biostars discussion of Fisher’s Exact vs. KS
Helpful PP with background on GSEA tests
Used TopGO output from above with GoSemSim to generate heatmap of enriched biological processes. Output file GO_groupings_0.92clustering.csv will want to add as supplementary file, and then determined labels for groups based on examination of groups with parent GO terms, previous research (from old version of this figure-see file in other folder as needed). Worked up and finalized list of appropriate labels in excel: GO_groupings_0.92clustering_orthologsLMK.

Enriched_GO_Fisher_orthogroup_comparisons&analyses.R (script combining results across TopGO results to generate file including all GO terms for all treatments, and then apply a filter to include only GO terms that were enriched at p<0.01/0.05 for at least one treatment, write out to file that goes into the next script with GoSemSim, etc.)
Species_orthogroup_compare_GOheatmap_BP.R
