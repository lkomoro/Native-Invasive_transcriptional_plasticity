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
  * Transdecoder to predit the translated peptide sequences within each transcriptome *N.B. detailed HackMD document has additional info of output etc*
    * *transdecoder_ORF_DS.sh, transdecoder_ORF_Mb.sh, transdecoder_predict_DS.sh, transdecoder_predict_Mb.sh*
  * Identified orthologues between the two transcriptomes using Orthofinder
    * *Orthofinder.sh*
  * Linked orthogroup output with the transcript count files for each species (filtered/rm2reads versions), roll up to orthogroup counts, combine to generate one combined count matrix for DGE (including all samples for both species, transcripts summed within the same orthogroup)
    * *Orthologue_assignment.R*
    
#### 2. Differential Expression Analyses
  * Differential expression analyses using the limma-voom transformation and a general linear model 
   * *Limma_temptimesplit_orthogroup_rm2rds.R*
 * Created graphs of LFC distribution shifts, DE reaction norms and heatmap for MS
   * *Sig_genes_heatmaps_orthogroups.R, Reaction_normgraph.R*
 * UpsetR to visualize shared and unique DE orthogroups between treatments
   * *UpsetR_siggenes_shared.R*

#### 3. Functional Analyses
 * Pre-processing: combined and removed duplicate GO terms for annotated transcripts within orthogroups
 * Assessed functional enrichment for a priori contrasts using Fisher’s exact tests as implemented in TopGO
    * *TopGO_Orthogroup.R*
  * Created similarity matrix using GoSemSim followed by hierarchical clustering for visualization of enriched biological processes, labels from Revigo and other database sources *(manual investigation for each subgroup, finalized list of appropriate labels in excel: GO_groupings_0.92clustering_orthologsLMK).*
    * *Enriched_GO_Fisher_orthogroup_comparisons&analyses.R, Species_orthogroup_compare_GOheatmap_BP.R
