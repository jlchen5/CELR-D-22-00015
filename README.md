# CELR-D-22-00015

- Repository that contains the analysis from Li et al. ***Cell Regeneration***, 2022 (*under review*).  

## I. Upstream analysis

- [RNA-seq_analysis.sh](https://github.com/jlchen5/CELR-D-22-00015/blob/main/RNA-seq_analysis.sh)    
RNA-seq data is applied quality and adapter trimming by [Trim Galore v0.6.4_dev](https://zenodo.org/badge/latestdoi/62039322) and mapped to the [mm9 reference genome](http://genome.ucsc.edu/cgi-bin/hgGateway?clade=mammal&org=Mouse&db=mm9) using [STAR v2.7.1a](https://github.com/alexdobin/STAR/releases/tag/2.7.1a). And we performed reads summarization for genomic features via [featureCounts v2.0.0](http://subread.sourceforge.net/featureCounts.html).

- [ChIP-seq_analysis.sh](https://github.com/jlchen5/CELR-D-22-00015/blob/main/ChIP-seq_analysis.sh)  
Adapters and low-quality bases were removed using [Trim Galore v0.6.4_dev](https://zenodo.org/badge/latestdoi/62039322). ChIP-seq reads were mapped to the [mm9 reference genome](http://genome.ucsc.edu/cgi-bin/hgGateway?clade=mammal&org=Mouse&db=mm9) using [Bowtie2 v2.3.5.1](https://github.com/BenLangmead/bowtie2/releases/tag/v2.3.5.1). [Samtools v1.10](https://github.com/samtools/samtools/releases/tag/1.10) was used to filter and convert file formats. To remove PCR duplicate reads, which might occur false positive results, [MarkDuplicates (v4.1.4.1, Picard)](https://gatk.broadinstitute.org/hc/en-us/sections/360007458971-4-1-4-1) was used to locate and tag duplicate reads. And sequencing duplicates were removed using the `--REMOVE_SEQUENCING_DUPLICATES` options. ChIP-seq peaks were called using [MACS2 v2.2.6](https://github.com/macs3-project/MACS/releases/tag/v2.2.6).  

## II. Downstream analysis

- The raw read counts normalization and differential expression were determined by [RNA-seq_rif1.R](https://github.com/jlchen5/CELR-D-22-00015/blob/main/RNA-seq_rif1.R).
 
  > [coding_genes.R](https://github.com/jlchen5/CELR-D-22-00015/blob/main/coding_genes.R) is for further analysis of genes  
  > [repeats.R](https://github.com/jlchen5/CELR-D-22-00015/blob/main/repeats.R) is for further analysis of repeats

- Correlation among different samples is analyzed via [CalculateCorrelation.R](https://github.com/jlchen5/CELR-D-22-00015/blob/main/CalculateCorrelation.R) and [plotCorrHeatmap.R](https://github.com/jlchen5/CELR-D-22-00015/blob/main/plotCorrHeatmap.R) is used to plot correlation heatmap.
- Principal Component Analysis (PCA) is calculated the foldchange (vs. control/wild type) between the samples with the depletion of the indicated genes to control (or wild type) samples via [PCAanalysis.R](https://github.com/jlchen5/CELR-D-22-00015/blob/main/PCAanalysis.R).
- [GOanalysis.R](https://github.com/jlchen5/CELR-D-22-00015/blob/main/GOanalysis.R) is used to GO analysis.


## III. Quality Control

All data in this article has been quality controlled:   
 - [ChIP-seq_multiqc_report](https://github.com/jlchen5/CELR-D-22-00015/blob/main/QC/ChIP-seq_multiqc_report.html)  
 - [RNA-seq_multiqc_report](https://github.com/jlchen5/CELR-D-22-00015/blob/main/QC/RNA-seq_multiqc_report.html)  
 - [Repli-seq_multiqc_report](https://github.com/jlchen5/CELR-D-22-00015/blob/main/QC/repli-seq_multiqc_report.html)
