# CELR-D-22-00015

[![DOI](https://img.shields.io/badge/DOI-10.1186%2Fs13619--022--00124--9-blue)](https://doi.org/10.1186/s13619-022-00124-9)


- Repository that contains the analysis from [Li et al. ***Cell Regeneration***, 2022](https://doi.org/10.1186/s13619-022-00124-9).  
- All ChIP-seq ([GSE203304](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203304)) and RNA-seq ([GSE203305](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203305)) data in this paper are accessible at GEO Accession viewer under the number: [GSE203306](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203306).

## I. Upstream analysis

- `RNA-seq_analysis` : [RNA-seq_analysis.sh](https://github.com/jlchen5/CELR-D-22-00015/blob/main/RNA-seq_analysis.sh)    
RNA-seq data is applied quality and adapter trimming by [Trim Galore v0.6.4_dev](https://zenodo.org/badge/latestdoi/62039322) and mapped to the [mm9 reference genome](http://genome.ucsc.edu/cgi-bin/hgGateway?clade=mammal&org=Mouse&db=mm9) using [STAR v2.7.2d](https://github.com/alexdobin/STAR/releases/tag/2.7.2d). And we performed reads summarization for genomic features via [featureCounts v2.0.0](http://subread.sourceforge.net/featureCounts.html).

- `ChIP-seq_analysis` : [ChIP-seq_analysis.sh](https://github.com/jlchen5/CELR-D-22-00015/blob/main/ChIP-seq_analysis.sh)  
Adapters and low-quality bases were removed using [Trim Galore v0.6.4_dev](https://zenodo.org/badge/latestdoi/62039322). ChIP-seq reads were mapped to the [mm9 reference genome](http://genome.ucsc.edu/cgi-bin/hgGateway?clade=mammal&org=Mouse&db=mm9) using [Bowtie2 v2.3.5.1](https://github.com/BenLangmead/bowtie2/releases/tag/v2.3.5.1). [Samtools v1.10](https://github.com/samtools/samtools/releases/tag/1.10) was used to filter and convert file formats. To remove PCR duplicate reads, which might occur false positive results, [MarkDuplicates (v4.1.4.1, Picard)](https://gatk.broadinstitute.org/hc/en-us/sections/360007458971-4-1-4-1) was used to locate and tag duplicate reads. And sequencing duplicates were removed using the `--REMOVE_SEQUENCING_DUPLICATES` options. ChIP-seq peaks were called using [MACS2 v2.2.6](https://github.com/macs3-project/MACS/releases/tag/v2.2.6).  

- `Repetitive Elements binding analysis`  
For the repetitive Elements binding analysis, the reads were then aligned to the mm9 genome assembly using [STAR v2.7.2d](https://github.com/alexdobin/STAR/releases/tag/2.7.2d) with the options `--alignIntronMax 1` and `--alignEndsType EndToEnd` as previously reported. The parameter `--outFilterMultimapNmax 1` was applied to include only the uniquely mapped reads. Duplicate reads were then removed using  [MarkDuplicates (v4.1.4.1, Picard)](https://gatk.broadinstitute.org/hc/en-us/sections/360007458971-4-1-4-1). Replicate samples were merged using the [Samtools v1.10](https://github.com/samtools/samtools/releases/tag/1.10).

## II. Downstream analysis

- The raw read counts normalization and differential expression were determined by [RNA-seq_rif1.R](https://github.com/jlchen5/CELR-D-22-00015/blob/main/RNA-seq_rif1.R).
 
  > [coding_genes.R](https://github.com/jlchen5/CELR-D-22-00015/blob/main/coding_genes.R) is for further analysis of genes  
  > [repeats.R](https://github.com/jlchen5/CELR-D-22-00015/blob/main/repeats.R) is for further analysis of repeats

- Correlation among different samples is analyzed via [CalculateCorrelation.R](https://github.com/jlchen5/CELR-D-22-00015/blob/main/CalculateCorrelation.R) and [plotCorrHeatmap.R](https://github.com/jlchen5/CELR-D-22-00015/blob/main/plotCorrHeatmap.R) is used to plot correlation heatmap.
- Principal Component Analysis (PCA) is calculated the foldchange (vs. control/wild type) between the samples with the depletion of the indicated genes to control (or wild type) samples via [PCAanalysis.R](https://github.com/jlchen5/CELR-D-22-00015/blob/main/PCAanalysis.R).
- [GOanalysis.R](https://github.com/jlchen5/CELR-D-22-00015/blob/main/GOanalysis.R) is used to perform GO analysis.
- Single locus repeats RNA-seq analysis needs single locus position annotation file to annotate the reads signal in genome. 
  ~~~
  featureCounts sample.repeats.Aligned.sortedByCoord.out.bam \
                -a ~/project/rif1/mm9_repeats_single_locus.gtf \
                -g gene_id2 -o sample.counts -T 40
  ~~~


## III. Quality Control

All data in this article has been quality controlled, the qc reports are as follows:   
 - [ChIP-seq_multiqc_report](https://htmlpreview.github.io/?https://github.com/jlchen5/CELR-D-22-00015/blob/main/QC/ChIP-seq_multiqc_report.html)  
 - [RNA-seq_multiqc_report](https://htmlpreview.github.io/?https://github.com/jlchen5/CELR-D-22-00015/blob/main/QC/RNA-seq_multiqc_report.html)  
 
