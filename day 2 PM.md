# Differential gene expression analyses

Pipeline to 

1. Obtaining expression values

2. Normalisation of gene expression

3. Identify differentially expressed genes

4. Visualisation of gene expression

## 1.Obtaining expression values

Discussion of different measures of gene expression: RPKM,FPKM,CPM. Need read counts to calculate these.

**Cufflinks output files: isoforms.fpkm_tracking & genes.fpkm_tracking**
Thes files contains the estimated isoform or gene-level expression values as Fragments Per Kilobase of exon model per Million mapped fragments (FPKM). It is possible to use this FPKM information for an initial assessment of gene expression but expression analyses may require further processing of this data (eg normalisation). Therefore, we do not recommend conducting analyses on these raw FPKM values.

Instead, we can use [HTSeq](https://htseq.readthedocs.io/en/release_0.11.1/) to extract read counts for downstream analyses.

* **Sort BAM files with [Samtools](http://www.htslib.org/doc/samtools-1.0.html).**

HTSeq requires BAM files to be sorted either by read name or by alignment position. 

        >samtools sort file.bam - o file.sorted.bam

The default sort is by coordinate. However, you can use the sort order (SO) flag in the BAM header to check if the file has been sorted. 

        >samtools view -H file.bam | head

* **Extract read counts.**

Use [htSeq-count](https://htseq.readthedocs.io/en/release_0.11.1/count.html) to extract read counts for each transcript.

        >htseq-count -f BAM/SAM file.bam merged.gtf -r order > file.htseq

* **Sort BAM files with [Samtools](http://www.htslib.org/doc/samtools-1.0.html).**

HTSeq requires BAM files to be sorted either by read name or by alignment position. 

        >samtools sort file.bam - o file.sorted.bam
    
**-r** For paired-end data, the alignment have to be sorted either by read name or by alignment position. If your data is not sorted, use the samtools sort function of samtools to sort it. Use this option, with name or pos for <order> to indicate how the input data has been sorted. The default is name.
    
**-s <yes/no/reverse>** Specify whether the data is from a strand-specific assay (default: yes)
    
htseq-count outputs a table with counts for each feature, followed by the special counters, which count reads that were not counted for any feature for various reasons. 

## 2. Normalisation of gene expression

## 3. Identify differentially expressed genes

## 4. Visualisation of gene expression
