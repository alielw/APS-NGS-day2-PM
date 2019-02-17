**_Advanced Data Analysis - Introduction to NGS data analysis_**<br>
*Department of Animal and Plant Sciences, University of Sheffield*

# Differential gene expression analyses
#### Alison Wright

The aim of this practical is to learn how to perform differential gene expression analyses. We will be using a dataset of expression data for 4 individuals of *Heliconius melpomene*. For each individual, two different wing regions have been sequenced. We will try to identify genes that are differentially expressed between wing regions.

## Table of contents 

1. Obtaining expression values

2. Normalisation of gene expression

3. Identify differentially expressed genes

4. Visualisation of gene expression

---

## Initial set up
First of all, this tutorial must be run using an interactive session in ShARC. You will also submit jobs to ShARC. For that, you should log in into ShARC with `ssh`, and then request an interactive session with `qrsh`. Your shell prompt should show `sharc-nodeXXX` (XXX being a number between 001 and 172) and not `@sharc-login1` nor `@sharc-login2`.

For this particular tutorial, we are going to create and work on a directory called `DE` in your /fastdata/$USER directory:

        cd /fastdata/$USER
        mkdir DE
        cd DE

Remember you can check your current working directory anytime with the command `pwd`.
It should show something like:

        pwd
        /fastdata/myuser/DE

---

## 1.Obtaining expression values

Discussion of different measures of gene expression: RPKM,FPKM,CPM. Need read counts to calculate these.

**Cufflinks output files: isoforms.fpkm_tracking & genes.fpkm_tracking**
Cufflinks generates files that contain the estimated isoform or gene-level expression values as Fragments Per Kilobase of exon model per Million mapped fragments (FPKM). It is possible to use this FPKM information for an initial assessment of gene expression but expression analyses may require further processing of this data (eg normalisation). Therefore, we do not recommend conducting analyses on these raw FPKM values.

Instead, we can use [HTSeq](https://htseq.readthedocs.io/en/release_0.11.1/) to extract read counts for downstream analyses.

* **Sort BAM files with [Samtools](http://www.htslib.org/doc/samtools-1.0.html).**

HTSeq requires BAM files to be sorted either by read name or by alignment position. 

        >samtools sort file.bam - o file.sorted.bam

The default sort is by coordinate. However, you can use the sort order (SO) flag in the BAM header to check if the file has been sorted. 

        >samtools view -H file.bam | head

* **Extract read counts.**

Use [htSeq-count](https://htseq.readthedocs.io/en/release_0.11.1/count.html) to extract read counts for each transcript.

        >htseq-count -f BAM/SAM file.bam merged.gtf -r order > file.htseq

**-r** For paired-end data, the alignment have to be sorted either by read name or by alignment position. If your data is not sorted, use the samtools sort function of samtools to sort it. Use this option, with name or pos for <order> to indicate how the input data has been sorted. The default is name.
    
**-s <yes/no/reverse>** Specify whether the data is from a strand-specific assay (default: yes)
    
htseq-count outputs a table with counts for each feature, followed by the special counters, which count reads that were not counted for any feature for various reasons. 

## a. PRACTICAL ACTIVITY

Extract read counts with htseq.

* First, copy the folder containing the BAM files for our 8 samples (4 individuals, 2 samples per individual) to your fastdata

        cp -r /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/NGS_data/Tophat_output /fastdata/$USER/DE

* Check that the BAM files are sorted, and if so are they sorted by read name or by alignment position. You only need to check a handful.

        samtools view -H file.bam | head

* Htseq-count takes a couple of hours to run, so lets just submit one job to ShARC as practice (96I). You can run this in interactive mode. Remember to specify with -r whether the bam file is sorted by read name (name) or alignnment position (pos)

        qrsh
        
        source /usr/local/extras/Genomics/.bashrc
        
        cd /fastdata/$USER/DE/Tophat_output

        htseq-count -f BAM 96I.bam /fastdata/$USER/align/Cufflinks_output -r [name/pos] > 96I.htseq 

* We have already generated read count files for all the samples for you. There should be 8 files, one for each sample.

        /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/NGS_data/htseqCounts
        
* Copy the read count folder into your fastdata folder.

        cp -r /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/NGS_data/htseqCounts /fastdata/$USER/DE/
        
* The downstream analyses we will perform next in R on the desktop. Copy the read count folder onto your desktop.

        cd Desktop

        scp -r /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/NGS_data/htseqCounts .

## 2. Normalisation of gene expression



## 3. Identify differentially expressed genes

## 4. Visualisation of gene expression
