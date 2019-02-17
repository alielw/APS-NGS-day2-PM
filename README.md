**_Advanced Data Analysis - Introduction to NGS data analysis_**<br>
*Department of Animal and Plant Sciences, University of Sheffield*

# Differential gene expression analyses
#### Alison Wright

The aim of this practical is to learn how to perform differential gene expression analyses. We will be using a dataset of expression data for 4 individuals of *Heliconius melpomene*. For each individual, two different wing regions have been sequenced. We will try to identify genes that are differentially expressed between wing regions. Samples are labelled I or A. I is the part of the wing that is iridescent, A is the top part of the wing, which is called the androchonial region.

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
        /fastdata/$USER/DE

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

* Check that the BAM files are sorted, and if so whether they are sorted by read name or by alignment position. You only need to check one file (96I).

        samtools view -H 96I.bam | head

* Htseq-count takes a couple of hours to run, so lets just submit one job to ShARC as practice (96I). You can run this in interactive mode. Remember to specify with -r whether the bam file is sorted by read name (`name`) or alignnment position (`pos`)

        qrsh
        
        source /usr/local/extras/Genomics/.bashrc
        
        cd /fastdata/$USER/DE/Tophat_output

        htseq-count -f BAM 96I.bam /fastdata/$USER/align/Cufflinks_output -r [name/pos] > 96I.htseq 

* We have already generated read count files for all the samples. There should be 8 files, one for each sample.

        /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/NGS_data/htseqCounts
        
* Copy the read count folder into your fastdata folder.

        cp -r /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/NGS_data/htseqCounts /fastdata/$USER/DE/
        
* We need to merge the read counts across samples so we have one file with all the read counts. You can do this using a custom python script 

        python script.py merged_readcounts.txt
        
* The downstream analyses we will perform next in R on the desktop. Copy the merged read count file onto your desktop. You need to open a new terminal window to do this.

        cd Desktop

        scp $USER@sharc.shef.ac.uk:/usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/NGS_data/htseqCounts/merged_readcounts.txt .

---

## 2.Introduction to EdgeR

We will use EdgeR for the next analyses ... in R

        data <- read.table("~/Desktop/merged_readcounts.txt",stringsAsFactors=F,header=T, row.names=1)

        names(data)

        dim(data)
        
Need to specify conditions. Include a table of samples names and id.

        conditions <- factor(c("APLFG","APLFG","APLFG","APLFG","APLFG","APLMG","APLMG","APLMG"))

Object

        expr <- DGEList(counts=data,group=conditions)
        
        $samples

---

## 3.Filter expression data

Explain cpm

        cpm_data <- cpm(expr)

Filter expression, at least 1 cpm in either sample. Recommended by EdgeR. Dont exclude sample-limited genes.

        keep <- rowSums(cpm(expr)>1) >=1
        expr_filtered <- expr[keep,,keep.lib.sizes=FALSE]
        dim(expr_filtered)

## 4. Normalisation of gene expression

Explain how cpm accounts for differences in library size.

RPKM accounts for differences in gene length

TMM for cross species differences in GC content etc

expr_norm = normalise(expr, conditions, species_name)

---

## 5. Identify differentially expressed genes

#Estimating dispersions
expr_norm <- estimateCommonDisp(expr_norm)
expr_norm <- estimateTagwiseDisp(expr_norm)

levels(expr_norm $samples$group)

et <- exactTest(object, pair=c(A,B))
print(et$comparison)

#Correction for multiple testing
p <- et$table$PValue
p_FDR <- p.adjust(p, method = c("fdr"), n = length(p))
print(length(p_FDR))

table_de <- et$table
table_de$Padj <- p_FDR

print("p<0.05")
print(length(which(table_de$Padj < 0.05)))
print("Male-biased-pvalue&fc")
print(length(which(table_de$Padj < 0.05 & table_de$logFC > 1)))
print("Male-biased-fc")
print(length(which(table_de$logFC > 1)))
print("Female-biased-pvalue&fc")
print(length(which(table_de$Padj < 0.05 & table_de$logFC < -1)))
print("Female-biased-fc")
print(length(which(table_de$logFC < -1)))
	
filename = paste("~/Desktop/Sex-bias-fdr",A,B,name, sep="-")
print(filename)
write.table(table_de, file=filename,quote=F, sep="\t")
}

---

## 4. Visualisation of gene expression

---
