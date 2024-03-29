**_Advanced Data Analysis - Introduction to NGS data analysis_**<br>
*Department of Animal and Plant Sciences, University of Sheffield*

# Differential gene expression analyses
#### Alison Wright, Nicola Nadeau

The aim of this practical is to learn how to perform differential gene expression analyses. We will be using a dataset of expression data for 4 individuals of *Heliconius melpomene*. For each individual, two different wing regions have been sequenced. We will try to identify genes that are differentially expressed between wing regions. Samples are labelled I or A. I is the part of the wing that is iridescent, A is the top part of the wing, which is called the androchonial region.

## Table of contents 

1. Obtaining expression values

2. Introduction to edgeR

3. Filter expression data

4. Normalisation of gene expression

5. Visualisation of gene expression

6. Identify differentially expressed genes

---

## PRACTICAL - Initial set up
First of all, this tutorial must be run using an interactive session in ShARC. You will also submit jobs to ShARC. For that, you should log in into ShARC with `ssh`, and then request an interactive session with `qrsh`. Your shell prompt should show `sharc-nodeXXX` (XXX being a number between 001 and 172) and not `@sharc-login1` nor `@sharc-login2`.

	qrsh

For this particular tutorial, we are going to create and work on a directory called `DE` in your /fastdata/$USER directory:

        cd /fastdata/$USER
        mkdir 3.DE
        cd 3.DE

Remember you can check your current working directory anytime with the command `pwd`.
It should show something like:

        pwd
        /fastdata/$USER/3.DE

---

## READ - 1. Obtaining expression values

Last session, you ran StringTie on one sample to assemble transcripts. This generated two files `gene_abund` and `gtf`. A `gene_abund` file looks like this:

![alt text](https://github.com/alielw/APS-NGS-day2-PM/blob/master/StringTie.png)

As you can see, StringTie also calculates coverage or expression of these transcripts. It calculates FPKM, which is Fragments Per Kilobase of exon model per Million mapped fragments. It is possible to use this FPKM information for an initial assessment of gene expression but expression analyses may require further processing of this data (eg normalisation). Therefore, we do not recommend conducting analyses on these raw FPKM values. Instead, we need to extract read counts for each gene for downstream analyses.

Link to [BLOG](https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/) explaining difference between FPKM and TPM.

A `gtf` file slooks like this and contains information on the physical location of each transcript.

![alt_test](https://github.com/alielw/APS-NGS-day2-PM/blob/master/gtf.png)

We have run StringTie on all samples for you. Specifically, we have StringTie output files for 4 individuals of *Heliconius melpomene*. For each individual, two different wing regions have been sequenced. We will try to identify genes that are differentially expressed between wing regions. Samples are labelled I or A. I is the part of the wing that is iridescent, A is the top part of the wing, which is called the androchonial region.

![alt text](https://github.com/alielw/APS-NGS-day2-PM/blob/master/Samples_id.png)

## PRACTICAL - 1. Obtaining expression values

* First let's download the StringTie output files.

		cp -r /usr/local/extras/Genomics/workshops/NGS_AdvSta_2020/NGS_data/StringTie_output /fastdata/$USER/3.DE
		
* Let's check what they look like

		cd StringTie_output

You can check a `gene_abund` file using the command below. The format should be similar to the image above.

	head 61A.trim_StringTie.gene_abund

You can check a `gtf` file using the command below. The format should be similar to the image above.

	head 61A.trim_StringTie.gtf

* How many StringTie output files are there? Does it correspond to the number of samples?

		ls *.gene_abund | wc -l

* How many transcripts are in each `gene_abund` file? Remember you can use the `wc -l` command from the last session to count the number of lines in a file.

			wc -l 61A.trim_StringTie.gene_abund

* Now, let's extract the names and paths of the `gtf` StringTie output file.
	
		for i in *.gtf; do echo ${i:0:3} `pwd`/$i; done

* Now, let's save these names and paths into a new file called `StringTie_filenames.txt`

		for i in *.gtf; do echo ${i:0:3} `pwd`/$i; done > StringTie_filenames.txt

* Check how many lines are in `StringTie_filenames.txt` using `wc -l`. Is it what you expect?

		wc -l StringTie_filenames.txt
		
* Next, we will extract raw read counts. StringTie provides a Python script `prep_DE.py` to extract this read count information directly from the files generated by StringTie (run with the -e parameter). Let's download this python script.

		cp /usr/local/extras/Genomics/workshops/NGS_AdvSta_2020/NGS_data/prep_DE.py /fastdata/$USER/3.DE
		
* Let's now run this python script to extract raw read counts. We also need to load python as well. This will generate generates two CSV files containing the count matrices for genes and transcripts. The parameter `-i` is specifying a text file containing a list of paths to GTF files, `-g` is where to output the gene count matrix and `-t` is where to output the transcript count matrix.

		cd /fastdata/$USER/3.DE
		
		module load apps/python/anaconda2-4.2.0

		python prep_DE.py -i /fastdata/$USER/3.DE/StringTie_output/StringTie_filenames.txt -g /fastdata/$USER/3.DE/StringTie_output/gene_count_matrix.csv -t /fastdata/$USER/3.DE/StringTie_output/transcript_count_matrix.csv -s HMEL

* Let's look at a bit of the `gene_count_matrix.csv` file. 

		cd /fastdata/$USER/3.DE/StringTie_output/
		head -100 gene_count_matrix.csv

* How many columns are there? Is this what you expect? What is each column? Can you identify any genes by eye that look to be differentially expressed between I and A samples?

* Download the `gene_count_matrix.csv` file to your Desktop.

If you are using mac, you need to open a new terminal window. If you use the following commands in this new window, your file will download onto your desktop.

	cd Desktop
	scp $USER@sharc.sheffield.ac.uk:/fastdata/$USER/3.DE/StringTie_output/gene_count_matrix.csv .


---

## PRACTICAL - 2. Introduction to edgeR

We will use [edgeR](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) to perform differential gene expression analyses. This is implemented in R on your Desktop.

* Load R.

* Install the [Bioconductor](https://www.bioconductor.org) package in R.

* Load the edgeR library

        library(edgeR)

* Read in the data to R.

        data <- read.csv("~/Desktop/gene_count_matrix.csv",stringsAsFactors=F,header=T, row.names=1)

        names(data)

        dim(data)
	
		head(data)

* We have two treatments or conditions: A and I. Find which samples correspond to which treatment from the file header.

* Specify the order of the treatments in the file by filling in the `""` below (e.g "A" or "I").

        treatment <- factor(c("","","","","","","",""))

* Check you specified this correctly? 

		treatment

* [edgeR](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) stores data in a simple list-based data object called a `DGEList`. This type of object is easy to use because it can be manipulated like any list in R. The function readDGE makes a DGEList object directly. If the table of counts is already available as a matrix or a data.frame, x say, then a DGEList object can be made by:

        expr <- DGEList(counts=data)

* We can add samples at the same time:

        expr <- DGEList(counts=data, group=treatment)

        expr$samples

---

## PRACTICAL - 3. Filter expression data

Genes with very low counts across all samples provide little evidence for differential expression. From a biological point of view, a gene must be expressed at some minimal level before it is likely to be translated into a protein or to be biologically important. Therefore, we need to filter genes with biologically irrelevant expression. 

The developers of edgeR recommend that gene is required to have a count of 5-10 in a library to be considered expressed in that library. However, users should filter with count-per-million (`CPM`) rather than filtering on the read counts directly, as the latter does not account for differences in library sizes between samples. Therefore, they recommend filtering on a CPM of 1.

* We can calculate count-per-million (`CPM`) using cpm(`DGEList`).

        cpm_data <- cpm(expr)
		head(cpm_data)

* Filter expression to remove lowly expressed genes. We do not want to exclude genes specific to one wing region (ie I or A-limited genes). Therefore, a sensible filtering approach is to filter genes that are not expressed > 1 CPM in at least half of the samples. We have 8 samples.

        keep <- rowSums(cpm(expr)>1) >=4
        expr_filtered <- expr[keep,,keep.lib.sizes=FALSE]

* How many genes remain? Does this seem a sensible number of expressed genes?
	
		dim(expr)
        	dim(expr_filtered)

* Try some different filtering thresholds (e.g. > 2 CPM). Does it have a big effect on the number of genes that are filtered?

## READ - 4. Normalisation of gene expression

We need to normalise gene expression across samples before conducting differential gene expression analyses. There are a number of sources of variation we need to account for.

* **Sequencing depth**
The most obvious technical factor that affects the read counts, other than gene expression levels, is the sequencing depth of each sample. edgeR controls for varying sequencing depths as represented by differing library sizes. This is part of the basic modeling procedure and flows automatically into downstream statistical analyses. It doesn’t require any user intervention.

* **RNA composition**
The second most important technical influence on differential expression is one that is RNA composition. RNA-seq provides a measure of the relative abundance of each gene in each RNA sample, but does not provide any measure of the total RNA output on a per-cell basis. This commonly becomes important when a small number of genes are very highly expressed in one sample, but not in another. The highly expressed genes can consume a substantial proportion of the total library size, causing the remaining genes to be under-sampled in that sample. Unless this RNA composition effect is adjusted for, the remaining genes may falsely appear to be down-regulated in that sample. This is normally a problem for cross-species comparisons. 

## PRACTICAL - 4. Normalisation of gene expression

* Variation in RNA composition can be accounted for by the calcNormFactors function. This normalizes for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes. The default method for computing these scale factors uses a trimmed mean of M- values (TMM) between each pair of samples. Further information can be found in the [edgeR manual](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) 

		expr_norm = calcNormFactors(expr_filtered)
		expr_norm

* Let's check that the normalisation has worked correctly. If so, each sample should have similar read counts.

		cpm_log <- cpm(expr_norm, log = TRUE)
		
		plot(density((cpm_log[,1])), col="red")
		lines(density((cpm_log[,2])), col="red")
		lines(density((cpm_log[,3])), col="red")
		lines(density((cpm_log[,4])), col="red")
		lines(density((cpm_log[,5])), col="red")
		lines(density((cpm_log[,6])), col="blue")
		lines(density((cpm_log[,7])), col="blue")
		lines(density((cpm_log[,8])), col="blue")

* Do the graphs look as you would expect? Are there similar read counts across samples?

---

## PRACTICAL - 5. Visualisation of gene expression

There are a number of ways to visualise gene expression data. Heatmaps and PCA plots are widely used approaches. We will cover heatmaps today.

* Calculate log CPM

        cpm_log <- cpm(expr_norm, log = TRUE)

* Install pvclust from the `Package Installer` in R.
	
* Perform clustering of expression
	
        library(pvclust)
	
        bootstraps = pvclust(cpm_log, method.hclust="average", method.dist="euclidean")

        plot(bootstraps)

How do the samples cluster? What can you conclude about the genomic architecture of wing iridescence? Is it likely that many genes with small expression differences or a few genes with big expression changes encode this phenotype?
	
* Plot heatmap
	
        library(pheatmap)

        pheatmap(cpm_log, show_colnames=T, show_rownames=F, clustering_distance_cols = "euclidean", clustering_method="average") 

---

## PRACTICAL - 6. Identify differentially expressed genes in a pairwise comparison

Next we can use edgeR to identify differentially expressed genes between the two groups of samples. We have a pairwise comparison between two groups. Further information about each step can be found in the [edgeR manual](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf). You can also find other approaches that might be more appropriate if you have more than two treatments. 

Briefly, edgeR uses the quantile-adjusted conditional maximum likelihood (qCML) method for experiments with single factor. The qCML method calculates the likelihood by conditioning on the total counts for each tag, and uses pseudo counts after adjusting for library sizes. The edgeR functions `estimateCommonDisp` and `estimateTagwiseDisp` produce a pseudo counts. We can then proceed with determining differential expression using the `exactTest` function. The exact test is based on the qCML methods. We can compute exact p-values by summing over all sums of counts that have a probability less than the probability under the null hypothesis of the observed sum of counts. The exact test for the negative binomial distribution has strong parallels with Fisher’s exact test.

* Estimating common dispersion.

        expr_filtered_norm<- estimateCommonDisp(expr_norm)

* Estimating tagwise dispersion.

        expr_filtered_norm<- estimateTagwiseDisp(expr_filtered_norm)

        levels(expr_filtered_norm$samples$group)

* Exact test to calculate logFC, logCPM and PValue for every gene.

        et <- exactTest(expr_filtered_norm, pair=c("A","I"))

        et

* Correct for multiple testing.

        p <- et$table$PValue

        p_FDR <- p.adjust(p, method = c("fdr"), n = length(p))

        print(length(p_FDR))

        table_de <- et$table

        table_de$Padj <- p_FDR

* Print results to a new file.

        filename = paste("~/Desktop/Wing-bias-fdr.txt")

        print(filename)

        write.table(table_de, file=filename,quote=F, sep="\t")
	
* Identify differentially expressed genes between I and A.

Normally, we use a fdr p-value threshold < 0.05. It is also important to consider imposing a fold change threshold eg logFC of 1. As we conducted `exactTest(expr, pair=c("A","I"))`, positive logFC means I > A (I-biased), negative logFC means A > I (A-biased).

How many genes are expressed more in the iridescent region of the wing?

		print(length(which(table_de$Padj < 0.05 & table_de$logFC > 1)))
	
		table_de[which(table_de$Padj < 0.05 & table_de$logFC > 1),]

How many genes are expressed more in the androchonial region of the wing?
	
		print(length(which(table_de$Padj < 0.05 & table_de$logFC < -1)))
	
		table_de[which(table_de$Padj < 0.05 & table_de$logFC < -1),]
	
* Lets look up and check whether any of these differentially expressed genes are annotated, and find out gene features. A note of caution here, the *Heliconius melpomene* is not well annotated relative to genomes of model species.

	Look up the gene id in [Lepbase](http://ensembl.lepbase.org/Heliconius_melpomene_melpomene_hmel2/Info/Index). You need to drop the transcript info from the gene name.
	
	What information can you find out about the gene? Does it have any orthologs and if so can you infer the function of this gene? You can use [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to help with this.
	
	What about other genes that you identified as differentially expressed? What can you find out about these? Are there any good candidates for iridescence?

* Let's try to relax the criteria for differential expression and remove the fold change threshold. (HINT: you can tweak the above commands to do this).

	How many genes are expressed more in the iridescent region of the wing?

	Do any of these genes look like good candidates for iridescence?
	
---

Return to the main course page: https://github.com/njnadeau/AdvDataAna-introNGS
