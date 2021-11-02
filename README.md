# 5hmc_Pathway_Analysis
## Bash scripts:
### 1- Fastq file processing
Location: ~/fastq_processing/Scripts/hmCSeal_v4.sh
   This script takes as input .fastq files and a list of sample names, uses trimgalor & bowtie 2 to trim & align reads to a refence genome. It outputs a readcounts file for use in DESeq2 and other sequencing QC data which is output to a log file.

### 2- Bowtie2 reference genome creation
Location: ~/fastq_processing/Reference_genomes/bowtie2/make_hg19.sh
   This script takes as input the hg19 reference genome, and creates the necessary reference genome files for use by bowtie2. This is a pre-requisite for script 1.
   
## R scripts:
### 3- RawCounts file creation
Location: ~/R Code/rawcounts_processing.R
This script takes as input raw counts file(s) outputted from script 1, as well as a sample conditions or phenotype file, and outputs a formatted ReadCounts file for use in DESeq2 analysis in script 5 as well as use in downstream tools such as GSEA and ssGSEA. It has several options to allow the combining of multiple ReadCounts files, as well as several normalization options.

### 4- Group randomizer 
Location: ~/R Code/dataset_group_randomizer.R
This script is optional. It allows for the splitting up of a formatted ReadCounts file from script 3 into two groups randomly with a pre-set ratio (e.g. 70%-group 1; 30%-group 2.)

### 5- Differential Gene Analysis using DESeq2 
Location: ~/R Code/DESeq2_analysis.R
This script takes as input the formatted files from script 3 and uses DESeq2 to output standard analysis charts such as: Heatmaps, PCA, UMAP, Volcano and result lists.

### 6- Gene list similarity comparison
Location: ~/R Code/GeneList_Comparison.R
This script takes as input two rank-ordered gene lists (ranked by log2foldchange) taken as output from script 5, and provides several statistics on how similar these ranked lists are.

### 7- GSEA results heatmap chart
Location: ~/R Code/GeneSet_heatmap.R
This script takes as input a list of gene sets, as a result of GSEA analysis, and provides a heatmap of the results with each row representing a specifi gene set.

### 8 - Elastic net regularization to the multivariate logistic regression models using glmnet
Location: ~/R Code/glmnet_prediction.R
This script takes as input the training and validation sets output by script 4, and uses the training group to generate a logistic regression model which is then tested against the validation set.

## References:
https://github.com/hbctraining/DGE_workshop/tree/master/lessons
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8
https://intobioinformatics.wordpress.com/2019/06/08/running-umap-for-data-visualisation-in-r/
https://bioconductor.org/packages/release/bioc/manuals/OrderedList/man/OrderedList.pdf
https://mirrors.nju.edu.cn/bioconductor/bioc/1.9/vignettes/OrderedList/inst/doc/tr_2006_01.pdf

## Copyright
These are open access materials distributed under the terms of the Creative Commons Attribution license (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.

## Author:
Yuri Malina
ymalina@gmail.com
