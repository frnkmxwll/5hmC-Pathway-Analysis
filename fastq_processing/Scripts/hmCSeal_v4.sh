#!/bin/bash

export NameSeal=`cat full_dataset_sample_names.txt`
#export NameSeal=`cat test_dataset_sample_names.txt`
export Ref="../reference_genomes"
export Project="../../fastq_processing"

LOG="$Project/logs/hmCSeal_bowtie2_run4.log"

# samtools sort -@ 16 -m 8G --input-fmt-option="filter=[XM]<=1" - | \
for i in $NameSeal;do
	echo $i
	echo "============================================= Start mapping ${i} by bowtie2 ============================================" >> $LOG
	trim_galore -o $Project/trimmed_data_bam --paired --length 15 --gzip $Project/raw_data_fastq/${i}_R1.fastq.gz $Project/raw_data_fastq/${i}_R2.fastq.gz
	bowtie2 -q -N 1 -p 48 -x $Ref/bowtie2/hg19 -1 $Project/trimmed_data_bam/${i}_R1_val_1.fq.gz -2 $Project/trimmed_data_bam/${i}_R2_val_2.fq.gz  2>> $LOG | \
	grep -v "XS:i:" | \
	samtools sort -@ 16 -m 8G - | \
	samtools rmdup - $Project/trimmed_data_bam/${i}_bowtie2.bam 2>> $LOG
	samtools index $Project/trimmed_data_bam/${i}_bowtie2.bam
done

featureCounts -a $Ref/hg19_genebody_fromZhangLab.saf -F SAF -p -Q 10 -C --ignoreDup -o $Project/DESeq2/ReadsCount_AutosomeGenes_ProteinCoding_v4run3.txt $Project/trimmed_data_bam/*_bowtie2.bam 2>>$Project/DESeq2/ReadsCount_AutosomeGenes_ProteinCoding_bowtie2_run3.log

featureCounts -a $Ref/hg19_genehancer.saf -F SAF -p -Q 10 -C --ignoreDup -o $Project/DESeq2/ReadsCount_AutosomeGenes_promoter_and_enhancer_v4run3.txt $Project/trimmed_data_bam/*_bowtie2.bam 2>>$Project/DESeq2/ReadsCount_AutosomeGenes_promoter_and_enhancer_bowtie2_run3.log
