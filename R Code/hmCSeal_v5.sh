#!/bin/bash

export NameSeal=`cat full_dataset_sample_names.txt`
#export NameSeal=`cat test_dataset_sample_names.txt`
export Ref="/home/turagapm/reference_genomes"
export Project="/home/turagapm/peritoneal_processing"
LOG="$Project/logs/hmCSeal_bowtie2_runv4_091823.log"

# for i in $NameSeal; do
# TRIM_REPORT="$Project/trimmed_data_bam_9_18_2023/${i}_R1.fastq.gz_trimming_report.txt"
# # Skip if trimming report already exists
# if [ ! -e $TRIM_REPORT ]; then
# echo $i
# echo "============================================= Start mapping ${i} by bowtie2 ============================================" >> $LOG
# trim_galore -o $Project/trimmed_data_bam_9_18_2023 --paired --length 15 --gzip $Project/raw_data_fastq/${i}_R1.fastq.gz $Project/raw_data_fastq/${i}_R2.fastq.gz --cores 4
# bowtie2 -q -N 1 -p 48 -x $Ref/bowtie2/hg19 -1 $Project/trimmed_data_bam_9_18_2023/${i}_R1_val_1.fq.gz -2 $Project/trimmed_data_bam_9_18_2023/${i}_R2_val_2.fq.gz 2>> $LOG | \
# grep -v "XS:i:" | \
# samtools sort -@ 16 -m 8G - | \
# samtools rmdup - $Project/trimmed_data_bam_9_18_2023/${i}_bowtie2.bam 2>> $LOG
# samtools index $Project/trimmed_data_bam_9_18_2023/${i}_bowtie2.bam
# fi
# done

OUTPUT_COUNT1="$Project/DESeq2/ReadsCount_AutosomeGenes_ProteinCoding_v4_091823.txt"
OUTPUT_COUNT2="$Project/DESeq2/ReadsCount_AutosomeGenes_promoter_and_enhancer_v4_091823.txt"
OUTPUT_COUNT3="$Project/DESeq2/ReadsCount_AutosomeGenes_macs2_pe5_040124.txt"

# Check and run featureCounts if output does not exist

if [ ! -e $OUTPUT_COUNT1 ]; then
 featureCounts -T 64 -a $Ref/hg19_genebody_fromZhangLab.saf -F SAF -p -Q 10 -C --ignoreDup -o $OUTPUT_COUNT1 $Project/trimmed_data_bam_9_18_2023/*_bowtie2.bam 2>>$LOG
fi

if [ ! -e $OUTPUT_COUNT2 ]; then
 featureCounts -T 64 -a $Ref/hg19_genehancer.saf -F SAF -p -Q 10 -C --ignoreDup -o $OUTPUT_COUNT2 $Project/trimmed_data_bam_9_18_2023/*_bowtie2.bam 2>>$LOG
fi

if [ ! -e $OUTPUT_COUNT3 ]; then
  featureCounts -T 64 -a $Ref/merged_peaks_filtered_pe5_040124.saf -F SAF -p -Q 10 -C --ignoreDup -o $OUTPUT_COUNT3 $Project/trimmed_data_bam_9_18_2023/*_bowtie2.bam 2>>$LOG
fi