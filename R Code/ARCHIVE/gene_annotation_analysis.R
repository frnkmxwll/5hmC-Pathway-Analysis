# List of packages
packages <- c("csaw", "GenomicRanges", "TxDb.Hsapiens.UCSC.hg19.knownGene", "ggplot2", "reshape2")

# Function to install missing packages
install_missing <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
}

# Use BiocManager for Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install missing CRAN packages
#install_missing(c("ggplot2", "reshape2"))

# Install missing Bioconductor packages
#BiocManager::install(c("csaw", "GenomicRanges", "TxDb.Hsapiens.UCSC.hg19.knownGene"))

library(csaw)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(reshape2)
library(Rsamtools)
library(rtracklayer)
library(GenomicAlignments)
library(CGPfunctions) #needed for plotxtabs2 

# 1. Read in BAM files and the conditions CSV

bam_dir <- "~/peritoneal_processing/trimmed_data_bam_9_18_2023"
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)

# Index each BAM file
lapply(bam_files, function(bam_file) {
  if (!file.exists(paste0(bam_file, ".bai"))) { # Check if index exists
    indexBam(bam_file)
  }
})

conditions_df <- read.csv("~/5hmC-Pathway-Analysis/Raw Input/Working Inputs/1-2-3-4vs8-9-11_fn.csv")
# define contrast groups
groups <- unique(meta[c("condition")])
contrast_groups <- c("condition",groups[1,1], groups[2,1])

# Only keep BAM files that are present in the conditions file
bam_files <- bam_files[basename(bam_files) %in% conditions_df$filename]

bam_conditions <- conditions_df$condition[match(basename(bam_files), conditions_df$filename)]
paste(bam_conditions)

# 2. Identify 5hmC peaks for each BAM file
# Directory to store MACS2 output
macs2_output_dir <- "~/5hmC-Pathway-Analysis/Output/macs2_output/"
dir.create(macs2_output_dir, showWarnings = FALSE)

# Modify the call_peaks_macs2 function to call peaks for individual bam files
call_peaks_macs2 <- function(bam_file, output_dir) {
  base_name <- basename(bam_file)
  output_name <- tools::file_path_sans_ext(base_name)
  
  # Construct the expected output file path
  expected_output_file <- file.path(output_dir, paste0(output_name, "_peaks.narrowPeak"))
  
  # Only call MACS2 if the expected output file doesn't exist
  if (!file.exists(expected_output_file)) {
    cmd <- sprintf("macs2 callpeak -t %s -f BAM -g hs -n %s --outdir %s", 
                   bam_file, 
                   output_name,
                   output_dir)
    system(cmd)
  } else {
    message(paste0("Output file for ", base_name, " already exists. Skipping MACS2 call."))
  }
  
  return(expected_output_file)
}

# Call peaks for individual bam files
peak_files <- lapply(bam_files, function(bam_file) {
  call_peaks_macs2(bam_file, macs2_output_dir)
})

# 3. Merge, filter, and intersect peaks
# Assuming you have bedtools installed on your system

# Merge peak files from all samples
merged_peak_file <- file.path(macs2_output_dir, "merged_peaks.bed")
cmd <- sprintf("cat %s | sort -k1,1 -k2,2n | bedtools merge > %s", 
               paste(peak_files, collapse=" "),
               merged_peak_file)
system(cmd)

# Filter merged peaks: retain peaks in >10 samples and <1000 bp
filtered_peak_file <- file.path(macs2_output_dir, "filtered_peaks.bed")
cmd <- sprintf("awk '($3-$2)<1000' %s | bedtools intersect -wa -c -a - -b %s | awk '$4>1' > %s", 
               merged_peak_file,
               paste(peak_files, collapse=" "),
               filtered_peak_file)
system(cmd)

# Filter out peaks from blacklisted genomic regions (assuming you have a file named "blacklist.bed" in the current directory)
filtered_blacklist_removed_peak_file <- file.path(macs2_output_dir, "filtered_blacklist_removed_peaks.bed")
cmd <- sprintf("bedtools subtract -a %s -b ~/data/peritoneal_processing/macs2_output/hg19-blacklist.v2.bed > %s", 
               filtered_peak_file, 
               filtered_blacklist_removed_peak_file)
system(cmd)

# Exclude peaks in chromosome X and Y
final_peak_file <- file.path(macs2_output_dir, "final_peaks.bed")
cmd <- sprintf("grep -v -P '^chrX|chrY' %s > %s", 
               filtered_blacklist_removed_peak_file, 
               final_peak_file)
system(cmd)

# Intersect individual peak files with final peak set
final_peak_files <- lapply(peak_files, function(peak_file) {
  output_file <- file.path(macs2_output_dir, paste0("final_", basename(peak_file)))
  cmd <- sprintf("bedtools intersect -a %s -b %s > %s", peak_file, final_peak_file, output_file)
  system(cmd)
  return(output_file)
})

### Source of blacklist:
# https://github.com/Boyle-Lab/Blacklist

# 4. Load genomic regions

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoters <- promoters(genes(txdb))
five_utrs <- fiveUTRsByTranscript(txdb, use.names=TRUE)
three_utrs <- threeUTRsByTranscript(txdb, use.names=TRUE)
exons <- exonsBy(txdb, by="tx")
introns <- intronsByTranscript(txdb)
# Flatten the CompressedGRangesList objects into GRanges objects
five_utrs_flat <- unlist(five_utrs)
three_utrs_flat <- unlist(three_utrs)
exons_flat <- unlist(exons)
introns_flat <- unlist(introns)

# Filter hg19_seqlengths to include only canonical chromosomes
canonical_chromosomes <- c(paste0("chr", 1:22))
hg19_seqlengths <- seqlengths(txdb)
filtered_seqlengths <- hg19_seqlengths[canonical_chromosomes]

# Create the total GRanges object
hg19 <- GRanges(seqnames = names(filtered_seqlengths),
                IRanges(start = 1, end = as.integer(filtered_seqlengths)))

# Concatenate and then reduce to merge overlapping regions
#all_features <- reduce(append(append(append(append(promoters, five_utrs_flat), three_utrs_flat), exons_flat), introns_flat))
all_features <- unlist(GRangesList(promoters=promoters, five_utrs=five_utrs_flat, three_utrs=three_utrs_flat, exons=exons_flat, introns=introns_flat))


intergenic <- GenomicRanges::subtract(GRanges(seqnames=names(filtered_seqlengths), 
                                              IRanges(start=1, end=as.integer(filtered_seqlengths))),
                                      all_features)

region_list <- list(total=hg19, intergenic=intergenic, promoter=promoters, five_utr=five_utrs, exon=exons, intron=introns, three_utr=three_utrs)

# 5. Get final peak filenames and strip unnecessary parts for matching with BAM filenames

final_peak_basenames <- basename(unlist(final_peak_files))
stripped_peak_basenames <- gsub("^final_|_peaks\\.narrowPeak$", "", final_peak_basenames)

# 6. Match with conditions_df

bam_basenames_without_extension <- gsub("\\.bam$", "", conditions_df$filename)
matched_conditions <- conditions_df$condition[match(stripped_peak_basenames, bam_basenames_without_extension)]


# 7. Create a DataFrame for counting
# Define a function to count peaks in a region for a given set of peaks
count_peaks_in_region <- function(peaks, region) {
  olap <- findOverlaps(peaks, region)
  return(length(unique(queryHits(olap))))
}

all_counts_list <- lapply(final_peak_files, function(peak_file) {
  peaks <- import(peak_file)
  sapply(region_list, count_peaks_in_region, peaks = peaks)
})

all_counts_df <- as.data.frame(do.call(rbind, all_counts_list))
all_counts_df$Sample <- stripped_peak_basenames
all_counts_df$Condition <- matched_conditions

# 8. Melt data for ggplot
df_melted <- melt(all_counts_df, id.vars = c("Sample", "Condition"))

# 9. Plot
white_theme <- theme(
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  panel.border = element_blank(),
  panel.grid.major = element_line(color = "grey90"),  # Light gray major grid lines
  panel.grid.minor = element_line(color = "grey95"),  # Even lighter minor grid lines
  axis.title = element_text(face = "bold"),
  axis.text = element_text(color = "black"),
  axis.ticks = element_line(color = "black"),
  plot.title = element_text(face = "bold"),
  legend.background = element_rect(fill = "white"),
  legend.key = element_rect(fill = "white", color = "white")
)

ggplot(df_melted, aes(x = variable, y = log10(value + 1), fill = Condition)) + # Adding 1 to avoid log10(0)
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(y = "log10(Peak Counts)", x = "Region Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Distribution of Peak Counts in Different Regions by Condition") + white_theme

### Plot enrichment over genebodies
# Select only the first two BAM files of each condition for faster testing
selected_bam_files <- c(head(bam_files[bam_conditions == "1-2-3-4PMpos"], 2), head(bam_files[bam_conditions == "8-9-11PMneg"], 2))
selected_conditions <- c(rep("1-2-3-4PMpos", 2), rep("8-9-11PMneg", 2))

pm_pos_files <- selected_bam_files[selected_conditions == "1-2-3-4PMpos"]
pm_neg_files <- selected_bam_files[selected_conditions == "8-9-11PMneg"]

# 1. Extract Gene Regions with Buffer
buffer <- 2000
genes <- genes(txdb)
genes <- genes[!seqnames(genes) %in% c("chrX", "chrY")]

export(genes, con=file.path(macs2_output_dir, "genes_noXY.bed"), format="bed")

expanded_genes <- resize(genes, width=width(genes) + 2*buffer, fix="center")

# Loop through each bam file and run bamCoverage
for (each_bam in pm_neg_files) {
  output_file <- file.path(dirname(each_bam), sub("\\.bam$", ".bw", basename(each_bam)))
  cmd1 <- paste("bamCoverage --bam", each_bam, "-o", output_file, "\
-p max
--binSize 10
--normalizeUsing RPKM
--effectiveGenomeSize 2150570000
--ignoreForNormalization chrX
--extendReads")
  system(cmd1)
}

cmd1 <- paste("computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R",file.path(macs2_output_dir, "genes_noXY.bed"), "-S", 
              paste(pm_pos_files, collapse = " "), 
              "-o matrix_pmPos.gz -p 16")
system(cmd1)


cmd1 <- paste("computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R genes_noXY.bed -S", 
              paste(bam_files_1, collapse = " "), "-o matrix_1.gz")
system(cmd1)


# 2. Compute Coverage
compute_coverage <- function(bam_file, regions) {
  reads <- GenomicAlignments::readGAlignmentPairs(bam_file)
  subset_reads <- subsetByOverlaps(reads, regions)
  coverage <- coverage(subset_reads)
  return(coverage)
}

coverage_1_2_3_4PMpos <- lapply(pm_pos_files, compute_coverage, regions = expanded_genes)
# Normalizing function: Normalize each gene's coverage to have a sum of 1

coverage_8_9_11PMneg <- lapply(pm_neg_files, compute_coverage, regions = expanded_genes)

# Normalizing function: Normalize each gene's coverage to have a sum of 1
normalize_coverage_list <- function(coverage_list) {
  lapply(coverage_list, function(coverage_rlelist) {
    # Apply normalization for each Rle inside RleList
    normalized_rles <- endoapply(coverage_rlelist, function(rle) {
      total_counts <- sum(as.integer(rle))
      if (total_counts != 0) {
        rle / total_counts
      } else {
        rle
      }
    })
    return(RleList(normalized_rles, compress=FALSE))
  })
}


# Apply normalization
normalized_coverage_1_2_3_4PMpos <- normalize_coverage_list(coverage_1_2_3_4PMpos)
normalized_coverage_8_9_11PMneg <- normalize_coverage_list(coverage_8_9_11PMneg)

# Create a data frame for plotting
# Compute average profile
average_profile <- function(normalized_coverages) {
  # Convert RleList objects to a list of numeric vectors
  numeric_coverages <- lapply(normalized_coverages, function(coverage) {
    return(as.numeric(runValue(coverage)[runLength(coverage) == max(runLength(coverage))]))
  })
  
  # Compute the average for each position
  avg_profile <- colMeans(do.call(rbind, numeric_coverages), na.rm = TRUE)
  return(avg_profile)
}

avg_coverage_1_2_3_4PMpos <- average_profile(normalized_coverage_1_2_3_4PMpos)
avg_coverage_8_9_11PMneg <- average_profile(normalized_coverage_8_9_11PMneg) # Assuming you've normalized this condition as well

# Create a data frame for plotting
coverage_df <- data.frame(
  position = 1:length(avg_coverage_1_2_3_4PMpos),
  `1-2-3-4PMpos` = avg_coverage_1_2_3_4PMpos,
  `8-9-11PMneg` = avg_coverage_8_9_11PMneg
)

# Melt the data for ggplot2
coverage_melted <- tidyr::gather(coverage_df, condition, value, -position)

# Plot
ggplot(coverage_melted, aes(x = position, y = value, color = condition)) +
  geom_line() +
  labs(title = "Average Enrichment across Gene Bodies", x = "Position", y = "Normalized Coverage") +
  theme_minimal()
  
  
  
### SAVE CROSS TABULATION CHARTS
output_xtabs=0
if(output_xtabs == 1){
  # Convert meta columns to factors (except for condition and age)
  meta_factor <- conditions_df
  #meta_factor[cols] <- lapply(meta_factor[cols],factor)
  summary(meta_factor)
  
  #Save PlotXTabs2 https://www.rdocumentation.org/packages/CGPfunctions/versions/0.6.3/topics/PlotXTabs2 to a grid
  
  plot_sex <- PlotXTabs2(meta_factor,condition,sex,title="Sex")
  plot_race <- PlotXTabs2(meta_factor,condition,race, title="Race")
  plot_batch <- PlotXTabs2(meta_factor,condition,batch, title = "Batch")
  #plot_primary <- PlotXTabs2(meta_factor,condition,primary_present, title = "Primary Present")
  plot_site <- PlotXTabs2(meta_factor,condition,primary_site, title = "Primary Site")
  plot_chemo <- PlotXTabs2(meta_factor,condition,chemo_6weeks, title = "Chemo <6weeks")
  plot_age <- ggboxplot(meta_factor, x = "condition", y = "age", 
                        color = "condition", palette = c("#00AFBB", "#E7B800"),
ylab = "Age", xlab = "Condition", title = "Age") + stat_compare_means(method = "t.test")
#  include other mets plot if applicable, and add "plot_other" to ggarrange.
#  plot_other <- PlotXTabs2(meta_factor,condition,other_mets)

png(paste("./Output/DESeq2/Xtabs/",contrast_groups[2],contrast_groups[3],"_",ver,".png", sep = ""),
    width=2600, height=1000, res = 144)
ggarrange(plot_sex, plot_race, plot_batch, plot_chemo, plot_site, plot_age, ncol = 3, nrow = 2)
#            labels = c("Sex", "Race", "Batch", "Primary Present","Primary Site", "Chemo 6 Weeks", "Age", "Other Mets (n/a"),
#            vjust	= 0.5,
dev.off()

svg(paste("./Output/DESeq2/Xtabs/",contrast_groups[2],contrast_groups[3],"_",ver,".svg", sep = ""),
    width=2600, height=1000)
ggarrange(plot_sex, plot_race, plot_batch, plot_chemo, plot_site, plot_age, ncol = 3, nrow = 2)
dev.off()

}