#Purpose of this script is to visualize differential gene expression heatmaps,
#PCA plot and other charts between two comparison groups.

### INSTALL LIBRARIES
# Setup, uncomment follow
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("augment")
# BiocManager::install("ggiraphExtra")
#devtools::install_github("cardiomoon/ggiraphExtra")

library(DESeq2) #for DESeq2 analysis
library(magrittr)
library(tibble)
library(dplyr)
library(tidyverse)
library(ggplot2) #for plotting
library(ggrepel)
#library(DEGreport)
library(RColorBrewer)
library(pheatmap) #for heatmaps
library(viridis)
library(colorspace)
#library(M3C)
library(ashr)
library(CGPfunctions) #needed for plotxtabs2 
library(ggpubr) # needed for whisker plot charting
library(aod) #for logistic regression
library(broom)
library(leaps) #for logistic regression optimizaton
library(pROC)
library(boot)
#library(ggfortify) #used for pca clustering
library(cluster) #used for pca clustering
#library(ggiraph)
#library(ggiraphExtra)
library(devtools) #used for complex heatmaps
library(ComplexHeatmap)#used for complex heatmaps
library(circlize)

### Raw counts
#counts_name_training <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_v3/1o3_pdONLYpos_8o11_metNEGtumNEG_training_rawcounts.csv"
#counts_name_validation <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_v3/1o3_pdONLYpos_8o11_metNEGtumNEG_validation_rawcounts.csv"

### Red in counts (Note: randomization is raw counts; script below normalizes raw counts.)
counts_name_training <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_final/1o3_pdONLYpos_8o11_metNEGtumNEG_training_rawcounts.csv"
meta_name_training <- "./Output/Randomization/1o3_pdONLYpos_8o11_metNEGtumNEG_DESeq2_final/1o3_pdONLYpos_8o11_metNEGtumNEG_training_conditions.csv"
counts_data_training <- read.csv(counts_name_training,row.names = 1)
meta_training <-  read.csv(meta_name_training,row.names=1)

# Correct for forbidden characters
rownames(counts_data_training) <- str_replace(c(rownames(counts_data_training)),'\\-',"_")
rownames(counts_data_training) <- str_replace(c(rownames(counts_data_training)),'\\.',"_")

#create conditions vector expected as input by glmnet (0,0,0,...,1,1,1,,...)
sample_count_training <- nrow(meta_training)
class1_name <- meta_training[2,1]
class2_name <- meta_training[nrow(meta_training),1]
class1_count_training <- sum(meta_training$condition == class1_name)
class2_count_training <- sum(meta_training$condition == class2_name)
conditions_vector=c(rep(0,class1_count_training),rep(1,class2_count_training))

###CONFIGURATION
#set working directory, select where you extracted folder
setwd("~/5hmC-Pathway-Analysis/")

#include to exclude specific sample
#counts_data <- subset(counts_data,select=-c(KT126))
#meta_training <- meta_training[!(row.names(meta_training) %in% c("KT126")),]

ver <- "training_genome_wide_vFINAL"
gene_number <- nrow(counts_data_training)

#Set desired outputs:
output_results_tables = 1
output_heatmap = 1
output_PCA = 1
output_config = 1

# define contrast groups
groups <- unique(meta_training[c("condition")])
contrast_groups <- c("condition",groups[1,1], groups[2,1])

if (file.exists(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",sep=""))) {
  cat("The folder already exists")
} else {
  dir.create(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",sep=""))
}

###VALIDATION
#check columns are equal
all(colnames(counts_data_training) %in% rownames(meta_training))
all(colnames(counts_data_training) == rownames(meta_training))

###CREATE DESeq2 OBJECT
#load data into DESeq2 object dds
design_formula <- ~ condition #sex + age + race + batch + condition 

dds <- DESeqDataSetFromMatrix(countData = counts_data_training, colData = meta_training, design = design_formula)

#load up size factors into dds object in order to normalize using median of ratios method
dds <- DESeq(dds)

#use counts 'normalized=true' function to pull out normalized counts
normalized_counts <- counts(dds,normalized=TRUE)

normalized_counts_tb <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

nsamples = nrow(meta_training)

normalized_counts_tb <- normalized_counts_tb[,c(1,2:nsamples)] %>% 
  data.frame() %>%
  column_to_rownames(var = "gene") 

###Genome wide regression

#determine number of samples and significant genes
nsamples <- ncol(normalized_counts_tb)
nsig <- nrow(normalized_counts_tb)

###LOGISTIC REGRESSION
#combine metadata and counts data

sig_gene_names <- normalized_counts_tb$gene
meta_sig_colnames = c(colnames(meta_training),sig_gene_names)
predictors_training = merge(meta_training,t(normalized_counts_tb),by.x=0,by.y=0,all.x=TRUE, all.y=TRUE)

#convert meta columns to factors
predictors_training[colnames(meta_training)] <- lapply(predictors_training[colnames(meta_training)],factor)

#replace outcome with 1=PM pos and 0 =PM neg (confirm that contrfast groups match)
predictors_training$condition = gsub(contrast_groups[2],'1',predictors_training$condition)

predictors_training$condition = gsub(contrast_groups[2],'1',predictors_training$condition)
predictors_training$condition = gsub(contrast_groups[3],'0',predictors_training$condition)

# remove "-" from gene names which cause problems in gene formula
colnames(predictors_training) <- gsub ("-","_",colnames(predictors_training))

predictors_training$condition <- as.numeric(predictors_training$condition)
predictors_training$age <- as.numeric(predictors_training$age)
predictors_training <- subset(predictors_training,select=-c(peritoneal_mets,other_mets,primary_present))
#predictors_training <- subset(predictors_training,select=-c(peritoneal_mets))


# Set the "Row.names" column as row names
rownames(predictors_training) <- predictors_training$Row.names

# Remove the "Row.names" column from the dataframe; and count number of predictors
predictors_training <- subset(predictors_training, select = -Row.names)
predictor_count <- ncol(predictors_training)

head(colnames(predictors_training),20)
head(rownames(predictors_training),20)


# Plot all predictors individually
for (single_predictor_index in 1:predictor_count) {
#  logit_variables = paste(colnames(predictors[-1]), collapse=" + ")
  logit_variables = colnames(predictors_training[single_predictor_index])
  logit_formula = as.formula(paste("condition ~ ",logit_variables, collapse=""))
  
  mylogit <- glm(logit_formula, data = predictors_training, family = "binomial")
  logistic_results <- summary(mylogit)
  #print(single_predictor_index)
  #print(summary(mylogit))
  if (!exists("individual_results")){
    individual_results = tidy(mylogit)
    } 
  else{
    individual_results <- rbind(individual_results,tidy(mylogit))
    #print(tidy_glm)
  }
  if(single_predictor_index == predictor_count){
    individual_results = individual_results[!(individual_results$term == "(Intercept)"),]
    write.csv(individual_results, paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/","individual_predictors.csv",sep=""))
  }
}

###Training heatmap
#Save heatmap based on logistic regression results:
#log_counts_data<- rownames(counts_data_training_subset)
heatmap_pvalue_cutoff = 0.01

log_counts_data <- counts_data_training
rownames(log_counts_data) <- str_replace(c(rownames(log_counts_data)),'\\.',"_")
rownames(log_counts_data) <- str_replace(c(rownames(log_counts_data)),'\\-',"_")
heatmap_title <- paste(contrast_groups[2],"/",contrast_groups[3],"logistic regression sig genes p-value <0.01")

#transpose
#log_counts_data <- t(log_counts_data) #uncomment for transpose

annotation <- meta_training %>% 
  select(condition)

# Identify susbset of genes from genome-wide logitic regression with pvalue <  cutoff
gene_names_u_cutoff = subset(individual_results, p.value < heatmap_pvalue_cutoff)$term

log_counts_data_mat = data.matrix(log_counts_data[rownames(log_counts_data) %in% gene_names_u_cutoff,], rownames.force = NA)
log_counts_data_mat = t(scale(t(log_counts_data_mat)))
ha = HeatmapAnnotation(condition=annotation$condition,col = list(condition = c("1o3_pdONLYpos" = "#ff9289", "8o11_metNEGtumNEG" = "#00dae0")))

svg(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"training_heatmap_eucl",".svg",sep=""))
#png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"training_heatmap_eucl",".png",sep=""), width = 900, height = 600)
#pheatmap(log_counts_data_mat,
#         main = heatmap_title,
##         clustering_distance_cols = "binary",
#         cluster_rows = T,
#         cluster_cols = T,
#         show_rownames = T,
#         top_annotation = ha,
##         annotation_row = annotation,  #uncomment for transpose
#         border_color = NA, 
#         fontsize = 10, 
#         scale = "row", 
##         scale = "column", #uncomment for transpose
#         fontsize_row = 10, 
#         height = 20
#)
Heatmap(log_counts_data_mat, 
        top_annotation = ha,
        #col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
        col = colorRamp2(c(-4, -2,-1,0,1, 2,4), c("#4676b5", "#82b1d3","#dbeff6","#fefebd","#fee395", "#fc9961","#d73027")),
        clustering_distance_rows = "euclidean",
        cluster_rows = TRUE,
        clustering_distance_columns = "euclidean",
)
dev.off()


###Whole dataset heatmap
#Save heatmap based on logistic regression results:
#log_counts_data<- rownames(counts_data_training_subset)
log_counts_data <- normalized_counts
rownames(log_counts_data) <- str_replace(c(rownames(log_counts_data)),'\\.',"_")
rownames(log_counts_data) <- str_replace(c(rownames(log_counts_data)),'\\-',"_")
log_counts_data <- log_counts_data[c(rownames(counts_data_training_subset)),]
heatmap_title <- paste(contrast_groups[2],"/",contrast_groups[3],"logistic regression sig genes p-value <0.01")

annotation <- meta %>% 
  select(condition)

#transpose
#log_counts_data <- t(log_counts_data)

png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/",contrast_groups[2],contrast_groups[3],ver,"wholeset100_def",".png",sep=""), width = 900, height = 1200)
pheatmap(log_counts_data,
         main = heatmap_title,
         #         clustering_distance_cols = "binary",
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         annotation = annotation,
         #         annotation_row = annotation,  #uncomment for transpose
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         #         scale = "column", #uncomment for transpose 
         fontsize_row = 10, 
         height = 20,
)
dev.off()

###GENERATE PCA PLOT
#See details for below operation in lesson 3 of DGE workshop
rld <- vst(dds, blind=TRUE)
res_genes <- row.names(res_table)
sig_genes <-  sig$gene
#p < 0.05
#sig_genes <- c("KRT34","SCGB1D4","FERD3L","HEXIM1","CLEC4M","OR5AC2","OR1L4","KRT19","KIF4B","IL10","OR1D2","TMEM236","CD1C","WNT8A","RPL29","APOF","HLA_DMB","LAIR2","METTL7B","PGLYRP3","KRTAP19_6","RPE65","ZNF816","CTD_2287O16_3","AL359195_1","CTD_2545G14_7","SBDS","CD300E","FKBP9","FCER1G","PCYT1A","PSAP","HCST","CTSS","IQCF5","TAAR5","ZNF572","SLC35B1","LGALS13","CD300LB","KRTAP21_2","PCDHB4","PRPS1L1","CLPS","S100A11","HIST1H4D","CAMP","MARCH1","KRT31","HRH4","CWF19L1","RP11_793H13_10","MS4A6E","FKBP1B","REG3A","OR7G3","PAX3","COMTD1","C21orf62","LRRC8C","FCRL5","HES3","ZNF134","ADGRE3","NLRC4","CLEC4G","PILRA","RP4_576H24_4","RP11_444E17_6","IL23A","NAIP","CEBPB","LITAF","LILRA4","SCGB1D2","LTB4R","SPRR2D","FOS","OTOP1","GJB7","RPL3","CTC_435M10_3","KLK3","HPD","SLAMF9","SLC9A4","STX10","LA16c_431H6_6","ANGPT4","TNFSF14","SPCS1","ARHGAP30","CXCR1","ITK","C6orf226","ITGAM","MEFV","LCE2C","SGPP2","SEMG2","CYP11B2","ISG20L2","C1orf54","LAMTOR4","HRH2","SLC22A14","ARRDC5","SV2A","MS4A1","ZNF284","CRIP3","SLC46A2","OR7E24","DDIT3","MS4A13","WDR46")
#sig_genes <- c("KRT34","SCGB1D4","FERD3L","HEXIM1","CLEC4M","OR5AC2","OR1L4","KRT19","KIF4B","IL10","OR1D2","TMEM236","CD1C","WNT8A","RPL29","APOF","HLA_DMB","LAIR2","METTL7B","PGLYRP3","KRTAP19_6","RPE65","ZNF816","CTD_2287O16_3","AL359195_1","CTD_2545G14_7","SBDS","CD300E","FKBP9","FCER1G","PCYT1A","PSAP","HCST","CTSS","IQCF5","TAAR5","ZNF572","SLC35B1","LGALS13","CD300LB","KRTAP21_2","PCDHB4","PRPS1L1","CLPS","S100A11","HIST1H4D","CAMP","44621","KRT31","HRH4","CWF19L1","RP11_793H13_10","MS4A6E","FKBP1B","REG3A","OR7G3","PAX3","COMTD1","C21orf62","LRRC8C","FCRL5","HES3","ZNF134","ADGRE3","NLRC4","CLEC4G","PILRA","RP4_576H24_4","RP11_444E17_6","IL23A","NAIP","CEBPB","LITAF","LILRA4","SCGB1D2","LTB4R","SPRR2D","FOS","OTOP1","GJB7","RPL3","CTC_435M10_3","KLK3","HPD","SLAMF9","SLC9A4","STX10","LA16c_431H6_6","ANGPT4","TNFSF14","SPCS1","ARHGAP30","CXCR1","ITK","C6orf226","ITGAM","MEFV","LCE2C","SGPP2","SEMG2","CYP11B2","ISG20L2","C1orf54","LAMTOR4","HRH2","SLC22A14","ARRDC5","SV2A","MS4A1","ZNF284","CRIP3","SLC46A2","OR7E24","DDIT3","MS4A13","WDR46","NXT1","PILRB","S100A12","SLC7A7","IL10RA","CST11","ZNF729","RNF175","KRTAP4_5","MS4A4A","OR5B17","AURKB","GATS","OR6P1","GJB6","CASP14","PENK","C19orf45","TMEM154","KIAA0226L","PADI2","KRTAP26_1","AC107021_1","HK3","AJUBA","NPIPA3","OR10J3","ZNF320","ACP6","TIGD2","LY75_CD302","STARD10","CLEC4E","RRS1","SLC28A3","CAPG","CEACAM4","RP11_257K9_8","MS4A6A","TMEM38A","chemo_6weeksYes","MYO1F","HCRT","SARS2","SIGLECL1","SLC2A14","WDR77","IL1RN","GAS7","CALB2","MRFAP1L1","KRT32","CRCP","ZNF333","APOBEC3A","BCL2L15","CXCR4","EXD1","CEACAM8","RAB39A","ZNF397","SPINT4","MNX1","S100A9","TIFAB","S100A1","CD200R1L","C1orf194","HSH2D","HES1","PCDH17","PHKG1","MPEG1","LRRD1","FOXE3","DEF8","HIST1H2AG","OR2B2","OR7G2","SP140","TMEM50B","PDPR","CCR3","KLF6","CTSK","AC010547_9","G6PC","MSRB1","KIF15","LGALS9B","RP11_599B13_6","IL36G","DUS2","MTNR1A","EEF1A1","IFNL2","OR2T10","RHOG","ADORA2A","KRT38","NUPR1","MZF1","ZNF486","ID1","HCK","MARCO","APMAP","RFPL3","ASIP","GOLGA8O","APEX1","RP11_49K24_6","MRM1","CX3CR1","EIF5A","CHRNA6","NLRP12","TYROBP","LYPD6B","LYZ","CD53","TREM1","MS4A14","FOLR3","C19orf38","CCR8","GRPEL1","TCP10L2","NWD1","CNR2","SMPDL3B","TUBA4A","SLC31A2","ALAS1","TKT","OR10G6","TMEM180","DHRS9","FOXB2","SYTL3","SLC26A8","CD14","HTR3A","JUN","MANF","CD209","SUPT4H1","MS4A5","AK5","TMEM45B","MT4","IFNA6","RAB31","TVP23A","FCRL6","S100A8","SERPINB8","C5orf38","TAS2R5","TOP3A","CKAP4","CHGB","KRTAP5_9","TMEM130","ADRB2","SLC18A1","DNAJC5","TRAPPC5","POLE3","RPL18","OR6V1","HOXA10","SNX11","CCDC125","CLEC17A","FPR1","SMN2","ID4","KRT33A","NLRP2","HLX","GNS","LILRB2","C16orf70","ARHGAP9","ALOX15","CLDN8","COX7C","COX6B1","SPRR1B","GABRR2","CCDC140","PRR13","OCM","MRPL11","ACOXL","FAM72A","FGD3","NGB","C14orf144","MRFAP1","CD44","MPZL3","FMOD","TAGAP","DEFB107B","SERPINA4","LHX4","COMMD3_BMI1","PROKR2","OR2V1","OR51E1","IL32","NDUFS3","OR14K1","GLI1","GLTP","HTR3B","RP4_583P15_15","NLRP3","EEF2KMT","CDHR2","MED11","RDH16","GMFG","LEP","C1orf162","OR1E1","POLR2G","CPLX2","RAB11FIP1","COPS6","AC018755_1","LAIR1","KIF22","AP001024_2","CEP97","SLC39A11","CCL27","PGD","FAM105A","ARSG","FLT3","AK2","FCAR","ICAM3","MAPKAPK3","FCRLB","FGF6","IFI27L2","SLC2A5","CDA","HIST1H3B","SEC31B","E2F2","NEURL2","SELV","NSUN5","LY75","KRT5","TMEM191B","TSTD1","KRT36","CEACAM3","SLC22A15","FCGR3A","MGARP","NFKBID","TMEM52B","KAT8","MS4A18","LTA4H","RASGRP1","RAB37","FABP1","RXRB","FAM27E2","SIRPD","SIGMAR1","CSNK1A1L","ARIH2OS","SURF1","CD48","DNAJC5B","OR52I1","NCR3LG1","TNP1","DAGLB","SERPINA1","PUS3","SNX10","PLAC8","DIO2","MPPE1","SLC35F2","LGALS3","PRB2","ADPRH","OR6X1","LINC00694","C10orf12","ZAR1L","NCR3","CEACAM1","HOXC10","GLIPR2","RAB3D","LMNB1","ECH1","FAM214B","CD68","PDAP1","FCGBP","MAP10","AP4B1","TMEM81","SLC22A7","SEMA4A","KRTAP19_1","TXNL4B","TDRKH","FCGR2B","DRC3","NUDT19","A4GNT","ZNF813","VAV1","RPL41","RRP12","PLXNC1","CTD_2116N17_1","OLAH","NFKBIA","RPL26L1","POLR2J3","PKD2L1","CYSLTR2","ZSCAN26","BORCS7","AADACL3","ZNF788","PFKFB4","OR2Z1","AL669831_1","TRADD","DPAGT1","RUNX1T1","HIGD1B","MYOF","FOXM1","PATE3","RBM8A","PBX2","OR6C70","PLA2G2D","PPP1R17","HAVCR2","THOC3","BRI3BP","CFL1","TMPRSS4","PCGF1","OR10G8","S100P","TTC14","APBB1IP","CORO2A","CBR3","PTAFR","TAS2R39","ADAT1","RHOH","PRMT5","PGLYRP4","COLGALT2","AICDA","C12orf76","CTB_54O9_9","C1orf106","UHRF1","RP3_369A17_6","PSMB8","IZUMO3","GMIP","CATIP","CCDC74A")
#sig_genes <- c("RPL3","ZNF134")
rownames(rld) <- str_replace(c(rownames(rld)),'\\.',"_")
rownames(rld) <- str_replace(c(rownames(rld)),'\\-',"_")

#save PCA plot to png
#In the below replace sig_genes with res_genes if you want to perform PCA analysis on all genes rather than just on significant genes.
png(paste("./Output/DESeq2/logit/",contrast_groups[2],contrast_groups[3],"_",ver,"/","log_sig_PCA_",contrast_groups[2],contrast_groups[3],ver,".png",sep=""), width = 900, height = 1200)

plotPCA_labels <- plotPCA(
  rld[sig_genes,], intgroup = c("condition")
  #rownames(c(row.names(meta)))
) + coord_fixed(ratio=1) + stat_ellipse(level = 0.9)
 
#reverse order of labels to match heatmap labelling
plotPCA_labels$data$group<-factor(plotPCA_labels$data$group,levels = rev(levels(plotPCA_labels$data$group)))
plotPCA_labels$data$condition<-factor(plotPCA_labels$data$condition,levels = rev(levels(plotPCA_labels$data$condition)))

plotPCA_labels 
#+ geom_text(aes(label = name),position=position_nudge(y = 0.07),) + ggtitle(heatmap_title) 

dev.off()

