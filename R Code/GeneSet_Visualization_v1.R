#PURPOSE OF THIS SCRIPT IS TO VISUALIZE DGE HEATMAPS & PCA PLOTS BETWEEN TWO COMPARISON GROUPS
#tutorial taken from here: https://github.com/hbctraining/DGE_workshop/tree/master/lessons

###LIBRARIES INSTALL ALL IF FIRST TIME RUNNING
#install.packages("package_name")
library(DESeq2)
library(magrittr)
library(tibble)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(colorspace)
library(M3C)
library(ashr)

###CONFIGURATION
#set working directory, select where you extracted folder
setwd("C:/Users/ymali/Google Drive/Personal Documents/Chuan Lab/Peritoneal Disease/Data Analysis/ssGSEA")


data_name <- "./Outputs/ssGSEA_sigonly_tpm_v3.csv"
meta_name <- "./Inputs/conditions_pm_v1.csv"
heatmap_data <- read.csv(data_name,row.names = 1)
#read in data, define what counts & conditions files
meta <-  read.csv(meta_name,row.names = 1)
meta
heatmap_data



#Select version for all output files (e.g. 1, 2, 3, ...)
ver <- "ssGSEA_sigonly_tpm_v2"



###GENERATE HEATMAP
#Annotate our heatmap (optional)
annotation <- meta %>% 
  select(condition)

annotation

#Save heatmap to png
heatmap_title <- paste("PM+ / PM- ssGSEA","FDR(BH) < 0.2")
png(paste("./Outputs/sig_heatmap_ssGSEA_sigonly","_v",ver,".png", sep = ""), width = 1500, height = 1200)
pheatmap(heatmap_data, 
         main = heatmap_title,
         color = diverging_hcl(15,"Blue-Red2"), 
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20,
)
dev.off()

