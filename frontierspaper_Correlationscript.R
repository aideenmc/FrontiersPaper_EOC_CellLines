#Froniters Paper Correlation Analysis R Script
#Calculating correlation of gene expression between cell lines and primary tumors

#Author: Aideen McCabe
#Affiliation: School of Biochemistry and Cell Biology, University College Cork, Cork, Ireland

####################################################################################################
### ------------------------------- library --------------------------------- ###

library(readr)
library(edgeR)
library(sva)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(corrplot)

####################################################################################################
### ------------------------------- Read data ------------------------------- ###

#Count matrix
cells_tissue_counts <- read_csv("~/PhD/Year2/Data/RNASeq/CellLines/cells_tissue_counts_noearly.csv")

#Associated metadata
cells_tissue_coldata <- read_csv("~/PhD/Year2/Data/RNASeq/CellLines/cells_tissue_coldata_noearly.csv")

#Isolate cell and tissue coldata seperately
cell_coldata <- cells_tissue_coldata[92:147,]

tissue_coldata <- cells_tissue_coldata[c(1:91),]

#2233 most variable transcripts from NMF analysis
mostvariable_transcripts <- read_delim("~/PhD/Year2/Data/RNASeq/CellLines/nmf_genes.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

####################################################################################################
### ------------------------------- Tidy data ------------------------------- ###

#Store gene names in vector
rows <- cells_tissue_counts$gene

#Remove gene column
cells_tissue_counts <- cells_tissue_counts[,-c(1)]

#Set rownames
rownames(cells_tissue_counts) <- rows

####################################################################################################
### ------------------- Normalization and batch correction ------------------ ###

#Create a digital gene expression data object from count matrix
my_DGElist <- DGEList(counts = cells_tissue_counts)

#Calculate scaling factors to convert raw library sizes into effective library sizes
my_DGElist_norm <- calcNormFactors(my_DGElist, method = 'upperquartile')

#Upper quartile normalization
normalized <- cpm(my_DGElist_norm, log = T, prior.count = 1)

#Batch correction
#Specify which datasets each cell line originated from (column in coldata)
batch <-  cells_tissue_coldata$batch

#Model matrix adjusting for subtype covariate
modcombat <- model.matrix(~subtype, data = cells_tissue_coldata)

#Batch correct counts
adjusted_counts <- ComBat(normalized, batch = batch, mod = modcombat, par.prior = T, prior.plots = F)

#Filter adjusted counts to retain most variable transcripts only
adjusted_filtered_counts <- adjusted_counts[rownames(adjusted_counts) %in% mostvariable_transcripts$x,]

####################################################################################################
### ------------------------- Correlation Analysis -------------------------- ###

#Specify cells
cells <- adjusted_filtered_counts[,c(92:147)]

#Specify tissue
tissues <- adjusted_filtered_counts[,c(1:91)]

#Calculate Spearman correlation between cells and tissues
cormat_combo <- cor(tissues, cells, method = "s")

#Max and min values
max(cormat_combo)
min(cormat_combo)

#Load Domcke and Yu's HGSOC rankings to calculate Spearmans Rho with ranking generated here
domcke_yu_correlation_ranks <- read_csv("~/PhD/Year2/Writing/FrontiersPaper/domcke_yu_correlation_ranks.csv")

#Domcke
corr_domcke <- cor.test(x=domcke_yu_correlation_ranks$Domcke_rank, y=domcke_yu_correlation_ranks$RevisedRank, method = 'spearman')

corr_yu <- cor.test(x=domcke_yu_correlation_ranks$Yu_rank, y=domcke_yu_correlation_ranks$RevisedRank, method = 'spearman')
  
corr_domcke
corr_yu

