#Froniters Paper NMF R Script
#Running NMF on gene counts from 56 OC cell lines to get consensus clusters

#Author: Aideen McCabe
#Affiliation: School of Biochemistry and Cell Biology, University College Cork, Cork, Ireland

####################################################################################################
### ------------------------------- library --------------------------------- ###

library(readr)
library(dplyr)
library(DESeq2)
library(circlize)
library(ComplexHeatmap)
library(NMF)
library(sva)

####################################################################################################
### ------------------------------- Read data ------------------------------- ###

counts <- read_csv("~/PhD/Year2/Data/RNASeq/CellLines/cell_line_counts_canceronly.csv")
coldata <- read_csv("~/PhD/Year2/Data/RNASeq/CellLines/cell_line_coldata_canceronly.csv")


####################################################################################################
### ------------------------------- Tidy data ------------------------------- ###

#Editing coldata and counts to be in the correct format for DESeq2

#Add rownames to coldata
rows <- coldata$name
rownames(coldata) <- rows

#Create a copy of counts data
all_counts <- counts

#Remove gene column and use as rownames instead
all_counts <- all_counts %>% dplyr::select(-gene)
rownames(all_counts) <- counts$gene
colnames(all_counts) <- rownames(coldata)

#Convert counts to matrix
counts <- as.matrix(all_counts)

#Quality check, ensure columns and rows match in the same order
all(rownames(coldata) == colnames(counts))

####################################################################################################
### ------------------------------- Data normalization ---------------------- ###

#Use DESeq2 to apply a variance stabilizing transformation  
vsd <- vst(counts, blind = T)

#Calculate median absolute deviation for each gene/row
mads <- apply(vsd, 1, mad)

#Select for genes with a median absolute deviation of >=1.5 (most variability)
mads_filtered <- mads[mads >= 1.5]

#Convert resulting vector of MAD for each gene to a data frame
mad_df <- as.data.frame(mads_filtered)

#Filter original matrix to retain only transcripts with mad >=1.5
filtered_vsd <- vsd[rownames(vsd) %in%  rownames(mad_df),]

####################################################################################################
### -------------------------------Batch Correction ------------------------- ###

#Specify which datasets each cell line originated from (column in coldata)
batch <- coldata$batch

#Create a null model since we are not adjusting for any other variables
modcombat <- model.matrix(~1, data = coldata)

#Carry out batch correction on normalized counts using ComBat
adjusted_counts <- ComBat(filtered_vsd, batch = batch, mod = modcombat, par.prior = T, prior.plots = F)


####################################################################################################
### -------------------- Estimating Factorization Rank (k) ------------------ ###

#Run NMF on normalized, filtered and batch corrected counts with k 2 to 8
#Specifying 50 runs and a seed
#NMF algorithm is brunet
estim.r.batch<- nmf(adjusted_counts, 2:8, method = "brunet", nrun = 50, seed = 123456)

#Shuffle the original data 
V.random.batch <- randomize(adjusted_counts)

#Estimate quality metrics from the shuffled data (to avoid overfitting)
estim.r.random.batch <- nmf(V.random.batch, 2:8, method = "brunet", nrun = 50, seed = 123456)

#Plot metrics
plot(estim.r.batch)
plot(estim.r.random.batch)

plot(estim.r.batch, estim.r.random.batch)
#k = 5 displays high quality metrics


#Calculate smallest value of r for which decrease in residuals is greater
#in data than decrease in random data
for(i in 1:6){
  diff_data <- estim.r.batch$measures$residuals[i]- estim.r.batch$measures$residuals[i+1]
  print(diff_data)
  
  diff_rand <- estim.r.random.batch$measures$residuals[i] - estim.r.random.batch$measures$residuals[i+1]
  print(diff_rand)
  
  print(diff_data>diff_rand)
  
  if(diff_data<diff_rand){
    print(i)
  }
}

#5 is optimal number of clusters, this will be used as factorization rank

####################################################################################################
### ----------------------------- Run NMF with k=5 -------------------------- ###

##Run nmf with 5 clusters
celllines_nmf_batch <- nmf(adjusted_counts, 5,  method = "brunet", nrun=200, seed = 198913)

#Consensus map for k=5
consensusmap(celllines_nmf_batch, 
             labCol=colnames(adjusted_counts), labRow=colnames(adjusted_counts))


#Basis
W <- basis(celllines_nmf_batch)


#Coefficients
H <- coef(celllines_nmf_batch)

#Find cluster with max coefficient score
clusters <-apply(H, 2, which.max)

#View clusters
clusters

