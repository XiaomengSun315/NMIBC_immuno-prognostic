rm(list=ls())
################################## BC progression ##################################
# install.packages("BiocManager")
# BiocManager::install(update=TRUE, ask=FALSE)
# BiocManager::install(version='3.13')
# BiocManager::install("impute", force=TRUE)
# BiocManager::install("PCAtools")
library(openxlsx)
library(preprocessCore)
library(impute)
library(RUVSeq)


# set working directory
data_dir <- "C:/0_xmsun/xmsun/Graduate/20210224_NMIBC/Data/"
result_dir <- "C:/0_xmsun/xmsun/Graduate/20210224_NMIBC/Results/"
dir.create(paste0(result_dir, "91_reviewer1_norm/"))


########## Function: PCA & RLE plot
pca_plot <- function(label, exp_data, pheno, save_path){

	########## data processing
	# expression matrix: genes in row, samples in column, numeric format
	exp_matrix <- apply(exp_data, 2, as.numeric)
	rownames(exp_matrix) <- rownames(exp_data)
	exp_matrix <- na.omit(exp_matrix)	
	exp_matrix_log <- log2(exp_matrix+1) # log2 transform

	# phenotype matirx: samples in row, phenotypes in column
	pheno_matrix <- pheno[match(colnames(exp_matrix_log), rownames(pheno)),]
	pheno_matrix$Source_dataset <- as.factor(pheno_matrix$Source_dataset)
	pheno_matrix$Platform <- as.factor(pheno_matrix$Platform)
	pheno_matrix$Staging <- factor(pheno_matrix$Staging, levels=c("T0", "Ta", "Ta-1", "T1", "Tis"))

	########## PCA plot
	p <- pca(exp_matrix_log, metadata=pheno_matrix, removeVar=0.1)
	scree_pca <- screeplot(p, components = getComponents(p, 1:20), axisLabSize = 18, titleLabSize = 22)
	pair_pca_platform <- pairsplot(p, components = getComponents(p, 1:3), colby='Platform', colkey=c("Microarray"="#d53e4f", "RNA-Seq"="#3288bd"))
	pair_pca_dataset <- pairsplot(p, components = getComponents(p, 1:3), colby='Source_dataset')
	pair_pca_stage <- pairsplot(p, components = getComponents(p, 1:3), colby='Staging', colkey=c("T0"="#d53e4f", "Ta"="#fc8d59", "Ta-1"="#e6f598", "T1"="#99d594", "Tis"="#3288bd"))

	pdf(paste0(save_path, label, ".pdf"))
	print(scree_pca)
	print(pair_pca_platform)
	print(pair_pca_dataset)
	print(pair_pca_stage)
	dev.off()
}

######################################### main #######################################
########## load data
# my original method, matrix included: exp_data_all, exp_data_all_house, exp_data_enough, exp_data_all_house_knn, exp_data_all_house_knn_quant, exp_data_all_house_knn_quant_log.
load(file=paste0(result_dir, "7_Housekeeping_normalization/7.99_ALL_norm_process.Rdata"))
# original phenotype data
pheno <- read.xlsx(paste0(data_dir, "combined/Clinicopathologic_v2.xlsx"), na.strings=c("","NA","Nan","#N/A"))
rownames(pheno) <- pheno$Sample_name
# RUVg, matrix included: set, set1, set2
load(file=paste0(result_dir, "7_Housekeeping_normalization_v3/7.99_RUVg_norm_process.Rdata"))

# PCA plots
pca_plot("0_non_norm", exp_data_all, pheno, paste0(result_dir, "91_reviewer1_norm/"))
pca_plot("1_house_norm", exp_data_all_house, pheno, paste0(result_dir, "91_reviewer1_norm/"))
pca_plot("2_xmsun_norm", exp_data_all_house_knn_quant, pheno, paste0(result_dir, "91_reviewer1_norm/"))
pca_plot("3.0_non_norm", counts(set), pData(set), paste0(result_dir, "91_reviewer1_norm/"))
pca_plot("3.1_RUVg_house", normCounts(set1), pData(set1), paste0(result_dir, "91_reviewer1_norm/"))
pca_plot("3.2_RUVg_empirical", normCounts(set2), pData(set2), paste0(result_dir, "91_reviewer1_norm/"))
