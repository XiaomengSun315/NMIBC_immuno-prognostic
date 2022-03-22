# modified according to reviewer1's suggestion: update with "6_Gene_Survival_v2.R" and "7_Housekeeping_normalization_v3.R"

rm(list=ls())
################################## BC progression ##################################
library(openxlsx)
library(ComplexHeatmap)
library(GSVA)
library(circlize)
library(RUVSeq)

# set working directory
data_dir <- "C:/0_xmsun/xmsun/Graduate/20210224_NMIBC/Data/"
result_dir <- "C:/0_xmsun/xmsun/Graduate/20210224_NMIBC/Results/"
dir.create(paste0(result_dir, "8_Cell_scores/"))

################################## Cell Score Heatmap ##################################

Cell_scores_heatmap <- function(dataset_name, top_genes, top_genes_convert){

	# expression data (without log)
	# load re-normalized matrixs: set (non-norm), set1 (seven housekeeping genes norm), set2 (empirical genes norm)
	sample_names <- names(read.csv(paste0(result_dir, "1_data_clean/", dataset_name, "_clean.txt"), sep="\t", row.names=1))
	exp_data <- normCounts(set2)[,colnames(normCounts(set2)) %in% sample_names]

	# top genes expression data
	exp_data_hit_nolog <- data.frame(exp_data[match(top_genes$gene_name, rownames(exp_data)),])
	exp_data_hit_log <- log2(exp_data_hit_nolog + 1)

	# 1) ssGSEA gene expression Heatmap (with/without log2)
	exp_DEG <- as.matrix(exp_data_hit_log)
	ht <- Heatmap(exp_DEG, na_col="grey", cluster_rows = FALSE, cluster_columns = TRUE, show_row_dend = FALSE, show_row_names = TRUE, show_column_names = FALSE, show_heatmap_legend = TRUE, row_names_gp = gpar(fontsize = 8))

	write.xlsx(data.frame(ht@matrix), paste0(result_dir, "8_Cell_scores/8.1_gene_exp_heatmap/", dataset_name, "_heatmap.xlsx"), rowNames=TRUE, overwrite=TRUE)

	pdf(paste0(result_dir, "8_Cell_scores/8.1_gene_exp_heatmap/", dataset_name, "_heatmap.pdf"))
	draw(ht)
	dev.off()

	jpeg(paste0(result_dir, "8_Cell_scores/8.1_gene_exp_heatmap/", dataset_name, "_heatmap.jpg"))
	print(ht)
	dev.off()

	# pptx <- paste0(result_dir, "8_Cell_scores/", dataset_name, "_heatmap.pptx")
	# heatmap <- draw(ht, padding=unit(c(50,10, 10, 20), "mm"))
	# topptx(heatmap, pptx)

	# 2) gene z-score Heatmap (without log2)
	exp_col_count <- ncol(exp_data_hit_nolog)
	exp_data_hit_nolog$mean <- apply(exp_data_hit_nolog[,1:exp_col_count], 1, mean)
	exp_data_hit_nolog$sd <- apply(exp_data_hit_nolog[,1:exp_col_count], 1, sd)

	gene_hit_z <- exp_data_hit_nolog[,1:exp_col_count]
	for(i in 1:nrow(exp_data_hit_nolog)){
		for(j in 1:exp_col_count){
			gene_hit_z[i,j] <- (exp_data_hit_nolog[i,j]-exp_data_hit_nolog[i,"mean"])/exp_data_hit_nolog[i,"sd"]
		}
	}

	exp_DEG <- as.matrix(gene_hit_z)
	ht <- Heatmap(exp_DEG, na_col="grey", cluster_rows = FALSE, cluster_columns = TRUE, show_row_dend = FALSE, show_row_names = TRUE, show_column_names = FALSE, show_heatmap_legend = TRUE, row_names_gp = gpar(fontsize = 8))

	write.xlsx(data.frame(ht@matrix), paste0(result_dir, "8_Cell_scores/8.2_gene_z-score_heatmap/", dataset_name, "_heatmap.xlsx"), rowNames=TRUE, overwrite=TRUE)

	pdf(paste0(result_dir, "8_Cell_scores/8.2_gene_z-score_heatmap/", dataset_name, "_heatmap.pdf"))
	draw(ht)
	dev.off()

	jpeg(paste0(result_dir, "8_Cell_scores/8.2_gene_z-score_heatmap/", dataset_name, "_heatmap.jpg"))
	print(ht)
	dev.off()

	# 3) ssGSEA Heatmap (with log2)
	library(GSVA)
	cell_gene_list <- lapply(as.list(top_genes_convert), function(x)x[!is.na(x)])
	ssgsea_score <- gsva(as.matrix(exp_data_hit_log), cell_gene_list, method="ssgsea")

	exp_DEG <- as.matrix(ssgsea_score)
	ht <- Heatmap(exp_DEG, na_col="grey", cluster_rows = FALSE, cluster_columns = TRUE, show_row_dend = FALSE, show_row_names = TRUE, show_column_names = FALSE, show_heatmap_legend = TRUE, row_names_gp = gpar(fontsize = 8)) # col=colorRamp2(c(1,0,-1), c("red", "white", "blue"))

	write.xlsx(data.frame(ssgsea_score), paste0(result_dir, "8_Cell_scores/8.3_ssGSEA_heatmap/", dataset_name, "_ssGSEA_score.xlsx"), rowNames=TRUE, overwrite=TRUE)

	pdf(paste0(result_dir, "8_Cell_scores/8.3_ssGSEA_heatmap/", dataset_name, "_ssGSEA_heatmap.pdf"))
	draw(ht)
	dev.off()

	jpeg(paste0(result_dir, "8_Cell_scores/8.3_ssGSEA_heatmap/", dataset_name, "_ssGSEA_heatmap.jpg"))
	print(ht)
	dev.off()

	# 4) Cell z-score Heatmap (without log2)
	cell_z <- data.frame(matrix(nrow=0, ncol=ncol(gene_hit_z)))
	names(cell_z) <- names(gene_hit_z)

	for(i in 1:ncol(top_genes_convert)){
		cell_gene_z <- na.omit(gene_hit_z[match(top_genes_convert[,i], rownames(gene_hit_z)),])
		cell_z[i,] <- apply(cell_gene_z, 2, sum)
		rownames(cell_z)[i] <- names(top_genes_convert)[i]
	}

	exp_DEG <- as.matrix(cell_z)
	ht <- Heatmap(exp_DEG, na_col="grey", cluster_rows = FALSE, cluster_columns = TRUE, show_row_dend = FALSE, show_row_names = TRUE, show_column_names = FALSE, show_heatmap_legend = TRUE, row_names_gp = gpar(fontsize = 8)) # col=colorRamp2(c(100,0,-100), c("red", "white", "blue"))

	write.xlsx(cell_z, paste0(result_dir, "8_Cell_scores/8.4_cell_z-score_heatmap/", dataset_name, "_cell_z-score.xlsx"), rowNames=TRUE, overwrite=TRUE)

	pdf(paste0(result_dir, "8_Cell_scores/8.4_cell_z-score_heatmap/", dataset_name, "_cell_z-score.pdf"))
	draw(ht)
	dev.off()

	jpeg(paste0(result_dir, "8_Cell_scores/8.4_cell_z-score_heatmap/", dataset_name, "_cell_z-score.jpg"))
	print(ht)
	dev.off()
}

######################################## main ########################################
# DEG for cell types
gene_survival <- read.xlsx(paste0(result_dir, "6_Gene_Survival/6.1_Gene_Survival.xlsx"))
cell_types <- unique(gene_survival$cell_type)

# original RNA expression matrix
datasets_all <- list.files(paste0(result_dir, "1_data_clean/"))
datasets_all <- datasets_all[grep(".txt$", datasets_all)]
datasets <- datasets_all[!grepl("^a|GSE88|GSE89|PMID", datasets_all)]

# 1) Find top genes for top cell types
top_genes <- gene_survival[gene_survival$KM_Pvalue < 0.05 & gene_survival$Forest_Pvalue < 0.05 & (gene_survival$HR_High < 0.5 | gene_survival$HR_High > 2.5) & gene_survival$HR_High != "not-converge", c("cell_type", "gene_name")]

top_genes_convert <- data.frame(matrix(nrow=0, ncol=0))
for(i in 1:length(cell_types)){
	cell_type <- cell_types[i]
	candidate_genes <- top_genes[top_genes$cell_type==cell_type,"gene_name"]
	top_genes_convert[1:length(candidate_genes),i] <- candidate_genes
	names(top_genes_convert)[i] <- cell_type
}

write.xlsx(top_genes, paste0(result_dir, "8_Cell_scores/Top_genes.xlsx"), rowNames=FALSE, overwrite=TRUE)
write.xlsx(top_genes_convert, paste0(result_dir, "8_Cell_scores/Top_genes_convert.xlsx"), rowNames=FALSE, overwrite=TRUE)

# 2) gene exp heatmap, gene z-score heatmap、ssGSEA heatmap、cell z-score heatmap VS survival
dir.create(paste0(result_dir, "8_Cell_scores/8.1_gene_exp_heatmap/"))
dir.create(paste0(result_dir, "8_Cell_scores/8.2_gene_z-score_heatmap/"))
dir.create(paste0(result_dir, "8_Cell_scores/8.3_ssGSEA_heatmap/"))
dir.create(paste0(result_dir, "8_Cell_scores/8.4_cell_z-score_heatmap/"))

# load normalized data: set (non-norm), set1 (seven housekeeping genes norm), set2 (empirical genes norm)
load(paste0(result_dir, "7_Housekeeping_normalization/7.99_RUVg_norm_process.Rdata"))

for(dataset in datasets){
	dataset_name <- gsub("_clean.txt", "", dataset)
	# Heatmap
	Cell_scores_heatmap(dataset_name, top_genes, top_genes_convert)
}
