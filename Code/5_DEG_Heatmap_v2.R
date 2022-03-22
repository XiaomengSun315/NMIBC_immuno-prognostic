# modified according to reviewer1's suggestion: use wFisher as combined p value method

rm(list=ls())
################################## BC progression ##################################
# install.packages("rJava")
# library(devtools)
# install_github("FedericoComoglio/rSymPy")
# install_github('unistbig/metapro', INSTALL_opts=c("--no-multiarch"))
library(openxlsx)
library(limma)
library(ggplot2)
library(ComplexHeatmap)
library(eoffice)
library(clusterProfiler)
library(dplyr)
library(rJava)
library(metapro)
library(rSymPy)
sympyStart()

# set working directory
data_dir <- "C:/0_xmsun/xmsun/Graduate/20210224_NMIBC/Data/"
result_dir <- "C:/0_xmsun/xmsun/Graduate/20210224_NMIBC/Results/"
dir.create(paste0(result_dir, "5_DEG_Heatmap_KEGG/"))

################################### subgroup by cells ###################################

groups <- function(cell_type, cell_alias, dataset_name, sources){

	cell_subgroup <- data.frame(matrix(nrow=0, ncol=2*length(sources)))
	names(cell_subgroup) <- c(paste0(sources, "_origin"), paste0(sources, "_level"))
	
	for(method_name in sources){

		# check if the method result exists in the dataset
		if(file.exists(paste0(result_dir, "2_immunedeconv/", dataset_name, "/", method_name, ".csv"))){

			# get immune infiltration score matrix
			score_matrix <- read.csv(paste0(result_dir, "2_immunedeconv/", dataset_name, "/", method_name, ".csv"), row.names=1, check.names=FALSE)
			hit_line <- grep(cell_alias, rownames(score_matrix))

			# name output sheets at first time
			if(length(rownames(cell_subgroup))==0){
				cell_subgroup[1:ncol(score_matrix),paste0(method_name, "_origin")] <- as.numeric(score_matrix[hit_line,])
				rownames(cell_subgroup) <- names(score_matrix)
			}else{
				# match scores and sample_names in following methods
				temp <- data.frame(t(score_matrix[hit_line,]))
				cell_subgroup[,paste0(method_name, "_origin")] <- as.numeric(temp[match(rownames(cell_subgroup), rownames(temp)),])
			}

			# subgroup by MEDIAN
			median <- median(as.numeric(score_matrix[hit_line,]))
			if(sum(summary((as.numeric(score_matrix[hit_line,]))))==0){
				# if MEDIAN is 0, delete corresponding subgroup
				delete_col <- grep(method_name, names(cell_subgroup))
				cell_subgroup <- cell_subgroup[,-delete_col]
			}else{
				cell_subgroup[cell_subgroup[,paste0(method_name, "_origin")] > median, paste0(method_name, "_level")] <- 1
				cell_subgroup[cell_subgroup[,paste0(method_name, "_origin")] <= median, paste0(method_name, "_level")] <- 0
			}
		}else{
			delete_col <- grep(method_name, names(cell_subgroup))
			cell_subgroup <- cell_subgroup[,-delete_col]
		}
	}

	# write cell_subgroup to cell_groups 
	current_line <- nrow(cell_groups) + 1
	if(current_line <= (current_line+nrow(cell_subgroup)-1)){
		cell_groups[current_line:(current_line+nrow(cell_subgroup)-1), "cell_type"] <<- cell_type
		cell_groups[current_line:(current_line+nrow(cell_subgroup)-1), "dataset_name"] <<- dataset_name
		cell_groups[current_line:(current_line+nrow(cell_subgroup)-1), "sample_name"] <<- rownames(cell_subgroup)
		for(col_name in names(cell_subgroup)){
			cell_groups[current_line:(current_line+nrow(cell_subgroup)-1), col_name] <<- cell_subgroup[,col_name]
		}
	}else{
		cell_groups[current_line, "cell_type"] <<- cell_type
		cell_groups[current_line, "dataset_name"] <<- dataset_name
		cell_groups[current_line, 3:ncol(cell_groups)] <<- ""
	}
	
	return(cell_subgroup)
}


######################################## DEG ########################################

DEG_func <- function(cell_type, dataset_name, cell_subgroup){

	# expression data
	exp_data <- read.table(paste0(result_dir, "1_data_clean/", dataset_name, "_clean.txt"), header=TRUE, row.names=1)

	# valid sources
	valid_sources <- unique(gsub("_[a-z]*$", "", names(cell_subgroup)))

	for(method_name in valid_sources){

		# 1) get DEG with limma (log2 expression as input data)
		design <- model.matrix( ~ cell_subgroup[,paste0(method_name, "_level")])
		rownames(design) <- rownames(cell_subgroup)
		colnames(design)[2] <- paste0(method_name, "_level")
		
		fit <- lmFit(exp_data, design=design)
		ebayes <- eBayes(fit, trend=TRUE)
		DEG_all <- topTable(ebayes, coef=2, n=Inf)
		write.xlsx(DEG_all, paste0(result_dir, "5_DEG_Heatmap_KEGG/", cell_type, "/", dataset_name, "_", method_name, "_DEG_all.xlsx"), rowNames=TRUE, overwrite=TRUE)

		DEG <- DEG_all[DEG_all$adj.P.Val < 0.05 & (DEG_all$logFC > 1 | DEG_all$logFC < -1),]
		write.xlsx(DEG_all, paste0(result_dir, "5_DEG_Heatmap_KEGG/", cell_type, "/", dataset_name, "_", method_name, "_DEG_all.xlsx"), rowNames=TRUE, overwrite=TRUE)
		write.xlsx(DEG, paste0(result_dir, "5_DEG_Heatmap_KEGG/", cell_type, "/", dataset_name, "_", method_name, "_DEG_filtered.xlsx"), rowNames=TRUE, overwrite=TRUE)

		# write DEG_all to DEGs
		current_line_DEGs <- nrow(DEGs) + 1
		if(nrow(DEG_all) > 0){
			DEGs[current_line_DEGs:(current_line_DEGs+nrow(DEG_all)-1), "cell_type"] <<- cell_type
			DEGs[current_line_DEGs:(current_line_DEGs+nrow(DEG_all)-1), "dataset_name"] <<- dataset_name
			DEGs[current_line_DEGs:(current_line_DEGs+nrow(DEG_all)-1), "method_name"] <<- method_name
			DEGs[current_line_DEGs:(current_line_DEGs+nrow(DEG_all)-1), "gene_name"] <<- rownames(DEG_all)
			DEGs[current_line_DEGs:(current_line_DEGs+nrow(DEG_all)-1), "logFC"] <<- DEG_all$logFC
			DEGs[current_line_DEGs:(current_line_DEGs+nrow(DEG_all)-1), "Pvalue"] <<- DEG_all$P.Value
			DEGs[current_line_DEGs:(current_line_DEGs+nrow(DEG_all)-1), "Padjust"] <<- DEG_all$adj.P.Val
		}else{
			DEGs[current_line_DEGs, "cell_type"] <<- cell_type
			DEGs[current_line_DEGs, "dataset_name"] <<- dataset_name
			DEGs[current_line_DEGs, "method_name"] <<- method_name
		}

		# 2) Heatmap
		column_df <- data.frame(cell_subgroup[,paste0(method_name, "_level")])
		rownames(column_df) <- rownames(cell_subgroup)
		colnames(column_df) <- "level"
		column_col <- list("level" = c("1" = "#1d91c0", "0" = "#f0f0f0"))
		column_ha <- HeatmapAnnotation(df = column_df, gp = gpar(col="white"), col = column_col, show_annotation_name = TRUE, show_legend = TRUE)

		exp_DEG <- as.matrix(exp_data[match(rownames(DEG), rownames(exp_data)),])
		ht <- Heatmap(exp_DEG, na_col="white", top_annotation = column_ha, cluster_rows = TRUE, cluster_columns = TRUE, show_row_dend = F, show_row_names = TRUE, show_column_names = F, show_heatmap_legend = TRUE, row_names_gp = gpar(fontsize = 8))

		write.xlsx(data.frame(ht@matrix), paste0(result_dir, "5_DEG_Heatmap_KEGG/", cell_type, "/", dataset_name, "_", method_name, "_heatmap.xlsx"), rowNames=TRUE, overwrite=TRUE)

		pdf(paste0(result_dir, "5_DEG_Heatmap_KEGG/", cell_type, "/", dataset_name, "_", method_name, "_heatmap.pdf"))
		draw(ht, padding=unit(c(50,10,10,20), "mm"))
		dev.off()

		jpeg(paste0(result_dir, "5_DEG_Heatmap_KEGG/", cell_type, "/", dataset_name, "_", method_name, "_heatmap.jpg"))
		print(ht)
		dev.off()

		# pptx <- paste0(result_dir, "5_DEG_Heatmap_KEGG/", cell_type, "/", dataset_name, "_", method_name, "_heatmap.pptx")
		# heatmap <- draw(ht, padding=unit(c(50,10, 10, 20), "mm"))
		# topptx(heatmap, pptx)
	}
}

################################ main: output format ################################
cells <- data.frame(matrix(nrow=11, ncol=3))
names(cells) <- c("cell_name", "sources", "cell_alias")
cells$cell_name <- c("Tcells_CD4+", "Tcells_CD8+", "DC", "Fibroblasts", "Macrophages", "Monocytes", "Tcells", "Endothelial", "Neutrophils", "Bcells", "NKcells")
cells$sources <- c("epic;quantiseq;timer", "epic;timer", "timer", "epic;mcp_counter", "mcp_counter;timer", "mcp_counter", "mcp_counter", "epic;mcp_counter", "mcp_counter;quantiseq;timer", "epic;mcp_counter;quantiseq;timer", "epic;mcp_counter;quantiseq")
cells$cell_alias <- c("T cell CD4+", "T cell CD8+", "Myeloid dendritic cell", "Cancer associated fibroblast", "Macrophage", "Monocyte", "T cell", "Endothelial cell", "Neutrophil", "B cell", "NK cell")

# output file 1: subgroups
cell_groups <- data.frame(matrix(nrow=0, ncol=11))
names(cell_groups) <- c("cell_type", "dataset_name", "sample_name", "epic_origin", "epic_level", "quantiseq_origin", "quantiseq_level", "mcp_counter_origin", "mcp_counter_level", "timer_origin", "timer_level")
# names(cell_groups) <- c("cell_type", "dataset_name", "method_name", "method_origin", "method_level")

# output file 2: DEG records
DEGs <- data.frame(matrix(nrow=0, ncol=7))
names(DEGs) <- c("cell_type", "dataset_name", "method_name", "gene_name", "logFC", "Pvalue", "Padjust")

# output file 3: top KEGG pathway enrichment 
KEGGs <- data.frame(matrix(nrow=0, ncol=7))
names(KEGGs) <- c("cell_type", "dataset_name", "method_name", "up_KEGG_order", "up_KEGG_name", "down_KEGG_order", "down_KEGG_name")


######################################## main ########################################
######################################## get DEG (limma)
for(cell_type in cells$cell_name){

	dir.create(paste0(result_dir, "5_DEG_Heatmap_KEGG/", cell_type))
	sources <- strsplit(cells[cells$cell_name==cell_type, "sources"], ";")[[1]]
	cell_alias <- cells[cells$cell_name==cell_type, "cell_alias"]
	datasets <- list.dirs(paste0(result_dir, "2_immunedeconv"), full.names=FALSE)[-1]
	# datasets <- datasets[-grep("GSE89", datasets)] # 其在计算DEG的logFC时容易过拟合。 # 后来GSE88、GSE89、PMID三个数据集直接删掉了。

	for(dataset_name in datasets){

		# DEG & Heatmap & KEGG enrichment
		cell_subgroup <- groups(cell_type, cell_alias, dataset_name, sources)
		DEG_func(cell_type, dataset_name, cell_subgroup)
	}
}

DEGs_filtered <- DEGs[DEGs$Padjust < 0.05 & (DEGs$logFC > 1 | DEGs$logFC < -1),]

write.xlsx(cell_groups, paste0(result_dir, "5_DEG_Heatmap_KEGG/5.1.1_cell_groups.xlsx"), rowNames=FALSE, overwrite=TRUE)
write.xlsx(DEGs, paste0(result_dir, "5_DEG_Heatmap_KEGG/5.1.2_DEGs_all.xlsx"), rowNames=FALSE, overwrite=TRUE)
write.xlsx(DEGs_filtered, paste0(result_dir, "5_DEG_Heatmap_KEGG/5.1.3_DEGs_filtered.xlsx"), rowNames=FALSE, overwrite=TRUE)
save(cell_groups, DEGs, DEGs_filtered, file=paste0(result_dir, "5_DEG_Heatmap_KEGG/5.1.4_DEG_process.Rdata"))

######################################## Select genes for each cell type (mean, ordmeta, wFisher)
# sample size in each dataset
sample_size <- data.frame(dataset=datasets, sample_size=c(86,460,84,159,194,47,213,69,58))

# output format
DEG_summary <- data.frame(matrix(nrow=0, ncol=12))
names(DEG_summary) <- c("Cell_type", "Gene_name", "Dataset_Count", "mean_logFC", "mean_Pvalue", "mean_Padjust", "ordmeta_combined_p", "ordmeta_optimal_rank", "ordmeta_minimum_marginal_p", "ordmeta_overall_eff_direction", "wFisher_combined_p", "wFisher_overall_eff_direction")

cell_types <- unique(DEGs$cell_type)
for(cell_type in cell_types){
	DEGs_cell <- DEGs[DEGs$cell_type==cell_type,]
	genes <- na.omit(unique(DEGs_cell$gene_name))
	current_line_DEG_summary <- nrow(DEG_summary) + 1
	DEG_summary[current_line_DEG_summary:(current_line_DEG_summary+length(genes)-1),"Cell_type"] <- cell_type
	for(i in 1:length(genes)){
		# remark
		DEG_summary[current_line_DEG_summary+i-1,"Gene_name"] <- genes[i]
		DEGs_gene <- na.omit(DEGs_cell[DEGs_cell$gene_name==genes[i],])
		DEG_summary[current_line_DEG_summary+i-1,"Dataset_Count"] <- length(unique(DEGs_gene$dataset_name))

		# mean
		DEG_summary[current_line_DEG_summary+i-1,"mean_logFC"] <- mean(DEGs_gene$logFC)
		DEG_summary[current_line_DEG_summary+i-1,"mean_Pvalue"] <- mean(DEGs_gene$Pvalue)
		DEG_summary[current_line_DEG_summary+i-1,"mean_Padjust"] <- mean(DEGs_gene$Padjust)

		# direction
		DEGs_gene$eff_direction <- ""
		DEGs_gene[DEGs_gene$logFC > 0,"eff_direction"] <- 1
		DEGs_gene[DEGs_gene$logFC < 0,"eff_direction"] <- -1

		# # ordmeta
		# ordmeta_p <- ordmeta(p=DEGs_gene$Pvalue, is.onetail=FALSE, eff.sign=DEGs_gene$eff_direction)

		# DEG_summary[current_line_DEG_summary+i-1,"ordmeta_combined_p"] <- ordmeta_p$p
		# DEG_summary[current_line_DEG_summary+i-1,"ordmeta_optimal_rank"] <- ordmeta_p$optimal_rank
		# DEG_summary[current_line_DEG_summary+i-1,"ordmeta_minimum_marginal_p"] <- ordmeta_p$MMP
		# DEG_summary[current_line_DEG_summary+i-1,"ordmeta_overall_eff_direction"] <- ordmeta_p$overall.eff.direction

		# wFisher
		DEGs_gene_size <- merge(DEGs_gene, sample_size, by.x="dataset_name", by.y="dataset", sort=FALSE)

		wfisher_p <- wFisher(p=DEGs_gene_size$Pvalue, weight=DEGs_gene_size$sample_size, is.onetail=FALSE, eff.sign=DEGs_gene_size$eff_direction)
		DEG_summary[current_line_DEG_summary+i-1,"wFisher_combined_p"] <- wfisher_p$p
		DEG_summary[current_line_DEG_summary+i-1,"wFisher_overall_eff_direction"] <- wfisher_p$overall.eff.direction
	}

	DEG_summary[current_line_DEG_summary:(current_line_DEG_summary+length(genes)-1),"wFisher_sig"] <- ""
	sig_hit <- intersect(rownames(DEG_summary[DEG_summary[current_line_DEG_summary:(current_line_DEG_summary+length(genes)-1),"wFisher_combined_p"] < 0.05/length(genes),]), current_line_DEG_summary:(current_line_DEG_summary+length(genes)-1))
	not_sig <- setdiff(current_line_DEG_summary:(current_line_DEG_summary+length(genes)-1), sig_hit)
	DEG_summary[sig_hit,"wFisher_sig"] <- "TRUE"
	DEG_summary[not_sig,"wFisher_sig"] <- "FALSE"
}

# DEG_summary <- arrange(DEG_summary, Cell_type, -Dataset_Count, -mean_logFC, mean_Padjust)
write.xlsx(DEG_summary, paste0(result_dir, "5_DEG_Heatmap_KEGG/5.1.5_DEG_summary.xlsx"), rowNames=FALSE, overwrite=TRUE)

DEG_summary_filtered <- DEG_summary[DEG_summary$Dataset_Count >= 3 & DEG_summary$wFisher_sig == "TRUE",]
write.xlsx(DEG_summary_filtered, paste0(result_dir, "5_DEG_Heatmap_KEGG/5.1.6_DEG_summary_filtered.xlsx"), rowNames=FALSE, overwrite=TRUE)
save(DEG_summary, DEG_summary_filtered, file=paste0(result_dir, "5_DEG_Heatmap_KEGG/5.1.7_DEG_summary_process.Rdata"))
