# modified according to reviewer1's suggestion: update with "5_DEG_Heatmap_v2.R"

rm(list=ls())
################################## BC progression ##################################
library(dplyr)
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)

# set working directory
args <- commandArgs(T)
data_dir <- paste0(args[1], "Data/")
result_dir <- paste0(args[1], "Results/")
dir.create(paste0(result_dir, "6_Gene_Survival/"))


################################### Function module ################################## 

######### Function: Survival Analysis
survival_pack_gene <- function(dataset_name, survival_type, survival_time, cell_type, candidate_genes){
	
	pheno_matrix <- pheno[pheno$Source_dataset==dataset_name,]
	exp_data <- read.table(paste0(result_dir, "1_data_clean/", dataset_name, "_clean.txt"), header=TRUE, row.names=1)
	exp_matrix <- na.omit(exp_data[match(candidate_genes, rownames(exp_data)),])

	output <- survival_plot(dataset_name, exp_matrix, pheno_matrix, survival_type, survival_time, cell_type, output)
	return(output)
}

################################## Function: median ##################################
survival_plot <- function(dataset_name, exp_matrix, pheno_matrix, survival_type, survival_time, cell_type, output){

	dir.create(paste0(result_dir, "6_Gene_Survival/", cell_type))
	current_line <- nrow(output)

	merged_data_all <- data.frame(matrix(nrow=0, ncol=3))
	names(merged_data_all) <- c("Sample_name", "survival_type", "survival_time")
	
	# survival data
	survival_data <- pheno_matrix[,c("Sample_name", survival_type, survival_time)]

	for(i in 1:nrow(exp_matrix)){
		
		gene_name <- rownames(exp_matrix)[i]
		print(paste("Start session:", survival_type, "-", survival_time, "-", dataset_name, "-", cell_type, "-", gene_name))

		score_data <- data.frame(t(exp_matrix[i,]), check.names=FALSE)
		score_data$gene_class <- ""
		score_data$gene_class[score_data[,gene_name] > median(score_data[,gene_name])] <- "High"
		score_data$gene_class[score_data[,gene_name] <= median(score_data[,gene_name])] <- "Low"
		score_data$Sample_name <- row.names(score_data)

		merged_data <- data.frame(merge(score_data, survival_data, by="Sample_name"))
		names(merged_data) <- c("Sample_name", gene_name, "gene_class", "survival_type", "survival_time")
		if(i == 1){
			current_col <- ncol(merged_data_all)
			merged_data_all[1:nrow(merged_data),"Sample_name"] <- merged_data$Sample_name
			merged_data_all$survival_type <- merged_data$survival_type
			merged_data_all$survival_time <- merged_data$survival_time
			merged_data_all$gene_score <- merged_data[,gene_name]
			merged_data_all$gene_class <- merged_data$gene_class
			names(merged_data_all)[(current_col+1):(current_col+2)] <- c(paste0(gene_name, "_score"), paste0(gene_name, "_level"))
		}else{
			current_col <- ncol(merged_data_all)
			merged_data_all[,(current_col+1)] <- merged_data[match(merged_data_all$Sample_name, merged_data$Sample_name),gene_name]
			merged_data_all[,(current_col+2)] <- merged_data[match(merged_data_all$Sample_name, merged_data$Sample_name),"gene_class"]
			names(merged_data_all)[(current_col+1):(current_col+2)] <- c(paste0(cell_type, "_score"), paste0(gene_name, "_level"))
		}

		######## 1) K-M curve
		fit <- survfit(Surv(survival_time, survival_type) ~ gene_class, data=merged_data)
		surv_plot <- ggsurvplot(fit, data=merged_data, suv.median.line="hv", conf.int=TRUE, risk.table=TRUE, pval=TRUE)
		cum_plot <- ggsurvplot(fit, data=merged_data, conf.int=TRUE, fun="cumhaz")

		if(grepl("/", gene_name)){
			gene_name <- gsub("/", " ", gene_name)
		}
		jpeg(paste0(result_dir, "6_Gene_Survival/", cell_type, "/", dataset_name, "_", gene_name, "_surv_plot.jpg"))
		print(surv_plot)
		dev.off()

		jpeg(paste0(result_dir, "6_Gene_Survival/", cell_type, "/", dataset_name, "_", gene_name, "_cum_plot.jpg"))
		print(cum_plot)
		dev.off()

		output[current_line+i, "survival_type"] <- survival_type
		output[current_line+i, "survival_time"] <- survival_time
		output[current_line+i, "dataset_name"] <- dataset_name
		output[current_line+i, "cell_type"] <- cell_type
		output[current_line+i, "gene_name"] <- gene_name
		output[current_line+i, "KM_Pvalue"] <- surv_pvalue(fit, merged_data)$pval

		######## 2) Single-varint: Cox analysis & Forest plot
		merged_data <- within(merged_data, {gene_class <- factor(gene_class, labels=c("Low", "Hign"))})
		fit2 <- tryCatch({coxph(Surv(survival_time, survival_type) ~ gene_class , data=merged_data)}, warning = function(w){print("warning")})
		if(fit2 != "warning"){
			forest_plot <- ggforest(model = fit2, data = merged_data)
			
			jpeg(paste0(result_dir, "6_Gene_Survival/", cell_type, "/", dataset_name, "_", gene_name, "_forest_plot.jpg"))
			print(forest_plot)
			dev.off()

			output[current_line+i, "HR_High"] <- summary(fit2)$conf.int[1]
			output[current_line+i, "Forest_Pvalue"] <- summary(fit2)$coefficients[5]
		}else{
			output[current_line+i, "HR_High"] <- "not-converge"
			output[current_line+i, "Forest_Pvalue"] <- ""
		}

		######## 3) Pearson Correlation: Gene score & Survival time
		# write.xlsx(merged_data, paste0(result_dir, "6_Gene_Survival/", cell_type, "/", dataset_name, "_", gene_name, "_merged_data.xlsx"), overwrite=TRUE)
	}

	write.xlsx(merged_data_all, paste0(result_dir, "6_Gene_Survival/", cell_type, "/", dataset_name, "_merged_data_all.xlsx"), overwrite=TRUE)

	return(output)
}


######################################## main ########################################
pheno <- read.xlsx(paste0(data_dir, "combined/Clinicopathologic_v2.xlsx"), na.strings=c("","NA","Nan","#N/A"))
datasets <- c("E-MTAB-4321", "GSE32894", "GSE13507")

# output format
output <- data.frame(matrix(nrow=0, ncol=8))
names(output) <- c("survival_type", "survival_time", "dataset_name", "cell_type", "gene_name", "KM_Pvalue", "HR_High", "Forest_Pvalue")
current_line <- nrow(output)

# Candidate genes for each cell type
# load data: DEG_summary, DEG_summary_filtered
load(paste0(result_dir, "5_DEG_Heatmap_KEGG/5.1.7_DEG_summary_process.Rdata"))
cell_types <- unique(DEG_summary_filtered$Cell_type)
cell_types <- cell_types[!grepl("^Tcells$", cell_types, perl=TRUE)]
for(cell_type in cell_types){
	if(cell_type == "NKcells"){
		candidate_genes <- DEG_summary_filtered[(DEG_summary_filtered$mean_logFC > 0.3 | DEG_summary_filtered$mean_logFC < -0.3) & DEG_summary_filtered$Cell_type == cell_type, "Gene_name"]
	}else{
		candidate_genes <- DEG_summary_filtered[(DEG_summary_filtered$mean_logFC > 0.5 | DEG_summary_filtered$mean_logFC < -0.5) & DEG_summary_filtered$Cell_type == cell_type, "Gene_name"]
	}

	# survival plot (function format: dataset_name, survival_type, survival_time)
	output <- survival_pack_gene("E-MTAB-4321", "Progression_beyond_T2_Progression", "Progression_free_survival", cell_type, candidate_genes)
	output <- survival_pack_gene("GSE32894", "Cancer_specific_satus_DOD_event", "Disease_specific_survival", cell_type, candidate_genes)
	output <- survival_pack_gene("GSE13507", "Vital_status", "OS", cell_type, candidate_genes)
}

# output results
write.xlsx(output, paste0(result_dir, "6_Gene_Survival/6.1_Gene_Survival.xlsx"), overwrite=TRUE)
