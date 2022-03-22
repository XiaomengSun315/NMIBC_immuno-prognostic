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
dir.create(paste0(result_dir, "4_Survival/"))

################################### Function module ################################## 

######### Function: Survival Analysis
survival_pack <- function(dataset_name, survival_type, survival_time){
	
	pheno_matrix <- pheno[pheno$Source_dataset==dataset_name,]

	for(method in c("xCell", "ESTIMATE", "immunedeconv")){
		if(method == "xCell" && file.exists(paste0(result_dir, "2_", method, "/", dataset_name, "_heatmap.csv"))){
			method_name <- method
			score_matrix <- read.csv(paste0(result_dir, "2_", method, "/", dataset_name, "_heatmap.csv"), row.names=1, check.names=FALSE)
			output <- survival_plot(dataset_name, method_name, score_matrix, pheno_matrix, survival_type, survival_time, output)
		}else if(method == "ESTIMATE" && file.exists(paste0(result_dir, "2_", method, "/", dataset_name, "_estimate_score.xlsx"))){
			method_name <- method
			score_matrix <- read.xlsx(paste0(result_dir, "2_", method, "/", dataset_name, "_estimate_score.xlsx"), rowNames=TRUE)
			output <- survival_plot(dataset_name, method_name, score_matrix, pheno_matrix, survival_type, survival_time, output)
		}else if(method == "immunedeconv"){
			sub_methods <- list.files(paste0(result_dir, "2_", method, "/", dataset_name))
			sub_methods <- gsub(".csv", "", sub_methods[grep(".csv", sub_methods)])
			for(sub_method in sub_methods){
				method_name <- paste0(method, "_", sub_method)
				score_matrix <- read.csv(paste0(result_dir, "2_", method, "/", dataset_name, "/", sub_method, ".csv"), row.names=1, check.names=FALSE)
				output <- survival_plot(dataset_name, method_name, score_matrix, pheno_matrix, survival_type, survival_time, output)
			}
		}
	}
	return(output)
}

######### Function: median
survival_plot <- function(dataset_name, method_name, score_matrix, pheno_matrix, survival_type, survival_time, output){

	dir.create(paste0(result_dir, "4_Survival/", survival_time))
	current_line <- nrow(output)

	# survival data
	survival_data <- pheno_matrix[,c("Sample_name", survival_type, survival_time)]

	for(i in 1:nrow(score_matrix)){
		
		cell_type <- rownames(score_matrix)[i]
		print(paste("Start session:", survival_type, "-", survival_time, "-", dataset_name, "-", method_name, "-", cell_type))

		score_data <- data.frame(t(score_matrix[i,]), check.names=FALSE)
		score_data$score_class <- ""
		score_data$score_class[score_data[,cell_type] > median(score_data[,cell_type])] <- "High"
		score_data$score_class[score_data[,cell_type] <= median(score_data[,cell_type])] <- "Low"
		score_data$Sample_name <- row.names(score_data)

		merged_data <- data.frame(merge(score_data, survival_data, by="Sample_name"))
		names(merged_data) <- c("Sample_name", cell_type, "score_class", "survival_type", "survival_time")
		fit <- survfit(Surv(survival_time, survival_type) ~ score_class, data=merged_data)
		surv_plot <- ggsurvplot(fit, data=merged_data, suv.median.line="hv", conf.int=TRUE, risk.table=TRUE, pval=TRUE)
		cum_plot <- ggsurvplot(fit, data=merged_data, conf.int=TRUE, fun="cumhaz")

		if(grepl("/", cell_type)){
			cell_type <- gsub("/", " ", cell_type)
		}
		jpeg(paste0(result_dir, "4_Survival/", survival_time, "/", dataset_name, "_", method_name, "_", cell_type, "_surv_plot.jpg"))
		print(surv_plot)
		dev.off()

		jpeg(paste0(result_dir, "4_Survival/", survival_time, "/", dataset_name, "_", method_name, "_", cell_type, "_cum_plot.jpg"))
		print(cum_plot)
		dev.off()

		output[current_line+i, "survival_type"] <- survival_type
		output[current_line+i, "survival_time"] <- survival_time
		output[current_line+i, "dataset_name"] <- dataset_name
		output[current_line+i, "method_name"] <- method_name
		output[current_line+i, "cell_type"] <- cell_type
		output[current_line+i, "Pvalue"] <- surv_pvalue(fit, merged_data)$pval
	}

	return(output)
}



######################################## main ########################################
pheno <- read.xlsx(paste0(data_dir, "combined/Clinicopathologic_v2.xlsx"), na.strings=c("","NA","Nan","#N/A"))
datasets <- c("E-MTAB-4321", "GSE32894", "GSE13507")

# output format
output <- data.frame(matrix(nrow=0, ncol=6))
names(output) <- c("survival_type", "survival_time", "dataset_name", "method_name", "cell_type", "Pvalue")
current_line <- nrow(output)

# survival plot
output <- survival_pack("E-MTAB-4321", "Progression_beyond_T2_Progression", "Progression_free_survival")
output <- survival_pack("GSE32894", "Cancer_specific_satus_DOD_event", "Disease_specific_survival")
output <- survival_pack("GSE13507", "Vital_status", "OS")

# output results
write.xlsx(output, paste0(result_dir, "4_Survival/survival.xlsx"), overwrite=TRUE)
