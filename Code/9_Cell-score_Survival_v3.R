# modified according to reviewer1's suggestion: update with "8_Cell_scores_v3.R" or "8_Cell_scores_v4.R". 

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
data_dir <- "C:/0_xmsun/xmsun/Graduate/20210224_NMIBC/Data/"
result_dir <- "C:/0_xmsun/xmsun/Graduate/20210224_NMIBC/Results/"
dir.create(paste0(result_dir, "9_Cell_Survival/"))

################################### Function module ################################## 

######### Function: Survival Analysis
survival_pack <- function(dataset_name, survival_type, survival_time){
	
	pheno_matrix <- pheno[pheno$Source_dataset==dataset_name,]

	for(method in c("ssGSEA", "Z-score")){
		if(method == "ssGSEA" && file.exists(paste0(result_dir, "8_Cell_scores/8.3_ssGSEA_heatmap/", dataset_name, "_ssGSEA_score.xlsx"))){
			method_name <- method
			score_matrix <- read.xlsx(paste0(result_dir, "8_Cell_scores/8.3_ssGSEA_heatmap/", dataset_name, "_ssGSEA_score.xlsx"), rowNames=TRUE)
			survival_plot(dataset_name, method_name, score_matrix, pheno_matrix, survival_type, survival_time)
		}else if(method == "Z-score" && file.exists(paste0(result_dir, "8_Cell_scores/8.4_cell_z-score_heatmap/", dataset_name, "_cell_z-score.xlsx"))){
			method_name <- method
			score_matrix <- read.xlsx(paste0(result_dir, "8_Cell_scores/8.4_cell_z-score_heatmap/", dataset_name, "_cell_z-score.xlsx"), rowNames=TRUE)
			survival_plot(dataset_name, method_name, score_matrix, pheno_matrix, survival_type, survival_time)
		}
	}
}

################################### Function: median ###################################
survival_plot <- function(dataset_name, method_name, score_matrix, pheno_matrix, survival_type, survival_time){

	dir.create(paste0(result_dir, "9_Cell_Survival/", survival_time))
	current_line <- nrow(output)

	merged_data_all <- data.frame(matrix(nrow=0, ncol=3))
	names(merged_data_all) <- c("Sample_name", "survival_type", "survival_time")
	
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
		if(i == 1){
			current_col <- ncol(merged_data_all)
			merged_data_all[1:nrow(merged_data),"Sample_name"] <- merged_data$Sample_name
			merged_data_all$survival_type <- merged_data$survival_type
			merged_data_all$survival_time <- merged_data$survival_time
			merged_data_all$cell_type <- merged_data[,cell_type]
			merged_data_all$score_class <- merged_data$score_class
			names(merged_data_all)[(current_col+1):(current_col+2)] <- c(paste0(cell_type, "_score"), paste0(cell_type, "_level"))
		}else{
			current_col <- ncol(merged_data_all)
			merged_data_all[,(current_col+1)] <- merged_data[match(merged_data_all$Sample_name, merged_data$Sample_name),cell_type]
			merged_data_all[,(current_col+2)] <- merged_data[match(merged_data_all$Sample_name, merged_data$Sample_name),"score_class"]
			names(merged_data_all)[(current_col+1):(current_col+2)] <- c(paste0(cell_type, "_score"), paste0(cell_type, "_level"))
		}

		######## 1) K-M curve
		fit <- survfit(Surv(survival_time, survival_type) ~ score_class, data=merged_data)
		surv_plot <- ggsurvplot(fit, data=merged_data, suv.median.line="hv", conf.int=TRUE, risk.table=TRUE, pval=TRUE)
		cum_plot <- ggsurvplot(fit, data=merged_data, conf.int=TRUE, fun="cumhaz")

		if(grepl("/", cell_type)){
			cell_type <- gsub("/", " ", cell_type)
		}
		jpeg(paste0(result_dir, "9_Cell_Survival/", survival_time, "/", dataset_name, "_", method_name, "_", cell_type, "_surv_plot.jpg"))
		print(surv_plot)
		dev.off()

		jpeg(paste0(result_dir, "9_Cell_Survival/", survival_time, "/", dataset_name, "_", method_name, "_", cell_type, "_cum_plot.jpg"))
		print(cum_plot)
		dev.off()

		output[current_line+i, "survival_type"] <<- survival_type
		output[current_line+i, "survival_time"] <<- survival_time
		output[current_line+i, "dataset_name"] <<- dataset_name
		output[current_line+i, "method_name"] <<- method_name
		output[current_line+i, "cell_type"] <<- cell_type
		output[current_line+i, "KM_Pvalue"] <<- surv_pvalue(fit, merged_data)$pval

		######## 2) Forest plot
		if(length(unique(merged_data$score_class))>1){
			merged_data <- within(merged_data, {score_class <- factor(score_class, labels=c("High", "Low"))})
		}else{
			merged_data <- within(merged_data, {score_class <- factor(score_class, labels="Low")})
		}

		fit2 <- tryCatch({coxph(Surv(survival_time, survival_type) ~ score_class , data=merged_data)}, warning = function(w){print("warning")}, error=function(e){print("error")})
		if(fit2 == "error"){
			output[current_line+i, "HR_High"] <<- "single_outcome"
			output[current_line+i, "Forest_Pvalue"] <<- "single_outcome"
		}else if(fit2 == "warning"){
			output[current_line+i, "HR_High"] <<- "not-converge"
			output[current_line+i, "Forest_Pvalue"] <<- "not-converge"
		}else{
			forest_plot <- ggforest(model = fit2, data = merged_data)

			jpeg(paste0(result_dir, "9_Cell_Survival/", survival_time, "/", dataset_name, "_", method_name, "_", cell_type, "_forest_plot.jpg"))
			print(forest_plot)
			dev.off()

			output[current_line+i, "HR_High"] <<- summary(fit2)$conf.int[1]
			output[current_line+i, "Forest_Pvalue"] <<- summary(fit2)$coefficients[5]
		}
	}
	write.xlsx(merged_data_all, paste0(result_dir, "9_Cell_Survival/", survival_time, "/", dataset_name, "_", method_name, "_merged_data_all.xlsx"), overwrite=TRUE)
}

######################################## main ########################################
pheno <- read.xlsx(paste0(data_dir, "combined/Clinicopathologic_v2.xlsx"), na.strings=c("","NA","Nan","#N/A"))
datasets <- c("E-MTAB-4321", "GSE32894", "GSE13507")

# output format
output <- data.frame(matrix(nrow=0, ncol=8))
names(output) <- c("survival_type", "survival_time", "dataset_name", "method_name", "cell_type", "KM_Pvalue", "HR_High", "Forest_Pvalue")
current_line <- nrow(output)

# survival plot 
# function format: dataset_name, survival_type, survival_time)
survival_pack("E-MTAB-4321", "Progression_beyond_T2_Progression", "Progression_free_survival")
survival_pack("GSE32894", "Cancer_specific_satus_DOD_event", "Disease_specific_survival")
survival_pack("GSE13507", "Vital_status", "OS")

# output results
write.xlsx(output, paste0(result_dir, "9_Cell_Survival/Cell_Survival.xlsx"), overwrite=TRUE)
