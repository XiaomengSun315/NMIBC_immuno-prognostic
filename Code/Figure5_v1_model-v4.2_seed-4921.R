# Figure 5: v4.2, seed 4921

rm(list=ls())
################################## BC progression ##################################
# load package
# install.packages("openxlsx", dependencies=TRUE)
# install.packages("ggsignif", dependencies=TRUE)
library(dplyr)
library(readxl)
library(openxlsx)
library(eoffice)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(glmnet)
library(pROC)
library(doBy)
library(nnet)
library(car)
library(questionr)
library(multiROC)
library(dummies)
library(plyr)
library(ggsci)

# set working directory
args <- commandArgs(T)
data_dir <- paste0(args[1], "Data/")
result_dir <- paste0(args[1], "Results/")

dir.create(paste0(result_dir, "10_Models_plot/"))

################################### Function module ##################################
######### Function: Box plot
box_plot <- function(current_info, merged_data, path, survival_type, method_name, random_seed){
	for(group_name in c("Sex", "Age", "Staging", "Grade73", "Grade98", "EORTC_risk_score_NMIBC", "Tumor_size", "Tumor_status", "CIS_in_disease_course", "Recurrence", "Progression_beyond_T2_Progression", "Cancer_specific_satus_DOD_event", "Vital_status", "Survival_risk_NMIBConly", "Survival_risk_combine_NMIBConly")){
		# Results directory for bar plots
		setwd(path)
		print(paste("Start BOX session:", survival_type, "-", method_name, "-", random_seed, "-", group_name))

		# box data
		box_data <- merged_data%>%filter(is.na(get(group_name))==FALSE, get(group_name)!="")
		if(nrow(box_data)!=0){

			rownames(box_data) <- box_data$Sample_name
			
			# relevel factors
			box_data$Sex <- factor(box_data$Sex, levels=c("Female", "Male"))
			box_data$Age <- factor(box_data$Age, levels=c("<=30", "31-50", "51-60", "61-70", "71-80", ">80"))
			box_data$Staging <- factor(box_data$Staging, levels=c("T0", "Ta", "Ta-1", "T1", "Tis"))
			box_data$Grade73 <- factor(box_data$Grade73, levels=c("G0", "G1", "G2", "G3", "Gx"))
			box_data$Grade98 <- factor(box_data$Grade98, levels=c("Low", "High"))
			box_data$EORTC_risk_score_NMIBC <- factor(box_data$EORTC_risk_score_NMIBC, levels=c("0", "1"))
			box_data$Tumor_size <- factor(box_data$Tumor_size, levels=c("<=3cm", ">3cm"))
			box_data$Tumor_status <- factor(box_data$Tumor_status, levels=c("Primary", "Recurrent", "Progress"))
			box_data$CIS_in_disease_course <- factor(box_data$CIS_in_disease_course, levels=c("CIS-", "CIS+"))
			box_data$Recurrence <- factor(box_data$Recurrence, levels=c("FALSE", "TRUE"))
			box_data$Progression_beyond_T2_Progression <- factor(box_data$Progression_beyond_T2_Progression, levels=c("FALSE", "TRUE"))
			box_data$Cancer_specific_satus_DOD_event <- factor(box_data$Cancer_specific_satus_DOD_event, levels=c("FALSE", "TRUE"))
			box_data$Vital_status <- factor(box_data$Vital_status, levels=c("FALSE", "TRUE"))
			box_data$Survival_risk_NMIBConly <- factor(box_data$Survival_risk_NMIBConly, levels=c("A", "B", "C", "D", "E"))
			box_data$Survival_risk_combine_NMIBConly <- factor(box_data$Survival_risk_combine_NMIBConly, levels=c("Low", "High"))

			# col_name <- names(box_data)[i]
			temp_data <- box_data[,c("pred_num", group_name)]
			names(temp_data) <- c("Scores", "Group")
			temp_data <- temp_data %>% group_by(Group) %>% filter(n() >= 3) # delete groups with less than 3 samples

			levels <- length(unique(temp_data$Group))
			if(levels > 1){
				compare <- compare_means(Scores ~  Group, data=temp_data, method="t.test")
				my_comparisons <- split(unlist(sapply(data.frame(compare[,c("group1", "group2")]),c)), seq(nrow(compare)))
	
				box <- ggplot(temp_data, aes(Group, Scores, fill=Group)) +
					geom_boxplot(width=0.5) +
					stat_compare_means(comparisons=my_comparisons, method="t.test", size=8) +
					# geom_jitter(aes(fill=Group), width=0.2, shape=21, size=1) +
					labs(x=group_name) +
					scale_fill_aaas() +
					theme_bw() +
					theme(plot.title=element_text(size = 25), 
						# plot.title.position="top",
						axis.text.x=element_text(size=15, angle=0), 
						axis.text.y=element_text(size=15), 
						axis.title.x=element_text(size = 23, vjust=-1), 
						axis.title.y=element_blank(), 
						panel.grid.major=element_blank(),
						panel.grid.minor = element_blank(),
						panel.background = element_blank(),
						axis.line = element_line(colour = "black"),
						legend.position="top",
						legend.key.size=unit(1,"cm"),
						legend.text = element_text(size=13),
						legend.title = element_text(size=15))

				pdf(paste0(path, "clinical_box/", group_name, "_Box.pdf"))
				print(box)
				dev.off()

				jpeg(paste0(path, "clinical_box/", group_name, "_Box.jpg"))
				print(box)
				dev.off()

				topptx(box, paste0(path, "clinical_box/", group_name, "_Box.pptx"))
			}
		}
	}
}

######################################## main ########################################
# phenotype
pheno <- read.xlsx(paste0(data_dir, "combined/Clinicopathologic_v2.xlsx"), na.strings=c("","NA","Nan","#N/A"))
pheno <- pheno[!is.na(pheno$Survival_risk_NMIBConly),]
pheno <- pheno[pheno$Survival_risk_NMIBConly!="ref",]

# candidata models from v4.2 & v4.3 & v4.4
candidate_models <- read.xlsx(paste0(result_dir, "10_Models_summary.xlsx"))

# output format
output <- data.frame(matrix(nrow=0, ncol=15))
names(output)[1:ncol(candidate_models)] <- names(candidate_models)
names(output)[(ncol(candidate_models)+1):15] <- c("plot_type", "KM_survival_plot_type", "KM_survival_plot_time", "KM_Pvalue")
# for(version in unique(candidate_models$version)){
for(version in c("v4.2")){

	version_models <- candidate_models[candidate_models$version==version,]
	
	for(random_seed in c("4921")){

		current_info <- version_models[version_models$random_seed==random_seed,]
		survival_type <- current_info[,"survival_type"]
		dir.create(paste0(result_dir, "10_Models_plot/", version, "_Figure5/", survival_type, "/seed_", random_seed, "/"), recursive=TRUE)
		path <- paste0(result_dir, "10_Models_plot/", version, "_Figure5/", survival_type, "/seed_", random_seed, "/")
		if(version!="v4.4"){
			dir.create(paste0(path, "clinical_box/"), recursive=TRUE)
		}

		################################# 0) load model
		pheno[,eval(survival_type)] <- factor(pheno[,eval(survival_type)])
		levels <- levels(pheno[,eval(survival_type)])
		method_name <- version_models[version_models$random_seed==random_seed,"method_name"]
		load(file=paste0(result_dir, "10_Models_", version, "/10.2_", survival_type, "_", method_name, "_random/seed_", random_seed, "/1_model.Rdata"))

		################################# 1) prepare data
		# score matrix of cells (only ssGSEA methods)
		score_matrix <- data.frame(t(read.xlsx(paste0(result_dir, "8_Cell_scores/8.3_", method_name, "_matrix.xlsx"), rowNames=TRUE)), check.names=FALSE)
		score_col <- ncol(score_matrix)
		if(version %in% c("v4.2")){
			names(score_matrix) <- paste0(names(score_matrix), "_score")
		}

		# predicted outcome by model
		score_matrix[,(score_col+1):(score_col+length(levels))] <- predict(model, newdata=score_matrix, type="prob", se=TRUE)
		names(score_matrix)[(score_col+1):(score_col+length(levels))] <- paste0(levels, "_pred_MN")
		score_matrix$pred <- predict(model, newdata=score_matrix)

		# numeric predictions by model
		if(version %in% c("v4.2", "v4.3")){
			
			score_matrix$pred_num <- ""
			intercept <- as.numeric(current_info[,"intercept"])
			markers <- gsub("`", "", strsplit(current_info[,"markers"], split="|", fixed=TRUE)[[1]])
			coefficients <- as.numeric(strsplit(current_info[,"coefficients"], split="|", fixed=TRUE)[[1]])

			for(row in 1:nrow(score_matrix)){
				markers_score <- as.numeric(score_matrix[row, match(markers, names(score_matrix))])
				sum <- intercept
				for(marker_count in 1:length(markers)){
					sum <- sum + coefficients[marker_count] * markers_score[marker_count]
				}
				score_matrix[row, "pred_num"] <- sum
			}
			score_matrix$pred_num <- as.numeric(score_matrix$pred_num)
		}else{
			coefficients <- data.frame(coef(model), check.names=FALSE)
			names(coefficients) <- gsub("`|\\(|\\)", "", names(coefficients))
			markers <- names(coefficients)[-1]

			for(row in 1:nrow(score_matrix)){
				markers_score <-  as.numeric(score_matrix[row, match(markers, names(score_matrix))])
				for(formula_count in 1:nrow(coefficients)){
					formula_name <- rownames(coefficients)[formula_count]
					sum <- coefficients[formula_count,"Intercept"]
					for(marker_count in 1:length(markers)){
						marker_name <- markers[marker_count]
						sum <- sum + coefficients[formula_count, names(coefficients)==marker_name] * markers_score[marker_count]
					}
					score_matrix[row, paste0("pred_num_", formula_name)] <- sum
				}
			}
		}

		# combine data
		rownames(pheno) <- pheno$Sample_name
		merged_data <- merge(pheno, score_matrix, by=0)
		merged_data <- merged_data[,-1]
		number_cols <- c("Age_ori", "Progression_free_survival", "Disease_specific_survival", "OS")
		merged_data[number_cols] <- lapply(merged_data[number_cols], as.numeric)

		# optimal cut-off with Youden Index
		if(version %in% c("v4.2", "v4.3")){
			roc_data <- roc(merged_data[,eval(survival_type)], merged_data$pred_num)
			roc_cutoff <- coords(roc_data, "best")$threshold
			merged_data$roc_cutoff <- roc_cutoff
			merged_data$pred_cutoff <- ""
			merged_data[merged_data$pred_num < roc_cutoff, "pred_cutoff"] <- "FALSE"
			merged_data[merged_data$pred_num > roc_cutoff, "pred_cutoff"] <- "TRUE"
		}

		write.xlsx(merged_data, paste0(result_dir, "10_Models_plot/", version, "_Figure5/", survival_type, "/seed_", random_seed, "/merged_data.xlsx"), overwrite=TRUE)


		################################# 3) box plot of predicted scores
		if(version %in% c("v4.2", "v4.3")){
			box_plot(current_info, merged_data, path, survival_type, method_name, random_seed)
		}

		################################# 4) survival plot
		# format: path, current_info, merged_data, survival_plot_type, survival_plot_time
		survival_combos <- data.frame(survival_plot_type=c("Progression_beyond_T2_Progression", "Cancer_specific_satus_DOD_event", "Vital_status"), survival_plot_time=c("Progression_free_survival","Disease_specific_survival","OS"))
		for(survival_plot_type in survival_combos$survival_plot_type){
			
			survival_plot_time <- survival_combos[survival_combos$survival_plot_type==survival_plot_type,"survival_plot_time"]
			
			# plot data
			merged_data$survival_plot_type <- as.numeric(merged_data[,eval(survival_plot_type)])
			merged_data$survival_plot_time <- as.numeric(merged_data[,eval(survival_plot_time)])
	
			# KM curve plot
			if(version=="v4.4"){
				fit_survival <- survfit(Surv(survival_plot_time, survival_plot_type) ~ pred, data=merged_data)
			}else{
				fit_survival <- survfit(Surv(survival_plot_time, survival_plot_type) ~ pred_cutoff, data=merged_data)
			}
			

			if(survival_plot_time=="Progression_free_survival"){
				survival_plot <- ggsurvplot(fit_survival, 
					pval=TRUE, 
					conf.int=FALSE,
					ggtheme = theme_bw(),
					palette="aaas",
					ylim=c(0.5,1),
					pval.coord=c(0, 0.53),
					title=survival_plot_time)$plot +
					theme(plot.title=element_text(hjust=0.5))
			}else if(survival_plot_time=="Disease_specific_survival"){
				survival_plot <- ggsurvplot(fit_survival, 
					pval=TRUE, 
					conf.int=FALSE,
					ggtheme = theme_bw(),
					palette="aaas",
					ylim=c(0.7,1),
					pval.coord=c(0, 0.73),
					title=survival_plot_time)$plot +
					theme(plot.title=element_text(hjust=0.5))
			}else{
				survival_plot <- ggsurvplot(fit_survival, 
					pval=TRUE, 
					conf.int=FALSE,
					ggtheme = theme_bw(),
					palette="aaas",
					title=survival_plot_time)$plot +
					theme(plot.title=element_text(hjust=0.5))
			}

			pdf(paste0(path, survival_plot_type, "_", method_name, "_KMcurve.pdf"), onefile=FALSE)
			print(survival_plot)
			dev.off()

			jpeg(paste0(path, survival_plot_type, "_", method_name, "_KMcurve.jpg"))
			print(survival_plot)
			dev.off()

			topptx(survival_plot, paste0(path, survival_plot_type, "_", method_name, "_KMcurve.pptx"))
		}
	}
}
