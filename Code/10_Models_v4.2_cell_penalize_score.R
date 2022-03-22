rm(list=ls())
################################## BC progression ##################################
library(dplyr)
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(glmnet)
library(pROC)
library(doBy)

# set working directory
args <- commandArgs(T)
data_dir <- paste0(args[1], "Data/")
result_dir <- paste0(args[1], "Results/")
dir.create(paste0(result_dir, "10_Models_v4/"))

################################### Function module ################################## 
######### Function: Survival Analysis
survival_pack <- function(dataset_name, survival_type, survival_time){
	
	pheno_matrix <- pheno[pheno$Source_dataset==dataset_name,]
	candidate_cells <- c("Bcells", "DC", "Endothelial", "Fibroblasts","Macrophages", "Neutrophils", "NKcells", "Tcells_CD4+", "Tcells_CD8+")

	# for(method in c("ssGSEA", "Z-score")){
	for(method in c("ssGSEA")){
		if(method == "ssGSEA" && file.exists(paste0(result_dir, "9_Cell_Survival/", survival_time, "/", dataset_name, "_", method, "_merged_data_all.xlsx"))){
			method_name <- method
			score_matrix <- read.xlsx(paste0(result_dir, "9_Cell_Survival/", survival_time, "/", dataset_name, "_", method, "_merged_data_all.xlsx"), rowNames=FALSE)
			survival_plot_lr(dataset_name, method_name, score_matrix, pheno_matrix, survival_type, survival_time)
			survival_plot_penalize(score_matrix, pheno_matrix, survival_type, survival_time, candidate_cells, method_name)
		}else if(method == "Z-score" && file.exists(paste0(result_dir, "9_Cell_Survival/", survival_time, "/", dataset_name, "_", method, "_merged_data_all.xlsx"))){
			method_name <- method
			score_matrix <- read.xlsx(paste0(result_dir, "9_Cell_Survival/", survival_time, "/", dataset_name, "_", method, "_merged_data_all.xlsx"), rowNames=FALSE)
			survival_plot_lr(dataset_name, method_name, score_matrix, pheno_matrix, survival_type, survival_time)
			survival_plot_penalize(score_matrix, pheno_matrix, survival_type, survival_time, candidate_cells, method_name)
		}
	}
}

########################## Function: Cox & Logistic Regression ##########################
survival_plot_lr <- function(dataset_name, method_name, score_matrix, pheno_matrix, survival_type, survival_time){

	# 1) Cox analysis
	score_matrix <- na.omit(score_matrix)
	names(score_matrix) <- gsub("\\+", "", names(score_matrix))
	cell_cols <- names(score_matrix)[grep("_score$", names(score_matrix))]
	cox_fit <- coxph(as.formula(paste("Surv(survival_time, survival_type) ~", paste0(cell_cols, collapse=" + "))), data=score_matrix)
	forest_plot <- ggforest(model = cox_fit, data = score_matrix)

	pdf(paste0(result_dir, "10_Models_v4/10.1_", survival_time, "_", dataset_name, "_", method_name, "_forest_plot.pdf"))
	print(forest_plot)
	dev.off()

	# 2) Regression
	score_matrix$early_event <- ""
	score_matrix$survival_time <- as.numeric(score_matrix$survival_time)
	score_matrix[score_matrix$survival_time < 12,"early_event"] <- "1"
	score_matrix[score_matrix$survival_time >= 12,"early_event"] <- "0"
	score_matrix$early_event <- as.numeric(score_matrix$early_event)

	glm_fit <- glm(early_event ~ ., data = score_matrix[, c(cell_cols, "early_event")], family = binomial(link="logit"))
	glm_roc <- roc (score_matrix[,c("early_event")], predict(glm_fit, newdata = score_matrix, type="response"))
	ggROC <- ggroc(glm_roc) + annotate("text", x=0.5, y=1, label=paste0("AUC = ", glm_roc$auc))

	pdf(paste0(result_dir, "10_Models_v4/10.1_", survival_time, "_", dataset_name, "_", method_name, "_glm_ROC.pdf"))
	print(ggROC)
	dev.off()
}

########################## Function: Penalized model ##########################
survival_plot_penalize <- function(score_matrix, pheno_matrix, survival_type, survival_time, candidate_cells, method_name){

	# data
	model_data <- merge(pheno_matrix, score_matrix, by="Sample_name")
	number_cols <- c("Age_ori", "Progression_free_survival", "Disease_specific_survival", "OS", "survival_time", "Bcells_score", "DC_score", "Endothelial_score", "Fibroblasts_score","Macrophages_score", "Neutrophils_score", "NKcells_score", "Tcells_CD4+_score", "Tcells_CD8+_score")
	factor_cols <- c("Source_dataset", "Sex","Age","Staging","Grade73","Grade98","invasiveness","EORTC_risk_score_NMIBC","Tumor_size","Tumor_status","Growth_pattern","BCG_treatment","Chemotherapy","Radiotherapy","CIS_in_disease_course","Cystectomy","pN_at_Cystectomy","Recurrence","Progression","Progression_beyond_T2_Progression","Cancer_specific_satus_DOD_event","Vital_status","survival_type", "Bcells_level", "DC_level", "Endothelial_level", "Fibroblasts_level","Macrophages_level", "Neutrophils_level", "NKcells_level", "Tcells_CD4+_level", "Tcells_CD8+_level")
	model_data[number_cols] <- lapply(model_data[number_cols], as.numeric)
	model_data[factor_cols] <- lapply(model_data[factor_cols], factor)
	rownames(model_data) <- model_data$Sample_name
	# model_data <- model_data[,-grep("Sample_name", names(model_data))]
	model_data <- model_data[!is.na(model_data$survival_type),]

	# random seeds
	random_seeds <- sample(1:10000, 5000, replace=FALSE)
	for(random_count in 1:length(random_seeds)){
		random_seed <- random_seeds[random_count]

		# randomly split cohort 63 to training and test dataset
		set.seed(random_seed)
		train_set <- sampleBy(formula = ~ survival_type, frac=0.5, data=model_data)
		rownames(train_set) <- train_set$Sample_name
		
		# test set 1: all left samples
		test_set_left <- model_data[setdiff(rownames(model_data), rownames(train_set)),]
		rownames(test_set_left) <- test_set_left$Sample_name

		# test set 2: equal positive and negative samples
		pos_set <- test_set_left[test_set_left$survival_type=="TRUE",]
		neg_set <- test_set_left[test_set_left$survival_type=="FALSE",]
		test_set_balance <- rbind(pos_set, neg_set[sample(rownames(neg_set),nrow(pos_set)),])

		# model
		model_pack(paste0(result_dir, "10_Models_v4/"), survival_type, method_name, random_seed, candidate_cells, train_set, test_set_left, test_set_balance)
	}
}


#################################### Function: model ####################################
model_pack <- function(path, survival_type, method_name, random_seed, candidate_cells, train_set, test_set_left, test_set_balance){

	options(stringsAsFactors = F)
	#0.0load library----
	library(pbapply)
	library(mclust)
	library(miRBaseVersions.db)
	library(miRBaseConverter)
	library(pROC)
	library(SimDesign)
	library(glmnet)
	library(splitstackshape)
	library(randomForest)
	library(ggplot2)
	library(eoffice)
	library(ROCit)
	library(precrec)
	library(caret)
	
	display.progress = function (index, totalN, breakN=20) {
		if ( index %% ceiling(totalN/breakN)  ==0   ) {
			cat(paste(round(index*100/totalN), "% ", sep=""))
		}
	}
	
	####### output models record: random seed 
	current_line <- nrow(models) + 1
	models[current_line, "survival_type"] <<- survival_type
	models[current_line, "method_name"] <<- method_name
	models[current_line, "random_seed"] <<- random_seed

	# result dir by random_seed
	dir.create(paste0(path, "10.2_", survival_type, "_", method_name, "_random/seed_", random_seed), recursive=TRUE)
	setwd(paste0(path, "10.2_", survival_type, "_", method_name, "_random/seed_", random_seed))

	# record trainging and test subgroup
	write.xlsx(train_set, paste0("./0_training_set.xlsx"), rowNames=FALSE, overwrite=TRUE)
	write.xlsx(test_set_left, paste0("./0_test_set_left.xlsx"), rowNames=FALSE, overwrite=TRUE)
	write.xlsx(test_set_balance, paste0("./0_test_set_balance.xlsx"), rowNames=FALSE, overwrite=TRUE)

		########################## modeling
		#intersect_univ_lasso_RF_features
		candidate_cells 
		
		#pred_LR-----
		train_pred <- NULL
		test_left_pred <- NULL
		test_balance_pred <- NULL
		
		#pred_lasso_RF_AUC
		train_pred_AUC<-as.numeric()
		test_left_pred_AUC<-as.numeric()
		test_balance_pred_AUC<-as.numeric()
		
		# Penalized Elastic net regression
		x <- model.matrix(survival_type ~ Bcells_score + DC_score + Endothelial_score + Fibroblasts_score + Macrophages_score + Neutrophils_score + NKcells_score + `Tcells_CD4+_score` + `Tcells_CD8+_score`, data=train_set)[,-1]
		y <- train_set$survival_type

		set.seed(random_seed)
		model <- train(survival_type ~ ., data = train_set[,c(paste0(candidate_cells,"_score"),"survival_type")], method="glmnet", trControl = trainControl("cv", number = 10), tuneLength=10)
		get_best_result = function(caret_fit) {
			best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
			best_result = caret_fit$results[best, ]
			rownames(best_result) = NULL
			best_result
		}
		get_best_result(model)

		calc_acc = function(actual, predicted) {
			mean(actual == predicted)
		}
		calc_acc(actual = train_set$survival_type, predicted = predict(model, newdata = train_set))

		####### output models record: model 
		final_model_coef <- coef(model$finalModel, model$bestTune$lambda)
		models[current_line, "intercept"] <<- as.numeric(final_model_coef[1,ncol(final_model_coef)])
		models[current_line, "markers"] <<- paste(names(final_model_coef[-1,]), collapse="|")
		models[current_line, "coefficients"] <<- paste(as.numeric(final_model_coef[-1,ncol(final_model_coef)]), collapse="|")

		save(model, file="./1_model.Rdata")
		
		########################## 3.1.2 combine results
		train_pred <- rbind.data.frame(train_pred, data.frame(prob = predict(model, newdata = train_set), Sample_Group = train_set$survival_type, stringsAsFactors = F), stringsAsFactors = F)
		
		test_left_pred <- rbind.data.frame(test_left_pred, data.frame(prob = predict(model, newdata = test_set_left), Sample_Group = test_set_left$survival_type, stringsAsFactors = F), stringsAsFactors = F)

		test_balance_pred <- rbind.data.frame(test_balance_pred, data.frame(prob = predict(model, newdata = test_set_balance), Sample_Group = test_set_balance$survival_type, stringsAsFactors = F), stringsAsFactors = F)
		
		# 3.2 AUROC
		# 3.2.1 AUROC
		# train-----
		train_pred_ROC<-roc(train_set[,c("survival_type")],predict(model, newdata = train_set, type="prob")[,1])
		train_pred_ROC
		train_pred_AUC<- append(train_pred_AUC,formattable::formattable(as.numeric(pROC::auc(train_pred_ROC)),3,format = "f"))
		train_pred_AUC
		# test lfet----
		test_left_pred_ROC<-roc(test_set_left[,c("survival_type")],predict(model, newdata = test_set_left,type="prob")[,1])
		test_left_pred_ROC
		test_left_pred_AUC<- append(test_left_pred_AUC,formattable::formattable(as.numeric(pROC::auc(test_left_pred_ROC)),3,format = "f"))
		test_left_pred_AUC
		# test balance ----
		test_balance_pred_ROC<-roc(test_set_balance[,c("survival_type")],predict(model, newdata = test_set_balance,type="prob")[,1])
		test_balance_pred_ROC
		test_balance_pred_AUC<- append(test_balance_pred_AUC,formattable::formattable(as.numeric(pROC::auc(test_balance_pred_ROC)),3,format = "f"))
		test_balance_pred_AUC

		#3.2.2 AUC----
		train_pred_ROC
		test_left_pred_ROC
		test_balance_pred_ROC

		roclist<-list(train_pred_ROC=train_pred_ROC,test_left_pred_ROC=test_left_pred_ROC,test_balance_pred_ROC=test_balance_pred_ROC)
			
		# train_ROC
		train_ROC<-paste("Training	Set = ",formattable::formattable(as.numeric(pROC::auc(train_pred_ROC)),3,format = "f"))
		train_ROC
		# test_left ROC
		test_left_ROC<-paste("Test Set = ",formattable::formattable(as.numeric(pROC::auc(test_left_pred_ROC)),3,format = "f"))
		test_left_ROC
		
		# test_balance ROC
		test_balance_ROC<-paste("Test Set = ",formattable::formattable(as.numeric(pROC::auc(test_balance_pred_ROC)),3,format = "f"))
		test_balance_ROC

		########################## 3.2.3.2 ROC curve
		ggROC<-ggroc(roclist,legacy.axes = TRUE,size=1.5)+
			theme_bw()+theme(legend.background =element_blank())+
			scale_color_manual(labels = c(train_ROC, test_left_ROC, test_balance_ROC), values = c("#7fc97f", "#beaed4", "#fdc086"))+
			geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 1.5,color="grey", linetype="dashed")+
			theme(panel.background=element_rect(fill="white",colour="black",size=0.5))+
			theme(legend.position = c(0.75,0.15),legend.key.width = unit(1, "cm"),
				  legend.key.height = unit(1.5, "cm"))+	#修改图例的item的间距
			theme(legend.title=element_blank())+
			theme(legend.text = element_text(size = rel(1.5),color="black"))+
			theme(axis.text=element_text(size=rel(2),colour ="black"),#加粗刻度标签
				axis.title=element_text(face="bold", size=rel(2), color="black"))+
			ylab("Sensitivity")+xlab("1-Specificity")

		ggsave(ggROC,width=9,height=9,file = "./1_model_ROC.svg")
		ggsave(ggROC,width=9,height=9,file ="./1_model_ROC.pdf")
		topptx(ggROC, "./1_model_ROC.pptx")

		####### output models record: model performance (train & test)
		models[current_line, "train_AUC"] <<- as.numeric(train_pred_ROC$auc)
		models[current_line, "test_left_AUC"] <<- as.numeric(test_left_pred_ROC$auc)
		models[current_line, "test_balance_AUC"] <<- as.numeric(test_balance_pred_ROC$auc)
}

######################################## main ########################################
# output file 1: model record
models <- data.frame(matrix(nrow=0, ncol=10))
names(models) <- c("survival_type", "method_name", "random_seed", "intercept", "markers", "coefficients", "train_AUC", "test_left_AUC", "test_balance_AUC", "vali_AUC")

pheno <- read.xlsx(paste0(data_dir, "combined/Clinicopathologic_v2.xlsx"), na.strings=c("","NA","Nan","#N/A"))
datasets <- c("E-MTAB-4321", "GSE32894", "GSE13507")

# function format: dataset_name, survival_type, survival_time)
survival_pack("E-MTAB-4321", "Progression_beyond_T2_Progression", "Progression_free_survival")
survival_pack("GSE32894", "Cancer_specific_satus_DOD_event", "Disease_specific_survival")
survival_pack("GSE13507", "Vital_status", "OS")

# write to file
write.xlsx(models, paste0(result_dir, "10_Models_v4/10.2_Modes_AUC.xlsx"), overwrite=TRUE)
