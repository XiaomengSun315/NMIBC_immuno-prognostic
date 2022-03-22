# 2022.03.20
# modified according to reviewer1's latest comment: use original expression matrix

# 2022.02.27
# modified according to reviewer1's suggestion: update with "9_Cell-score_Survival_v3.R"

# 2021.10.20
# Penalized Elastic net regression: use cell enrichment score instead of low/high level

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
library(InformationValue)

# set working directory
data_dir <- "C:/0_xmsun/xmsun/Graduate/20210224_NMIBC/Data/"
result_dir <- "C:/0_xmsun/xmsun/Graduate/20210224_NMIBC/Results/"
dir.create(paste0(result_dir, "10_Models_v4/"))

################################### Function module ################################## 
######### Function: Survival Analysis
survival_pack <- function(dataset, survival_type, survival_time){
	
	# retrive score matrix with available survival data
	pheno_matrix_type <- pheno[!is.na(pheno[,c(survival_type)]),]
	pheno_matrix_both <- pheno[(!is.na(pheno[,c(survival_type)]))&(!is.na(pheno[,c(survival_time)])),]

	score_matrix_all <- read.xlsx(paste0(result_dir, "8_Cell_scores/8.3_ssGSEA_matrix.xlsx"))
	pheno_matrix_type <- pheno_matrix_type[pheno_matrix_type$Sample_name %in% names(score_matrix_all),]
	pheno_matrix_both <- pheno_matrix_both[pheno_matrix_both$Sample_name %in% names(score_matrix_all),]

	# score matirx with only survival type
	score_matrix_type <- data.frame(t(score_matrix_all[,match(pheno_matrix_type$Sample_name, names(score_matrix_all))]))
	score_matrix_type <- cbind(score_matrix_type, pheno_matrix_type$Sample_name, pheno_matrix_type[,survival_type], pheno_matrix_type[,survival_time])
	names(score_matrix_type) <- c(paste0(score_matrix_all[,1], "_score"), "Sample_name", "survival_type", "survival_time")
	
	# score matirx with survival type & survival time
	score_matrix_both <- data.frame(t(score_matrix_all[,match(pheno_matrix_both$Sample_name, names(score_matrix_all))]))
	score_matrix_both <- cbind(score_matrix_both, pheno_matrix_both[,survival_type], pheno_matrix_both[,survival_time])
	names(score_matrix_both) <- c(paste0(score_matrix_all[,1], "_score"), "survival_type", "survival_time")

	candidate_cells <- c("Bcells", "DC", "Endothelial", "Fibroblasts","Macrophages", "Neutrophils", "NKcells", "Tcells_CD4+", "Tcells_CD8+")

	# for(method in c("ssGSEA", "Z-score")){
	for(method in c("ssGSEA")){
		if(method == "ssGSEA"){
			method_name <- method
			survival_plot_lr(method_name, score_matrix_both, pheno_matrix_both, survival_type, survival_time)
			survival_plot_penalize(dataset, score_matrix_type, pheno_matrix_type, survival_type, survival_time, candidate_cells, method_name)
		}else if(method == "Z-score"){
			method_name <- method
			survival_plot_lr(method_name, score_matrix_both, pheno_matrix_both, survival_type, survival_time)
			survival_plot_penalize(dataset, score_matrix_type, pheno_matrix_type, survival_type, survival_time, candidate_cells, method_name)
		}
	}
}

########################## Function: Cox & Logistic Regression ##########################
survival_plot_lr <- function(method_name, score_matrix, pheno_matrix, survival_type, survival_time){

	# 1) Cox analysis
	score_matrix <- na.omit(score_matrix)
	names(score_matrix) <- gsub("\\+", "", names(score_matrix))
	cell_cols <- names(score_matrix)[grep("_score$", names(score_matrix))]
	cox_fit <- coxph(as.formula(paste("Surv(survival_time, survival_type) ~", paste0(cell_cols, collapse=" + "))), data=score_matrix)
	forest_plot <- ggforest(model = cox_fit, data = score_matrix)

	pdf(paste0(result_dir, "10_Models_v4/10.1_", survival_time, "_", method_name, "_forest_plot.pdf"))
	print(forest_plot)
	dev.off()

	# 2) Regression
	score_matrix$early_event <- ""
	score_matrix$survival_time <- as.numeric(score_matrix$survival_time)
	score_matrix[score_matrix$survival_time < 12,"early_event"] <- "1"
	score_matrix[score_matrix$survival_time >= 12,"early_event"] <- "0"
	score_matrix$early_event <- as.numeric(score_matrix$early_event)

	glm_fit <- glm(survival_type ~ ., data = score_matrix[, c(cell_cols, "survival_type")], family = binomial(link="logit"))

	# models[, "intercept"] <- as.numeric(glm_fit$coefficients[1])
	# models[, "markers"] <- paste(names(glm_fit$coefficients[-1]), collapse="|")
	# models[, "coefficients"] <- paste(as.numeric(glm_fit$coefficients[-1]), collapse="|")

	glm_roc <- roc(score_matrix[,c("survival_type")], predict(glm_fit, newdata = score_matrix, type="response"))
	ggROC <- ggroc(glm_roc) + annotate("text", x=0.5, y=1, label=paste0("AUC = ", glm_roc$auc))

	pdf(paste0(result_dir, "10_Models_v4/10.1_", survival_time, "_", method_name, "_glm_ROC.pdf"))
	print(ggROC)
	dev.off()
}

########################## Function: Penalized model ##########################
survival_plot_penalize <- function(dataset, score_matrix, pheno_matrix, survival_type, survival_time, candidate_cells, method_name){

	# data
	model_data <- merge(pheno_matrix, score_matrix, by="Sample_name")
	number_cols <- c("Age_ori", "Progression_free_survival", "Disease_specific_survival", "OS", "survival_time", "Bcells_score", "DC_score", "Endothelial_score", "Fibroblasts_score","Macrophages_score", "Neutrophils_score", "NKcells_score", "Tcells_CD4+_score", "Tcells_CD8+_score", "Progression_beyond_T2_Progression_num", "Cancer_specific_satus_DOD_event_num", "Vital_status_num")
	factor_cols <- c("Source_dataset", "Sex","Age","Staging","Grade73","Grade98","invasiveness","EORTC_risk_score_NMIBC","Tumor_size","Tumor_status","Growth_pattern","BCG_treatment","Chemotherapy","Radiotherapy","CIS_in_disease_course","Cystectomy","pN_at_Cystectomy","Recurrence","Progression","Progression_beyond_T2_Progression","Cancer_specific_satus_DOD_event","Vital_status","survival_type")
	model_data[number_cols] <- lapply(model_data[number_cols], as.numeric)
	model_data[factor_cols] <- lapply(model_data[factor_cols], factor)
	model_data$Staging <- factor(model_data$Staging, levels=c("Ta", "Ta-1", "T1", "Tis"))
	rownames(model_data) <- model_data$Sample_name
	# model_data <- model_data[,-grep("Sample_name", names(model_data))]
	model_data <- model_data[!is.na(model_data$survival_type),]
	model_data <- model_data[model_data$survival_time>0 | is.na(model_data$survival_tim),]

	# use pre-defined dataset to build and validate PFS/DFS/OS models
	model_data <- model_data[model_data$Source_dataset==dataset,]

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
		model_pack(paste0(result_dir, "10_Models_v4/"), survival_type, method_name, random_seed, candidate_cells, model_data, train_set, test_set_left, test_set_balance)
	}
}

#################################### Function: ROC plot ####################################
ROC_plot <- function(work_dir, model, cvfit, model_type, train_set, test_set_left, test_set_balance, x_bino_train, x_bino_test_left, x_bino_test_balance){

	setwd(work_dir)

	train_pred_AUC <- as.numeric()
	test_left_pred_AUC <- as.numeric()
	test_balance_pred_AUC <- as.numeric()
	
	train_pred_ROC <- roc(train_set[,c("survival_type")], as.numeric(predict(model, newx = x_bino_train, s=cvfit$lambda.min, type="response")))
	train_pred_ROC
	train_pred_AUC <- append(train_pred_AUC,formattable::formattable(as.numeric(pROC::auc(train_pred_ROC)),3,format = "f"))
	train_pred_AUC
	
	# test lfet----
	test_left_pred_ROC <- roc(test_set_left[,c("survival_type")], as.numeric(predict(model, newx = x_bino_test_left, s=cvfit$lambda.min, type="response")))
	test_left_pred_ROC
	test_left_pred_AUC <- append(test_left_pred_AUC,formattable::formattable(as.numeric(pROC::auc(test_left_pred_ROC)),3,format = "f"))
	test_left_pred_AUC
	
	# test balance ----
	test_balance_pred_ROC <- roc(test_set_balance[,c("survival_type")], as.numeric(predict(model, newx = x_bino_test_balance, s=cvfit$lambda.min, type="response")))
	test_balance_pred_ROC
	test_balance_pred_AUC <- append(test_balance_pred_AUC,formattable::formattable(as.numeric(pROC::auc(test_balance_pred_ROC)),3,format = "f"))
	test_balance_pred_AUC

	#3.2.2查看AUC----
	train_pred_ROC
	test_left_pred_ROC
	test_balance_pred_ROC

	roclist<-list(train_pred_ROC=train_pred_ROC,test_left_pred_ROC=test_left_pred_ROC,test_balance_pred_ROC=test_balance_pred_ROC)
			
	# train_ROC
	train_ROC <- paste("Training	Set = ",formattable::formattable(as.numeric(pROC::auc(train_pred_ROC)),3,format = "f"))
	train_ROC
	# test_left ROC
	test_left_ROC <- paste("Test Set = ",formattable::formattable(as.numeric(pROC::auc(test_left_pred_ROC)),3,format = "f"))
	test_left_ROC
		
	# test_balance ROC
	test_balance_ROC <- paste("Test Set = ",formattable::formattable(as.numeric(pROC::auc(test_balance_pred_ROC)),3,format = "f"))
	test_balance_ROC

	########################## 3.2.3.2 plot ROC curve
	ggROC <- ggroc(roclist,legacy.axes = TRUE,size=1.5) +
		theme_bw() +
		theme(legend.background =element_blank()) +
		scale_color_manual(labels = c(train_ROC, test_left_ROC, test_balance_ROC), values = c("#7fc97f", "#beaed4", "#fdc086")) +
		#geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), size = 1.5,color="grey", linetype="dashed") +
		geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 1.5,color="grey", linetype="dashed") +
		theme(panel.background=element_rect(fill="white",colour="black",size=0.5)) +
		theme(legend.position = c(0.75,0.15),legend.key.width = unit(1, "cm"), legend.key.height = unit(1.5, "cm")) +
		theme(legend.title=element_blank()) +
		theme(legend.text = element_text(size = rel(1.5),color="black")) +
		theme(axis.text=element_text(size=rel(2),colour ="black"), axis.title=element_text(face="bold", size=rel(2), color="black")) +
		ylab("Sensitivity") +
		xlab("1-Specificity")

	# ggsave(ggROC,width=9,height=9,file = paste0("./1_model_", model_type, "_ROC.svg"))
	ggsave(ggROC,width=9,height=9,file = paste0("./1_model_", model_type, "_ROC.pdf"))
	topptx(ggROC, paste0("./1_model_", model_type, "_ROC.pptx"))

	return(c(train_pred_AUC, test_left_pred_AUC, test_balance_pred_AUC))
}

#################################### Function: model ####################################
model_pack <- function(path, survival_type, method_name, random_seed, candidate_cells, model_data, train_set, test_set_left, test_set_balance){

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

	####### output models record: random seed 
	current_line <- nrow(models) + 1
	models[current_line, "survival_type"] <<- survival_type
	models[current_line, "method_name"] <<- method_name
	models[current_line, "random_seed"] <<- random_seed

	####### result dir by random_seed
	work_dir <- paste0(path, "10.2_", survival_type, "_", method_name, "_random/seed_", random_seed)
	dir.create(work_dir, recursive=TRUE)
	setwd(work_dir)

	########################## Lasso/Ridge/Penalized Elastic net regression models
	####### x & y matrix for binomail y regression
	x <- model.matrix(survival_type ~ Bcells_score + DC_score + Endothelial_score + Fibroblasts_score + Macrophages_score + Neutrophils_score + NKcells_score + `Tcells_CD4+_score` + `Tcells_CD8+_score`, data=model_data)[,-1]

	x_bino <- x
	x_bino_train <- x_bino[match(rownames(train_set), rownames(x_bino)),]
	x_bino_test_left <- x_bino[match(rownames(test_set_left), rownames(x_bino)),]
	x_bino_test_balance <- x_bino[match(rownames(test_set_balance), rownames(x_bino)),]
	
	y_bino <- model_data[,paste0(survival_type, "_num")]
	names(y_bino) <- rownames(model_data)
	y_bino_train <- y_bino[match(rownames(train_set), names(y_bino))]
	y_bino_test_left <- y_bino[match(rownames(test_set_left), names(y_bino))]
	y_bino_test_balance <- y_bino[match(rownames(test_set_balance), names(y_bino))]
	
	# x, y for cox regression
	y_surv <- data.frame(time=as.numeric(model_data$survival_time), status=as.character(model_data$survival_type))
	y_surv[y_surv$status=="TRUE", "status"] <- 1
	y_surv[y_surv$status=="FALSE", "status"] <- 0
	y_surv <- apply(y_surv, 2, as.numeric)
	rownames(y_surv) <- rownames(model_data)
	y_surv <- na.omit(y_surv)

	x_surv <- x[rownames(x) %in% rownames(y_surv),]
	x_surv_train <- x_surv[rownames(x_surv) %in% rownames(train_set),]
	x_surv_test_left <- x_surv[rownames(x_surv) %in% rownames(test_set_left),]
	x_surv_test_balance <- x_surv[rownames(x_surv) %in% rownames(test_set_balance),]
	
	y_surv_train <- y_surv[match(rownames(x_surv_train), rownames(y_surv)),]
	y_surv_test_left <- y_surv[match(rownames(x_surv_test_left), rownames(y_surv)),]
	y_surv_test_balance <- y_surv[match(rownames(x_surv_test_balance), rownames(y_surv)),]

	####### models
	set.seed(random_seed)
	model_bino_elnet <- glmnet(x_bino_train, y_bino_train, family="binomial", alpha=0.5)
	model_bino_lasso <- glmnet(x_bino_train, y_bino_train, family="binomial", alpha=1)
	model_bino_ridge <- glmnet(x_bino_train, y_bino_train, family="binomial", alpha=0)
	# model_cox_elnet <- glmnet(x_surv_train, y_surv_train, family="cox", alpha=0.5)
	# model_cox_lasso <- glmnet(x_surv_train, y_surv_train, family="cox", alpha=1)
	# model_cox_ridge <- glmnet(x_surv_train, y_surv_train, family="cox", alpha=0)

	cvfit_bino_elnet <- cv.glmnet(x_bino_train, y_bino_train, family="binomial", alpha=0.5, nfold=5)
	cvfit_bino_lasso <- cv.glmnet(x_bino_train, y_bino_train, family="binomial", alpha=1, nfold=5)
	cvfit_bino_ridge <- cv.glmnet(x_bino_train, y_bino_train, family="binomial", alpha=0, nfold=5)
	# cvfit_cox_elnet <- cv.glmnet(x_surv_train, y_surv_train, family="cox", alpha=0.5, nfold=5)
	# cvfit_cox_lasso <- cv.glmnet(x_surv_train, y_surv_train, family="cox", alpha=1, nfold=5)
	# cvfit_cox_ridge <- cv.glmnet(x_surv_train, y_surv_train, family="cox", alpha=0, nfold=5)

	coef(cvfit_bino_elnet, s="lambda.min")
	coef(cvfit_bino_lasso, s="lambda.min")
	coef(cvfit_bino_ridge, s="lambda.min")
	# coef(cvfit_cox_elnet, s="lambda.min")
	# coef(cvfit_cox_lasso, s="lambda.min")
	# coef(cvfit_cox_ridge, s="lambda.min")

	# actual_y <- y_bino
	# predict_y <- predict(model_bino_test, newx = x_bino, s=cvfit_bino_test$lambda.min, type="response")
	# misClassError(actual_y, predict_y)
	# plotROC(actual_y, as.numeric(predict_y))

	######################### ROC plot for the final model
	# input value: work_dir, model, model_type, train_set, test_set_left, test_set_balance
	ridge_AUCs <- ROC_plot(work_dir, model_bino_ridge, cvfit_bino_ridge, "ridge", train_set, test_set_left, test_set_balance, x_bino_train, x_bino_test_left, x_bino_test_balance)
	lasso_AUCs <- ROC_plot(work_dir, model_bino_lasso, cvfit_bino_lasso, "lasso", train_set, test_set_left, test_set_balance, x_bino_train, x_bino_test_left, x_bino_test_balance)
	elnet_AUCs <- ROC_plot(work_dir, model_bino_elnet, cvfit_bino_elnet, "elnet", train_set, test_set_left, test_set_balance, x_bino_train, x_bino_test_left, x_bino_test_balance)

	########################## output models record: model performance (train & test)
	models[current_line, "bino_elnet_train_AUC"] <<- elnet_AUCs[1]
	models[current_line, "bino_elnet_test_left_AUC"] <<- elnet_AUCs[2]
	models[current_line, "bino_elnet_test_balance_AUC"] <<- elnet_AUCs[3]
	models[current_line, "bino_elnet_misClass"] <<- misClassError(y_bino, predict(model_bino_elnet, newx = x_bino, s=cvfit_bino_elnet$lambda.min, type="response"))
	pdf("./ROC_bino_elnet.pdf")
	plotROC(y_bino, predict(model_bino_elnet, newx = x_bino, s=cvfit_bino_elnet$lambda.min, type="response"))
	dev.off()

	models[current_line, "bino_lasso_train_AUC"] <<- lasso_AUCs[1]
	models[current_line, "bino_lasso_test_left_AUC"] <<- lasso_AUCs[2]
	models[current_line, "bino_lasso_test_balance_AUC"] <<- lasso_AUCs[3]
	models[current_line, "bino_lasso_misClass"] <<- misClassError(y_bino, predict(model_bino_lasso, newx = x_bino, s=cvfit_bino_lasso$lambda.min, type="response"))
	pdf("./ROC_bino_lasso.pdf")
	plotROC(y_bino, predict(model_bino_lasso, newx = x_bino, s=cvfit_bino_lasso$lambda.min, type="response"))
	dev.off()

	models[current_line, "bino_ridge_train_AUC"] <<- ridge_AUCs[1]
	models[current_line, "bino_ridge_test_left_AUC"] <<- ridge_AUCs[2]
	models[current_line, "bino_ridge_test_balance_AUC"] <<- ridge_AUCs[3]
	models[current_line, "bino_ridge_misClass"] <<- misClassError(y_bino, predict(model_bino_ridge, newx = x_bino, s=cvfit_bino_ridge$lambda.min, type="response"))
	pdf("./ROC_bino_ridge.pdf")
	plotROC(y_bino, predict(model_bino_ridge, newx = x_bino, s=cvfit_bino_ridge$lambda.min, type="response"))
	dev.off()

	# models[current_line, "cox_elnet_train_C"] <<- max(assess.glmnet(model_cox_elnet, newx=x_surv_train, newy=y_surv_train)$C)
	# models[current_line, "cox_elnet_test_left_C"] <<- max(assess.glmnet(model_cox_elnet, newx=x_surv_test_left, newy=y_surv_test_left)$C)
	# models[current_line, "cox_elnet_test_balance_C"] <<- max(assess.glmnet(model_cox_elnet, newx=x_surv_test_balance, newy=y_surv_test_balance)$C)
	# models[current_line, "cox_elnet_misClass"] <<- misClassError(y_surv, as.numeric(predict(model_cox_elnet, newx = x_surv, s=cvfit_cox_elnet$lambda.min, type="response")))
	# pdf("./ROC_cox_elnet.pdf")
	# plotROC(y_surv, as.numeric(predict(model_cox_elnet, newx = x_surv, s=cvfit_cox_elnet$lambda.min, type="response")))
	# dev.off()

	# models[current_line, "cox_lasso_train_C"] <<- max(assess.glmnet(model_cox_lasso, newx=x_surv_train, newy=y_surv_train)$C)
	# models[current_line, "cox_lasso_test_left_C"] <<- max(assess.glmnet(model_cox_lasso, newx=x_surv_test_left, newy=y_surv_test_left)$C)
	# models[current_line, "cox_lasso_test_balance_C"] <<- max(assess.glmnet(model_cox_lasso, newx=x_surv_test_balance, newy=y_surv_test_balance)$C)
	# models[current_line, "cox_lasso_misClass"] <<- misClassError(y_surv, as.numeric(predict(model_cox_lasso, newx = x_surv, s=cvfit_cox_lasso$lambda.min, type="response")))
	# pdf("./ROC_cox_lasso.pdf")
	# plotROC(y_surv, as.numeric(predict(model_cox_lasso, newx = x_surv, s=cvfit_cox_lasso$lambda.min, type="response")))
	# dev.off()

	# models[current_line, "cox_ridge_train_C"] <<- max(assess.glmnet(model_cox_ridge, newx=x_surv_train, newy=y_surv_train)$C)
	# models[current_line, "cox_ridge_test_left_C"] <<- max(assess.glmnet(model_cox_ridge, newx=x_surv_test_left, newy=y_surv_test_left)$C)
	# models[current_line, "cox_ridge_test_balance_C"] <<- max(assess.glmnet(model_cox_ridge, newx=x_surv_test_balance, newy=y_surv_test_balance)$C)
	# models[current_line, "cox_ridge_misClass"] <<- misClassError(y_surv, as.numeric(predict(model_cox_ridge, newx = x_surv, s=cvfit_cox_ridge$lambda.min, type="response")))
	# pdf("./ROC_cox_ridge.pdf")
	# plotROC(y_surv, as.numeric(predict(model_cox_ridge, newx = x_surv, s=cvfit_cox_ridge$lambda.min, type="response")))
	# dev.off()

	# datasets and predicted results: training, test_left, test_balance
	current_col_train <- ncol(train_set)
	train_set <- cbind(train_set, 
		data.frame(prob_numb = predict(model_bino_elnet, newx = x_bino_train, s=cvfit_bino_elnet$lambda.min, type="response"), prob = predict(model_bino_elnet, newx = x_bino_train, s=cvfit_bino_elnet$lambda.min, type="class")),
		data.frame(prob_numb = predict(model_bino_lasso, newx = x_bino_train, s=cvfit_bino_lasso$lambda.min, type="response"), prob = predict(model_bino_lasso, newx = x_bino_train, s=cvfit_bino_lasso$lambda.min, type="class")),
		data.frame(prob_numb = predict(model_bino_ridge, newx = x_bino_train, s=cvfit_bino_ridge$lambda.min, type="response"), prob = predict(model_bino_ridge, newx = x_bino_train, s=cvfit_bino_ridge$lambda.min, type="class")))
	# cox_train <- cbind(
	# 	data.frame(prob_numb = predict(model_cox_elnet, newx = x_surv_train, s=cvfit_cox_elnet$lambda.min, type="response")),
	# 	data.frame(prob_numb = predict(model_cox_lasso, newx = x_surv_train, s=cvfit_cox_lasso$lambda.min, type="response")),
	# 	data.frame(prob_numb = predict(model_cox_ridge, newx = x_surv_train, s=cvfit_cox_ridge$lambda.min, type="response")))
	# train_set <- merge(train_set, cox_train, by=0)
	# names(train_set)[(ncol(train_set)-8):ncol(train_set)] <- c("bino_elnet_prob", "bino_elnet_pred", "bino_lasso_prob", "bino_lasso_pred", "bino_ridge_prob", "bino_ridge_pred", "cox_elnet_prob", "cox_lasso_prob", "cox_ridge_prob")
	names(train_set)[(ncol(train_set)-5):ncol(train_set)] <- c("bino_elnet_prob", "bino_elnet_pred", "bino_lasso_prob", "bino_lasso_pred", "bino_ridge_prob", "bino_ridge_pred")

	current_col_test_left <- ncol(test_set_left)
	test_set_left <- cbind(test_set_left, 
		data.frame(prob_numb = predict(model_bino_elnet, newx = x_bino_test_left, s=cvfit_bino_elnet$lambda.min, type="response"), prob = predict(model_bino_elnet, newx = x_bino_test_left, s=cvfit_bino_elnet$lambda.min, type="class")),
		data.frame(prob_numb = predict(model_bino_lasso, newx = x_bino_test_left, s=cvfit_bino_lasso$lambda.min, type="response"), prob = predict(model_bino_lasso, newx = x_bino_test_left, s=cvfit_bino_lasso$lambda.min, type="class")),
		data.frame(prob_numb = predict(model_bino_ridge, newx = x_bino_test_left, s=cvfit_bino_ridge$lambda.min, type="response"), prob = predict(model_bino_ridge, newx = x_bino_test_left, s=cvfit_bino_ridge$lambda.min, type="class")))
	# cox_left <- cbind(
	# 	data.frame(prob_numb = predict(model_cox_elnet, newx = x_surv_test_left, s=cvfit_cox_elnet$lambda.min, type="response")),
	# 	data.frame(prob_numb = predict(model_cox_lasso, newx = x_surv_test_left, s=cvfit_cox_lasso$lambda.min, type="response")),
	# 	data.frame(prob_numb = predict(model_cox_ridge, newx = x_surv_test_left, s=cvfit_cox_ridge$lambda.min, type="response")))
	# test_set_left <- merge(test_set_left, cox_left, by=0)
	# names(test_set_left)[(ncol(test_set_left)-8):ncol(test_set_left)] <- c("bino_elnet_prob", "bino_elnet_pred", "bino_lasso_prob", "bino_lasso_pred", "bino_ridge_prob", "bino_ridge_pred", "cox_elnet_prob", "cox_lasso_prob", "cox_ridge_prob")
	names(test_set_left)[(ncol(test_set_left)-5):ncol(test_set_left)] <- c("bino_elnet_prob", "bino_elnet_pred", "bino_lasso_prob", "bino_lasso_pred", "bino_ridge_prob", "bino_ridge_pred")

	current_col_test_balance <- ncol(test_set_balance)
	test_set_balance <- cbind(test_set_balance, 
		data.frame(prob_numb = predict(model_bino_elnet, newx = x_bino_test_balance, s=cvfit_bino_elnet$lambda.min, type="response"), prob = predict(model_bino_elnet, newx = x_bino_test_balance, s=cvfit_bino_elnet$lambda.min, type="class")),
		data.frame(prob_numb = predict(model_bino_lasso, newx = x_bino_test_balance, s=cvfit_bino_lasso$lambda.min, type="response"), prob = predict(model_bino_lasso, newx = x_bino_test_balance, s=cvfit_bino_lasso$lambda.min, type="class")),
		data.frame(prob_numb = predict(model_bino_ridge, newx = x_bino_test_balance, s=cvfit_bino_ridge$lambda.min, type="response"), prob = predict(model_bino_ridge, newx = x_bino_test_balance, s=cvfit_bino_ridge$lambda.min, type="class")))
	# cox_balance <-cbind(
	# 	data.frame(prob_numb = predict(model_cox_elnet, newx = x_surv_test_balance, s=cvfit_cox_elnet$lambda.min, type="response")),
	# 	data.frame(prob_numb = predict(model_cox_lasso, newx = x_surv_test_balance, s=cvfit_cox_lasso$lambda.min, type="response")),
	# 	data.frame(prob_numb = predict(model_cox_ridge, newx = x_surv_test_balance, s=cvfit_cox_ridge$lambda.min, type="response")))
	# test_set_balance <- merge(test_set_balance, cox_balance, by=0)
	# names(test_set_balance)[(ncol(test_set_balance)-8):ncol(test_set_balance)] <- c("bino_elnet_prob", "bino_elnet_pred", "bino_lasso_prob", "bino_lasso_pred", "bino_ridge_prob", "bino_ridge_pred", "cox_elnet_prob", "cox_lasso_prob", "cox_ridge_prob")
	names(test_set_balance)[(ncol(test_set_balance)-5):ncol(test_set_balance)] <- c("bino_elnet_prob", "bino_elnet_pred", "bino_lasso_prob", "bino_lasso_pred", "bino_ridge_prob", "bino_ridge_pred")
	
	########################## output record
	# record trainging and test datasets
	write.xlsx(train_set, paste0("./0_training_set.xlsx"), rowNames=FALSE, overwrite=TRUE)
	write.xlsx(test_set_left, paste0("./0_test_set_left.xlsx"), rowNames=FALSE, overwrite=TRUE)
	write.xlsx(test_set_balance, paste0("./0_test_set_balance.xlsx"), rowNames=FALSE, overwrite=TRUE)

	# save(x_bino, x_bino_train, x_bino_test_left, x_bino_test_balance, x_surv, x_surv_train, x_surv_test_left, x_surv_test_balance, y_bino, y_bino_train, y_bino_test_left, y_bino_test_balance, y_surv, y_surv_train, y_surv_test_left, y_surv_test_balance, model_bino_elnet, cvfit_bino_elnet, model_bino_lasso, cvfit_bino_lasso, model_bino_ridge, cvfit_bino_ridge, model_cox_elnet, cvfit_cox_elnet, model_cox_lasso, cvfit_cox_lasso, model_cox_ridge, cvfit_cox_ridge, file="./1_models.Rdata")
	save(x_bino, x_bino_train, x_bino_test_left, x_bino_test_balance, x_surv, x_surv_train, x_surv_test_left, x_surv_test_balance, y_bino, y_bino_train, y_bino_test_left, y_bino_test_balance, y_surv, y_surv_train, y_surv_test_left, y_surv_test_balance, model_bino_elnet, cvfit_bino_elnet, model_bino_lasso, cvfit_bino_lasso, model_bino_ridge, cvfit_bino_ridge, file="./1_models.Rdata")
}

######################################## main ########################################
# output file 1: model record
# models <- data.frame(matrix(nrow=0, ncol=27))
# names(models) <- c("survival_type", "method_name", "random_seed", "bino_elnet_train_AUC", "bino_elnet_test_left_AUC", "bino_elnet_test_balance_AUC", "bino_elnet_misClass", "bino_lasso_train_AUC", "bino_lasso_test_left_AUC", "bino_lasso_test_balance_AUC", "bino_lasso_misClass", "bino_ridge_train_AUC", "bino_ridge_test_left_AUC", "bino_ridge_test_balance_AUC", "bino_ridge_misClass", "cox_elnet_train_C", "cox_elnet_test_left_C", "cox_elnet_test_balance_C", "cox_elnet_misClass", "cox_lasso_train_C", "cox_lasso_test_left_C", "cox_lasso_test_balance_C", "cox_lasso_misClass", "cox_ridge_train_C", "cox_ridge_test_left_C", "cox_ridge_test_balance_C", "cox_ridge_misClass")
models <- data.frame(matrix(nrow=0, ncol=15))
names(models) <- c("survival_type", "method_name", "random_seed", "bino_elnet_train_AUC", "bino_elnet_test_left_AUC", "bino_elnet_test_balance_AUC", "bino_elnet_misClass", "bino_lasso_train_AUC", "bino_lasso_test_left_AUC", "bino_lasso_test_balance_AUC", "bino_lasso_misClass", "bino_ridge_train_AUC", "bino_ridge_test_left_AUC", "bino_ridge_test_balance_AUC", "bino_ridge_misClass")

pheno <- read.xlsx(paste0(data_dir, "combined/Clinicopathologic_v2.xlsx"), na.strings=c("","NA","Nan","#N/A"))
datasets <- c("E-MTAB-4321", "GSE32894", "GSE13507")

# function format: dataset, survival_type, survival_time)
survival_pack("E-MTAB-4321", "Progression_beyond_T2_Progression", "Progression_free_survival")
# survival_pack("GSE32894", "Cancer_specific_satus_DOD_event", "Disease_specific_survival")
# survival_pack("GSE13507", "Vital_status", "OS")

# write to file
write.xlsx(models, paste0(result_dir, "10_Models_v4/10.2_Modes_AUC.xlsx"), overwrite=TRUE)
