rm(list=ls())
################################## BC progression ##################################
# load package
# install.packages("openxlsx", dependencies=TRUE)
# install.packages("ggsignif", dependencies=TRUE)
library(dplyr)
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggpubr)

# set working directory
data_dir <- "C:/0_xmsun/xmsun/Graduate/20210224_NMIBC/Data/"
result_dir <- "C:/0_xmsun/xmsun/Graduate/20210224_NMIBC/Results/"



################################### Function module ################################## 

######### Function: Data matrix for box plot
compare_box_pack <- function(pheno, dataset_name){
	
	pheno_matrix <- pheno[pheno$Source_dataset==dataset_name,]

	for(method in c("xCell", "ESTIMATE", "immunedeconv")){
		if(method == "xCell" && file.exists(paste0(result_dir, "2_", method, "/", dataset_name, "_heatmap.csv"))){
			method_name <- method
			score_matrix <- read.csv(paste0(result_dir, "2_", method, "/", dataset_name, "_heatmap.csv"), row.names=1, check.names=FALSE)
			output <- box_plot(pheno, dataset_name, method_name, score_matrix, pheno_matrix, output)
		}else if(method == "ESTIMATE" && file.exists(paste0(result_dir, "2_", method, "/", dataset_name, "_estimate_score.xlsx"))){
			method_name <- method
			score_matrix <- read.xlsx(paste0(result_dir, "2_", method, "/", dataset_name, "_estimate_score.xlsx"), rowNames=TRUE)
			output <- box_plot(pheno, dataset_name, method_name, score_matrix, pheno_matrix, output)
		}else if(method == "immunedeconv"){
			sub_methods <- list.files(paste0(result_dir, "2_", method, "/", dataset_name))
			sub_methods <- gsub(".csv", "", sub_methods[grep(".csv", sub_methods)])
			for(sub_method in sub_methods){
				method_name <- paste0(method, "_", sub_method)
				score_matrix <- read.csv(paste0(result_dir, "2_", method, "/", dataset_name, "/", sub_method, ".csv"), row.names=1, check.names=FALSE)
				output <- box_plot(pheno, dataset_name, method_name, score_matrix, pheno_matrix, output)
			}
		}
	}
	return(output)
}

######### Function: Box plot
box_plot <- function(pheno, dataset_name, method_name, score_matrix, pheno_matrix, output){
	current_line <- nrow(output)

	score_matrix <- data.frame(t(score_matrix))
	score_matrix$sample_name <- rownames(score_matrix)

	for(group_name in c("Sex", "Age", "Staging", "Grade73", "Grade98", "EORTC_risk_score_NMIBC", "Tumor_size", "Tumor_status", "CIS_in_disease_course", "Recurrence", "Progression_beyond_T2_Progression", "Cancer_specific_satus_DOD_event", "Vital_status")){
		# Results directory for box plots
		dir.create(paste0(result_dir, "3_Cell_Clin_Box/", group_name, "/"), recursive = TRUE)
    print(paste("Start session:", dataset_name, "-", method_name, "-", group_name))
    
		# box data
		box_data <- merge(score_matrix, pheno_matrix, by.x="sample_name", by.y="Sample_name")
		box_data <- box_data%>%filter(is.na(get(group_name))==FALSE, get(group_name)!="")
		if(nrow(box_data)!=0){

			rownames(box_data) <- box_data$sample_name
			box_data <- box_data[,-1]
			
			# relevel factors
			box_data$Sex <- factor(box_data$Sex, levels=c("Female", "Male"))
			box_data$Age <- factor(box_data$Age, levels=c("<=30", "31-50", "51-60", "61-70", "71-80", ">80"))
			box_data$Staging <- factor(box_data$Staging, levels=c("T0", "Ta", "T1", "Ta-1", "T2", "T3", "T4", "T2-4", "Tx"))
			box_data$Grade73 <- factor(box_data$Grade73, levels=c("G0", "G1", "G2", "G3", "Gx"))
			box_data$Grade98 <- factor(box_data$Grade98, levels=c("Low", "High"))
			box_data$EORTC_risk_score_NMIBC <- factor(box_data$EORTC_risk_score_NMIBC, levels=c("0", "1"))
			box_data$Tumor_size <- factor(box_data$Tumor_size, levels=c("<=3cm", ">3cm"))
			box_data$Tumor_status <- factor(box_data$Tumor_status, levels=c("Primary", "Recurrent"))
			box_data$CIS_in_disease_course <- factor(box_data$CIS_in_disease_course, levels=c("CIS-", "CIS+"))
			box_data$Recurrence <- factor(box_data$Recurrence, levels=c("FALSE", "TRUE"))
			box_data$Progression_beyond_T2_Progression <- factor(box_data$Progression_beyond_T2_Progression, levels=c("FALSE", "TRUE"))
			box_data$Cancer_specific_satus_DOD_event <- factor(box_data$Cancer_specific_satus_DOD_event, levels=c("Survival", "Death"))
			box_data$Vital_status <- factor(box_data$Vital_status, levels=c("Survival", "Death"))

			# t test between groups
			end_col <- grep("Source_dataset", names(box_data))-1
			for(i in 1:end_col){
				col_name <- names(box_data)[i]
				temp_data <- box_data[,c(col_name, group_name)]
				names(temp_data) <- c("Scores", "Group")
				temp_data <- temp_data %>% group_by(Group) %>% filter(n() >= 3) # delete groups with less than 3 samples

				levels <- length(unique(temp_data$Group))
				if(levels==1){
					output[current_line+1, "Source_dataset"] <- dataset_name
					output[current_line+1, "Source_method"] <- method_name
					output[current_line+1, "Cell_name"] <- col_name
					output[current_line+1, "Group_name"] <- group_name
				}else{
					# t test
					compare <- compare_means(Scores ~  Group, data=temp_data, method="t.test")
					my_comparisons <- split(unlist(sapply(data.frame(compare[,c("group1", "group2")]),c)), seq(nrow(compare)))
	
					p <- ggplot(temp_data, aes(Group, Scores, fill=Group)) +
						geom_boxplot(width=0.5) +
						stat_compare_means(comparisons=my_comparisons, method="t.test", size=8) +
						labs(x=group_name, y=col_name) +
						theme(plot.title=element_text(size = 25), 
							axis.text.x=element_text(size=15,angle=0), 
							axis.text.y=element_text(size=15), 
							axis.title.x=element_text(size = 23), 
							axis.title.y=element_text(size = 23))
	
					jpeg(paste0(result_dir, "3_Cell_Clin_Box/", group_name, "/", dataset_name, "_", method_name, "_", col_name, ".jpg"))
					print(p)
					dev.off()

					n <- nrow(compare)
					means <- aggregate(temp_data$Scores, by=list(temp_data$Group), mean)
					names(means) <- c("Group", "Mean")
	
					if(n==0){
						output[current_line+1, "Source_dataset"] <- dataset_name
						output[current_line+1, "Source_method"] <- method_name
						output[current_line+1, "Cell_name"] <- col_name
						output[current_line+1, "Group_name"] <- group_name
					}else{
						output[(current_line+1):(current_line+n), "Source_dataset"] <- dataset_name
						output[(current_line+1):(current_line+n), "Source_method"] <- method_name
						output[(current_line+1):(current_line+n), "Cell_name"] <- col_name
						output[(current_line+1):(current_line+n), "Group_name"] <- group_name
						for(j in 1:n){
							output[current_line+j, "groupA_name"] <- compare$group1[j]
							output[current_line+j, "groupB_name"] <- compare$group2[j]
							output[current_line+j, "Pvalue"] <- compare$p[j]
							output[current_line+j, "Padjust"] <- compare$p.adj[j]
							output[current_line+j, "groupA_mean"] <- means[means$Group==compare$group1[j],"Mean"]
							output[current_line+j, "groupB_mean"] <- means[means$Group==compare$group2[j],"Mean"]
						}
					}
				}
				current_line <- nrow(output)
			}
		}
	}
	return(output)
}



######################################## main ########################################
pheno <- read.xlsx(paste0(data_dir, "combined/Clinicopathologic_v2.xlsx"), na.strings=c("","NA","Nan","#N/A"))
datasets_all <- list.files(paste0(result_dir, "1_data_clean/"))
datasets_all <- datasets_all[grep(".txt$", datasets_all)]
datasets <- datasets_all[!grepl("^a|GSE88|GSE89|PMID", datasets_all)]

# output format
output <- data.frame(matrix(nrow=0, ncol=10))
names(output) <- c("Source_dataset", "Source_method", "Cell_name", "Group_name", "groupA_name", "groupA_mean", "groupB_name", "groupB_mean", "Pvalue", "Padjust")

for(dataset in datasets){
	dataset_name <- gsub("_clean.txt", "", dataset)
	output <- compare_box_pack(pheno, dataset_name)
}

write.xlsx(output, paste0(result_dir, "3_Cell_Clin_Box/3.1_Compare.xlsx"), overwrite=TRUE)

########## select candidate cell types
cell_types <- read.xlsx(paste0(result_dir, "3_Cell_Clin_Box/3.0_Cell_types.xlsx"))
cell_stat_P <- cell_types
cell_stat_Padj <- cell_types
cell_stat_P$stat <- ""
cell_stat_Padj$stat <- ""
i = 1
j = 3
for(i in 1:nrow(cell_types)){
	for(j in 3:ncol(cell_types)){
		if(!is.na(cell_types[i,j])){
			method_name <- names(cell_types)[j]
			cell_name <- cell_types[i,j]
			output_refine <- output[output$Cell_name==cell_name & output$Source_method==method_name, ]
			cell_stat_P[i,j] <- length(grep("TRUE", output_refine$Pvalue<0.05))/nrow(na.omit(output_refine))
			cell_stat_Padj[i,j] <- length(grep("TRUE", output_refine$Padjust<0.05))/nrow(na.omit(output_refine))
		}
	}
}

cell_stat_P$stat <- rowMeans(apply(cell_stat_P[,3:ncol(cell_types)], 2, as.numeric), na.rm=TRUE)
cell_stat_Padj$stat <- rowMeans(apply(cell_stat_Padj[,3:ncol(cell_types)], 2, as.numeric), na.rm=TRUE)

write.xlsx(cell_stat_P, paste0(result_dir, "3_Cell_Clin_Box/3.1_Cell_stat_P.xlsx"), overwrite=TRUE, showNA=FALSE, na.string="")
write.xlsx(cell_stat_Padj, paste0(result_dir, "3_Cell_Clin_Box/3.1_Cell_stat_Padj.xlsx"), overwrite=TRUE, keepNA=FALSE, na.string="")
