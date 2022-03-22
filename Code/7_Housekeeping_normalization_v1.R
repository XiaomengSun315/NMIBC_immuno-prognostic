rm(list=ls())
################################## BC progression ##################################
# install.packages("BiocManager")
# BiocManager::install(update=TRUE, ask=FALSE)
# BiocManager::install(version='3.13')
# BiocManager::install("impute", force=TRUE)
library(openxlsx)
library(preprocessCore)
library(impute)


# set working directory
args <- commandArgs(T)
data_dir <- paste0(args[1], "Data/")
result_dir <- paste0(args[1], "Results/")
dir.create(paste0(result_dir, "7_Housekeeping_normalization/"))


######################################## main 1: housekeeping normalization 
# phenotype information
pheno <- read.xlsx(paste0(data_dir, "combined/Clinicopathologic_v2.xlsx"), na.strings=c("","NA","Nan","#N/A"))

# output file 1: housekeeping genes count for each dataset
housekeeping_record <- data.frame(matrix(nrow=0, ncol=3))
names(housekeeping_record) <- c("dataset_name", "housekeeping_count", "housekeeping_genes")

# output file 2: NMIBCexp all
exp_data_all <- data.frame(matrix(nrow=0, ncol=0))

# housekeeping genes
housekeeping_genes <- read.xlsx(paste0(result_dir, "7_Housekeeping_normalization/7.0_Housekeeping_genes.xlsx"))
my_housekeeping <- housekeeping_genes[housekeeping_genes$inBladder_BladderCancer=="yes", "Gene_name"]
GEP_housekeeping <- housekeeping_genes[housekeeping_genes$Source=="GEP", "Gene_name"]

# datasets
datasets_all <- list.files(paste0(result_dir, "1_data_clean/"))
datasets_all <- datasets_all[grep(".txt$", datasets_all)]
datasets <- datasets_all[!grepl("^a|GSE88|GSE89|PMID", datasets_all)]

for(i in 1:length(datasets)){
	dataset <- datasets[i]
	dataset_name <- gsub("_clean.txt", "", dataset)
	exp_data_log <- read.csv(paste0(result_dir, "1_data_clean/", dataset), sep="\t", row.names=1)
	exp_data <- 2^exp_data_log - 1
	exp_data_house <- na.omit(exp_data[match(my_housekeeping, rownames(exp_data)),])
	
	# housekeeping record
	housekeeping_count <- nrow(exp_data_house)
	housekeeping_record[i,"dataset_name"] <- dataset_name
	housekeeping_record[i,"housekeeping_count"] <- housekeeping_count
	housekeeping_record[i,"housekeeping_genes"] <- paste(rownames(exp_data_house),collapse=",")

	# housekeeping normalization: by substracting log2(mean(housekeeping_genes))
	exp_data_house[(housekeeping_count+1),] <- apply(exp_data_house, 2, mean)
	rownames(exp_data_house)[housekeeping_count+1] <- "mean"
	exp_data_norm <- exp_data
	col = 1
	for(col in 1:ncol(exp_data)){
		exp_data_norm[,col] <- exp_data[,col] / exp_data_house["mean",col] * 1000
	}

	write.xlsx(exp_data_norm, paste0(result_dir, "7_Housekeeping_normalization/", dataset_name, "_house_norm.xlsx"), rowNames=TRUE, overwrite=TRUE)

	if(i == 1){
		exp_data_all <- exp_data_norm
	}else{
		exp_data_all$name <- rownames(exp_data_all)
		exp_data_norm$name <- rownames(exp_data_norm)
		exp_data_all_temp <- merge(exp_data_all, exp_data_norm, by="name", all=TRUE, sort=TRUE)
		rownames(exp_data_all_temp) <- exp_data_all_temp$name
		exp_data_all <- exp_data_all_temp
		exp_data_all <- exp_data_all[,-grep("name", names(exp_data_all))]
	}
}

write.xlsx(housekeeping_record, paste0(result_dir, "7_Housekeeping_normalization/7.1_Housekeeping_records.xlsx"), overwrite=TRUE)
write.xlsx(exp_data_all, paste0(result_dir, "7_Housekeeping_normalization/7.1_ALL_house.xlsx"), rowNames=TRUE, overwrite=TRUE)

# NA values statistics
na_stat <- data.frame(matrix(nrow=nrow(exp_data_all), ncol=2))
names(na_stat) <- c("NA_values_count", "NA_values_ratio")
na_stat$NA_values_count <- rowSums(is.na(exp_data_all))
na_stat$NA_values_ratio <- na_stat$NA_values_count/ncol(exp_data_all)
rownames(na_stat) <- rownames(exp_data_all)
na_stat <- na_stat[order(na_stat$NA_values_ratio),]
write.xlsx(na_stat, paste0(result_dir, "7_Housekeeping_normalization/7.1_Na_values_stat.xlsx"), rowNames=TRUE, overwrite=TRUE)


######################################## main 2: missing value imputation (KNN)
# 2.1) delete columns with more than 80% missing values (n=78): PMID31286493 (n=16), GSE88 (n=31), GSE89 (n=31)
exp_data_enough <- exp_data_all[, which(colMeans(!is.na(exp_data_all)) > 0.2)]

# 2.2) KNN impute
exp_data_all_knn <- impute.knn(as.matrix(exp_data_enough), k=10, rowmax=0.5, colmax=0.8, maxp=1500, rng.seed=362436069) # 23846 rows with more than 50 % entries missing

write.csv(exp_data_all_knn$data, paste0(result_dir, "7_Housekeeping_normalization/7.2_ALL_house_knn.csv"))


######################################## main 3: quantile normalization & log2
exp_data_all_knn_quant <- normalize.quantiles(x=as.matrix(exp_data_all_knn$data))
rownames(exp_data_all_knn_quant) <- rownames(exp_data_all_knn$data)
colnames(exp_data_all_knn_quant) <- colnames(exp_data_all_knn$data)

# log 2
exp_data_all_knn_quant_log <- log2(1 + exp_data_all_knn_quant)

write.csv(exp_data_all_knn_quant_log, paste0(result_dir, "7_Housekeeping_normalization/7.3_ALL_house_knn_quantile_log.csv"))

save(exp_data_enough, exp_data_all_knn, exp_data_all_knn_quant, exp_data_all_knn_quant_log, file=paste0(result_dir, "7_Housekeeping_normalization/7.99_ALL_norm_process.Rdata"))
