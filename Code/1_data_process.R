rm(list=ls())
################################## NMIBC:  ##################################
# load new packages
library(openxlsx)

# set path
result_dir <- "C:/0_xmsun/xmsun/Graduate/20210224_NMIBC/Results/"


################################ Download public dataset ################################
dir.create(paste0(result_dir, "0_data/"))

######## GEO data
download_GEO_data <- function(GEO_ID, save_path) {
	# load packages
	library(GEOquery)

	# download data
	gse_data <- getGEO(GEO_ID, destdir = save_path, getGPL = FALSE)

	# expressiion data
	exp_data <- exprs(gse_data[[1]])

	# phenotype data
	pheno_data <- pData(gse_data[[1]])

	write.csv(exp_data, paste0(save_path, GEO_ID, "_exp.csv"))
	write.csv(pheno_data, paste0(save_path, GEO_ID, "_pheno.csv"))
}

# download_GEO_data("GSE13507", paste0(result_dir, "0_data/"))
# download_GEO_data("GSE88", paste0(result_dir, "0_data/"))
# download_GEO_data("GSE89", paste0(result_dir, "0_data/"))
# download_GEO_data("GSE120736", paste0(result_dir, "0_data/"))
# download_GEO_data("GSE32894", paste0(result_dir, "0_data/")) # Unavailable. Process with suppl/GSE32894_non-normalized_308UCsamples.txt instead ### awk 'BEGIN{FS=OFS=","} {print $1, $2, $4, $6, $8, $10, $12, $14, $16, $18, $20, $22, $24, $26, $28, $30, $32, $34, $36, $38, $40, $42, $44, $46, $48, $50, $52, $54, $56, $58, $60, $62, $64, $66, $68, $70, $72, $74, $76, $78, $80, $82, $84, $86, $88, $90, $92, $94, $96, $98, $100, $102, $104, $106, $108, $110, $112, $114, $116, $118, $120, $122, $124, $126, $128, $130, $132, $134, $136, $138, $140, $142, $144, $146, $148, $150, $152, $154, $156, $158, $160, $162, $164, $166, $168, $170, $172, $174, $176, $178, $180, $182, $184, $186, $188, $190, $192, $194, $196, $198, $200, $202, $204, $206, $208, $210, $212, $214, $216, $218, $220, $222, $224, $226, $228, $230, $232, $234, $236, $238, $240, $242, $244, $246, $248, $250, $252, $254, $256, $258, $260, $262, $264, $266, $268, $270, $272, $274, $276, $278, $280, $282, $284, $286, $288, $290, $292, $294, $296, $298, $300, $302, $304, $306, $308, $310, $312, $314, $316, $318, $320, $322, $324, $326, $328, $330, $332, $334, $336, $338, $340, $342, $344, $346, $348, $350, $352, $354, $356, $358, $360, $362, $364, $366, $368, $370, $372, $374, $376, $378, $380, $382, $384, $386, $388, $390, $392, $394, $396, $398, $400, $402, $404, $406, $408, $410, $412, $414, $416, $418, $420, $422, $424, $426, $428, $430, $432, $434, $436, $438, $440, $442, $444, $446, $448, $450, $452, $454, $456, $458, $460, $462, $464, $466, $468, $470, $472, $474, $476, $478, $480, $482, $484, $486, $488, $490, $492, $494, $496, $498, $500, $502, $504, $506, $508, $510, $512, $514, $516, $518, $520, $522, $524, $526, $528, $530, $532, $534, $536, $538, $540, $542, $544, $546, $548, $550, $552, $554, $556, $558, $560, $562, $564, $566, $568, $570, $572, $574, $576, $578, $580, $582, $584, $586, $588, $590, $592, $594, $596, $598, $600, $602, $604, $606, $608, $610, $612, $614, $616}' GSE32894_exp_pre.csv > GSE32894_exp.csv
# download_GEO_data("GSE32549", paste0(result_dir, "0_data/")) # Part of GSE32894
# download_GEO_data("GSE48075", paste0(result_dir, "0_data/"))
# download_GEO_data("GSE83586", paste0(result_dir, "0_data/"))
# download_GEO_data("GSE128959", paste0(result_dir, "0_data/"))
# download_GEO_data("GSE3167", paste0(result_dir, "0_data/"))
# download_GEO_data("GSE5479", paste0(result_dir, "0_data/")) # (old dataset in 2007) no gene/probe names # clinical info in supplemeantary file
# download_GEO_data("GSE154261", paste0(result_dir, "0_data/")) # exp data need manually work
# download_GEO_data("GSE163209", paste0(result_dir, "0_data/")) # probe-gene transfer manually


######## ArrayExpress data
# # install packages
# a = rownames(installed.packages())
# install_bioc <- c("Biobase","oligoClasses","ArrayExpress", "pd.hugene.1.0.st.v1", 
#                   "hugene10sttranscriptcluster.db", "oligo", "arrayQualityMetrics",
#                   "limma", "topGO", "ReactomePA", "clusterProfiler", "gplots", "ggplot2",
#                   "geneplotter", "RColorBrewer", "pheatmap", "dplyr","stringr","genefilter")
# for(i in install_bioc) {if(! i %in% a) BiocManager::install(i, update=F)}

download_ArrayExpress_data <- function(AE_ID, save_path){
	# load packages
	library(ArrayExpress)

	# download data (set local to FALSE when download new datasets)
	meta_AE <- getAE(AE_ID, path=save_path, type="full", local=TRUE)

	# phenotype data
	pheno_data <- read.delim(paste0(save_path, meta_AE$sdrf), check.names=FALSE)
	write.csv(pheno_data, paste0(save_path, AE_ID, "_pheno.csv"), row.names=FALSE, quote=FALSE)

	# expression data
	if(!is.null(meta_AE$processedFiles)){
		# processed data
		exp_data <- read.delim(paste0(save_path, meta_AE$processedFiles))
		write.csv(exp_data, paste0(save_path, AE_ID, "_exp.csv"), row.names=FALSE, quote=FALSE)
	}else if(!is.null(meta_AE$rawFiles)){
		# raw CEL data
		library(oligo)
		library(affyio)
		library(Biobase)

		raw_data <- oligo::read.celfiles(paste0(save_path, meta_AE$rawFiles))
		stopifnot(validObject(raw_data))
		exp_raw <- oligo::rma(raw_data)
		exp_data <- Biobase::exprs(exp_raw)
		getNetAffx(exp_data)
		write.csv(exp_data, paste0(save_path, AE_ID, "_exp.csv"), quote=FALSE)
	}
}

# download_ArrayExpress_data("E-MTAB-1940", paste0(result_dir, "0_data/"))
# download_ArrayExpress_data("E-MTAB-4321", paste0(result_dir, "0_data/"))
# download_ArrayExpress_data("E-MTAB-1803", paste0(result_dir, "0_data/")) # two exp files; # all MIBC


################################## Data cleaning ##################################
data_cleaning_array <- function(exp_data, dataset_name, bioc_package, mibc_col_name, mibc_value) {
	# # 1) Expression gene annotation
	# # 1.1) acquire gene symbol annotation from bioconductor (ex. GPL6102)
	# # Method 1: use GEOmetadb
	# # BiocManager::install("GEOmetadb")
	# library(GEOmetadb)
	# sqlfile <- getSQLiteFile()
	# con <- dbConnect(SQLite(), sqlfile)
	# GPL6102 <- dbGetQuery(con,"select gpl,title,bioc_package from gpl where gpl='GPL6102'")
	# dim(GPL6102)
	
	# Method 2: directly install from bioconductor 
	# # bioconductor package for chip information
	# # a) GPL_ID to bioc_package: blog.csdn.net/weixin_40739969/article/details/103186027
	# BiocManager::install("illuminaHumanv2.db") # GPL6102
	# BiocManager::install("illuminaHumanv3.db") # GPL6947, GSE48075
	# BiocManager::install("illuminaHumanv4.db") # GPL10558
	# BiocManager::install("hu6800.db") # GPL80
	# BiocManager::install("hgu133plus2.db") # A-AFFY-44
	# BiocManager::install("hugene10sttranscriptcluster.db") # GSE83586,GSE128959
	# BiocManager::install("hgu133a.db") # GSE3167

	# # b) GPL_ID to bioc_package: https://rdrr.io/bioc/GEOmetadb/man/getBiocPlatformMap.html
	# library(GEOmetadb)
	# if( !file.exists("GEOmetadb.sqlite") ) {
	# 	demo_sqlfile <- getSQLiteFile(destdir = getwd(), destfile = "GEOmetadb.sqlite.gz", type = "demo")
	# } else {
	# 	demo_sqlfile <- "GEOmetadb.sqlite"
	# }
	# con <- dbConnect(SQLite(), demo_sqlfile)
	# getBiocPlatformMap(con, bioc=c('GPL4060', 'GPL20301', 'GPL17586'))
	# dbDisconnect(con)


	# annotation with BiocManager database
	library(illuminaHumanv2.db)
	library(illuminaHumanv3.db)
	library(illuminaHumanv4.db)
	library(hu6800.db)
	library(hgu133plus2.db)
	library(hugene10sttranscriptcluster.db)
	library(hgu133a.db)
	library(annotate)
	probe2gene <- toTable(eval(parse(text=paste0(bioc_package, "SYMBOL"))))

	# # Method 3: use AnnoProbe from jmzeng1314
	# # library(devtools)
	# # install_github("jmzeng1314/AnnoProbe")
	# # setwd("D:/xmsun/Resource/GitHub_package/")
	# # devtools::install_local('AnnoProbe-master.zip')
	# library(AnnoProbe)
	# probe2gene_GPL6102_pipe <- idmap("GPL6102", type="pipe")
	# probe2gene_GPL6102_soft <- idmap("GPL6102", type="soft")
	# probe2gene_GPL6102_bioc <- idmap("GPL6102", type="bioc")

	# 1.2) annotation, keep probes with gene_symbol only
	# # To EntrezID
	# all_probe=eval(parse(text = paste0('mappedkeys(','illuminaHumanv2','ENTREZID)')))
	# EGID <- as.numeric(lookUp(all_probe, "illuminaHumanv2.db", "ENTREZID"))
	# To symbol
	exp_valid <- exp_data[rownames(exp_data) %in% probe2gene$probe_id,]
	probe2gene <- probe2gene[match(rownames(exp_valid), probe2gene$probe_id),]
	
	# 2) only keep the probe with the max expression value for the same gene
	tmp <- by(exp_valid, probe2gene$symbol, function(x) rownames(x)[which.max(rowMeans(x))])
	probes <- as.character(tmp)
	dim(exp_valid)
	exp_clean <- exp_valid[rownames(exp_valid) %in% probes,]
	dim(exp_clean)
	
	# 3) transfer probe names to gene names
	rownames(exp_clean) <- probe2gene[match(rownames(exp_clean), probe2gene$probe_id), 2]
	exp_clean <- cbind(rownames(exp_clean), exp_clean)
	names(exp_clean)[1] <- "Symbol"

	# extra log2: suitable for GSE3167, GSE32894, GSE88, GSE89
	if(dataset_name %in% c("GSE32894", "GSE3167", "GSE88", "GSE89")){
		exp_clean[,2:ncol(exp_clean)] <- log2(1 + exp_clean[,2:ncol(exp_clean)])
	}
	
	# write out clean expression dataset
	write.table(exp_clean, paste0(result_dir, "1_data_clean/all_", dataset_name, "_clean.txt"), quote=FALSE, sep="\t", row.names=FALSE)


	############################ NMIBC samples selection ############################
	# 4) delete T2-4 MIBC samples
	if(mibc_value==""){
		exp_clean_nmibc <- exp_clean
		deleted_samples_count[dataset_name, "deleted_samples_count"] <<- 0
	}else{
		pheno_data <- read.csv(paste0(result_dir, "0_data/", dataset_name, "_pheno.csv"), check.names=FALSE)
		mibc_col <- grepl(paste(mibc_value_all, collapse="|"), pheno_data[,names(pheno_data)==mibc_col_name])
		deleted_samples_count[dataset_name, "deleted_samples_count"] <<- length(grep("TRUE", mibc_col))
		
		nmibc_col <- !mibc_col
		nmibc_sample_id <- pheno_data[nmibc_col,1]
		exp_clean_nmibc <- exp_clean[,match(c("Symbol", nmibc_sample_id), names(exp_clean))]
	}

	# write out clean expression dataset (NMIBC only)
	write.table(exp_clean_nmibc, paste0(result_dir, "1_data_clean/", dataset_name, "_clean.txt"), quote=FALSE, sep="\t", row.names=FALSE)
}

data_cleaning_rna_seq <- function(exp_data, dataset_name, mibc_col_name, mibc_value){
	# keep max expression probe for each gene
	symbol_column <- grep("gene.name|symbol", names(exp_data))
	probe2gene <- data.frame(rownames(exp_data), exp_data[,symbol_column])

	exp_valid <- exp_data[,-grep("gene", names(exp_data))]
	tmp <- by(exp_valid, probe2gene[,2], function(x) rownames(x)[which.max(rowMeans(x))])
	probes <- as.character(tmp)
	exp_clean <- exp_valid[rownames(exp_valid) %in% probes,]
	
	# set gene symbols as rowname
	row.names(exp_clean) <- probe2gene[match(rownames(exp_clean), probe2gene[,1]), 2]
	exp_clean <- cbind(rownames(exp_clean), exp_clean)
	names(exp_clean)[1] <- "Symbol"

	# log2: only E-MTAB-4321, PMID31286493	
	if(dataset_name %in% c("E-MTAB-4321", "PMID31286493")){
		exp_clean[,2:ncol(exp_clean)] <- log2(1 + exp_clean[,2:ncol(exp_clean)])
	}

	# write out clean expression dataset
	write.table(exp_clean, paste0(result_dir, "1_data_clean/all_", dataset_name, "_clean.txt"), quote=FALSE, sep="\t", row.names=FALSE)

	############################ NMIBC samples selection ############################
	# 4) delete T2-4 MIBC samples
	if(mibc_value==""){
		exp_clean_nmibc <- exp_clean
		deleted_samples_count[dataset_name, "deleted_samples_count"] <<- 0
	}else{
		pheno_data <- read.csv(paste0(result_dir, "0_data/", dataset_name, "_pheno.csv"), check.names=FALSE)
		mibc_col <- grepl(paste(mibc_value_all, collapse="|"), pheno_data[,names(pheno_data)==mibc_col_name])
		deleted_samples_count[dataset_name, "deleted_samples_count"] <<- length(grep("TRUE", mibc_col))
		
		nmibc_col <- !mibc_col
		nmibc_sample_id <- pheno_data[nmibc_col,1]
		exp_clean_nmibc <- exp_clean[,match(c("Symbol", nmibc_sample_id), names(exp_clean))]
	}
	
	# write out clean expression dataset (NMIBC only)
	write.table(exp_clean_nmibc, paste0(result_dir, "1_data_clean/", dataset_name, "_clean.txt"), quote=FALSE, sep="\t", row.names=FALSE)
}



################################## main ##################################
dir.create(paste0(result_dir, "1_data_clean/"))

platforms <- read.csv(paste0(result_dir, "0_data/platforms.csv"))
datasets_exp <- list.files(paste0(result_dir, "0_data/"))
datasets_exp <- datasets_exp[which(grepl("_exp.csv", datasets_exp))]
datasets_exp <- gsub("_exp.csv", "", datasets_exp)
datasets_yes <- platforms[grep("yes", platforms$included), "dataset"]
datasets <- intersect(datasets_exp, datasets_yes)

# globle varible set
mibc_value_all <- c("Tx", "T2", "T3", "T4", "T2-4", "N/A", "Metastasis")

# deleted MIBC samples in each dataset
deleted_samples_count <- data.frame(matrix(nrow=length(datasets), ncol=1))
names(deleted_samples_count) <- "deleted_samples_count"
rownames(deleted_samples_count) <- datasets

for(dataset in datasets){
	dataset_name <- dataset
	platform_id <- platforms[platforms$dataset==dataset_name,"platform_id"]
	bioc_package <- platforms[platforms$dataset==dataset_name,"bioc_package"]
	tech_info <- platforms[platforms$dataset==dataset_name,"tech"]	
	mibc_col_name <- platforms[platforms$dataset==dataset_name,"mibc_col_name"]
	mibc_value <- platforms[platforms$dataset==dataset_name,"mibc_value"]
	exp_data <- read.csv(paste0(result_dir, "0_data/", dataset, "_exp.csv"), row.names=1)

	if(dataset_name=="PMID31286493"){
		exp_data$gene.name <- rownames(exp_data)
	}

	if(tech_info=="rna-seq"){
		data_cleaning_rna_seq(exp_data, dataset_name, mibc_col_name, mibc_value)
	}else if(tech_info=="array"){
		data_cleaning_array(exp_data, dataset_name, bioc_package, mibc_col_name, mibc_value)
	}
}

write.xlsx(deleted_samples_count, paste0(result_dir, "1_data_clean/deleted_samples_count.xlsx"), row.names=TRUE, overwrite=TRUE)
