rm(list=ls())
################################## BC progression ##################################
# load packages
library(readxl)
library(openxlsx)

# set path
args <- commandArgs(T)
result_dir <- paste0(args[1], "Results/")

################################### Function module ################################## 

####### Funtion: xCell
# devtools::install_github('dviraran/xCell')
run_xcell <- function(dataset_name, exp_data, save_path){
	dir.create(save_path)

	# 1) xCell 
	library(xCell)
	exp_data <- try(xCellAnalysis(exp_data))
	
	# 2) plot heatmap
	library(ComplexHeatmap)
	
	# main plot
	heatmap_data <- Heatmap(exp_data, column_names_gp = gpar(fontsize = 1), row_names_gp = gpar(fontsize = 8))
	pdf(paste0(save_path, dataset_name, "_heatmap.pdf"))
	draw(heatmap_data)
	dev.off()
	
	# cluster sheet
	write.csv(heatmap_data@matrix, paste0(save_path, dataset_name, "_heatmap.csv"))
}

####### Function: immnedeconv
# library(devtools)
# devtools::install_local('GfellerLab-EPIC.zip')
# devtools::install_local('immunedeconv-master.zip')
immunedeconv_pack <- function(dataset_name, exp_data, save_path, timer_cancer, array) {
	library(immunedeconv)
	dir.create(paste0(save_path, dataset_name, "/"))

	##### Function: visualize results as a stacked bar: EPIC (comparable between cell-types)
	stacked_bar_plot <- function(deconv_data, method) {
		library(dplyr)
		library(ggplot2)
		library(tidyr)
		library(immunedeconv)
		library(tibble)
		
		p <- deconv_data %>%
			gather(sample, fraction, -cell_type) %>%
			# plot as stacked bar chart
			ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
			geom_bar(stat='identity') +
			coord_flip() +
			scale_fill_brewer(palette="Paired") +
			scale_x_discrete(limits = rev(levels(deconv_data)))
		pdf(paste0(save_path, dataset_name, "/", method, "_stack.pdf"))
		print(p)
		dev.off()
	}

	##### Function: visulize results by cell-types: MCP-counter (ONLY comparable between sampels)
	cell_type_plot <- function(deconv_data, method) {
		library(dplyr)
		library(ggplot2)
		library(tidyr)
		library(immunedeconv)
		library(tibble)

		p <- deconv_data %>%
			gather(sample, score, -cell_type) %>%
			ggplot(aes(x=sample, y=score, color=cell_type)) +
			geom_point(size=4) +
			facet_wrap(~cell_type, scales="free_x", ncol=3) +
			scale_color_brewer(palette="Paired", guide=FALSE) +
			coord_flip() +
			theme_bw() +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
			theme(axis.text.y = element_text(size = 1))
		pdf(paste0(save_path, dataset_name, "/", method, "_cell_type_plot.pdf"))
		print(p)
		dev.off()
	}

	##### run immunedeconv
	quantiseq <- try(deconvolute(exp_data, "quantiseq", arrays=array))
	mcp_counter <- try(deconvolute(exp_data, "mcp_counter"))
	xcell <- try(deconvolute(exp_data, "xcell"))
	epic <- try(deconvolute(exp_data, "epic"))
	timer <- try(deconvolute(exp_data, "timer", indications=timer_cancer))
	# cibersort <- immunedeconv::deconvolute(exp_data, "cibersort", arrays=array)
	# cibersort_abs <- immunedeconv::deconvolute(exp_data, "cibersort_abs", arrays=array)
	
	knitr::kable(quantiseq, digits=2)

	if(exists("quantiseq") && (!"try-error" %in% class(quantiseq))){
		write.csv(quantiseq, paste0(save_path, dataset_name, "/quantiseq.csv"), row.names=FALSE)
		stacked_bar_plot(quantiseq, "quantiseq")
		cell_type_plot(quantiseq, "quantiseq")
	}
	if(exists("mcp_counter") && (!"try-error" %in% class(mcp_counter))){
		write.csv(mcp_counter, paste0(save_path, dataset_name, "/mcp_counter.csv"), row.names=FALSE)
		stacked_bar_plot(mcp_counter, "mcp_counter")
		cell_type_plot(mcp_counter, "mcp_counter")
	}
	if(exists("xcell") && (!"try-error" %in% class(xcell))){
		write.csv(xcell, paste0(save_path, dataset_name, "/xcell.csv"), row.names=FALSE)
		stacked_bar_plot(xcell, "xcell")
		cell_type_plot(xcell, "xcell")
	}
	if(exists("epic") && (!"try-error" %in% class(epic))){
		write.csv(epic, paste0(save_path, dataset_name, "/epic.csv"), row.names=FALSE)
		stacked_bar_plot(epic, "epic")
		cell_type_plot(epic, "epic")
	}
	if(exists("timer") && (!"try-error" %in% class(timer))){
		write.csv(timer, paste0(save_path, dataset_name, "/timer.csv"), row.names=FALSE)
		stacked_bar_plot(timer, "timer")
		cell_type_plot(timer, "timer")
	}
	
	##### combine socres in one cell type tree
	library(dplyr)
	quantiseq_comb <- deconvolute(immunedeconv::dataset_racle$expr_mat, "quantiseq") %>% map_result_to_celltypes(c("T cell CD4+"), "quantiseq")
	knitr::kable(quantiseq_comb, digits=2)

	rm(quantiseq)
	rm(mcp_counter)
	rm(xcell)
	rm(epic)
	rm(timer)
}

####### Function: ESTIMATE
# library(utils)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos=rforge, dependencies=TRUE)
estimate <- function(dataset_name, exp_data, save_path, platform){
	dir.create(save_path)

	# 1) ESTIMATE
	library(estimate)
	write.table(exp_data, file=paste0(save_path, dataset_name, "_estimate_input.txt"), sep="\t", quote=FALSE)
	filterCommonGenes(input.f=paste0(save_path, dataset_name, "_estimate_input.txt"), output.f=paste0(save_path, dataset_name, "_estimate_gene.gct"), id="GeneSymbol")
	estimateScore(input.ds=paste0(save_path, dataset_name, "_estimate_gene.gct"),  output.ds=paste0(save_path, dataset_name, "_estimate_score.gct"), platform=platform)
	
	# write ESTIMATE scores to excel file
	estimate_score <- read.table(paste0(save_path, dataset_name, "_estimate_score.gct"), skip=2, header=TRUE)
	rownames(estimate_score)=estimate_score[,1]
	estimate_score=estimate_score[,3:ncol(estimate_score)]
	write.xlsx(estimate_score, paste0(save_path, dataset_name, "_estimate_score.xlsx"), rowNames=TRUE, overwrite=TRUE)

	rm(estimate_score)
}



################################### main ################################## 
datasets_all <- list.files(paste0(result_dir, "1_data_clean/"))
datasets_all <- datasets_all[grep(".txt$", datasets_all)]
datasets <- datasets_all[!grepl("^a|GSE88|GSE89|PMID", datasets_all)]
platforms <- read.csv(paste0(result_dir, "0_data/platforms.csv"))

for(dataset in datasets){
	dataset_name <- gsub("_clean.txt", "", dataset)
	tech_info <- platforms[platforms$dataset==dataset_name, "tech"]
	platform <- platforms[platforms$dataset==dataset_name, "platform"]
	exp_data <- read.csv(paste0(result_dir, "1_data_clean/", dataset), sep="\t", row.names=1)

	# ##### 1) xCell
	run_xcell(dataset_name, exp_data, paste0(result_dir, "2_xCell/"))

	# ##### 2) immunedeconv
	if(tech_info=="array"){
		array <- TRUE
	}else if(tech_info=="rna-seq"){
		array <- FALSE
	}
	immunedeconv_pack(dataset_name, exp_data, paste0(result_dir, "2_immunedeconv/"), rep("BLCA", ncol(exp_data)), TRUE)

	##### 3) ESTIMATE
	# For RNA-Seq datasets: set platform as "illumina"
	if(platform=="rna-seq"){
		platform <- "illumina"
	}
	estimate(dataset_name, exp_data, paste0(result_dir, "2_ESTIMATE/"), platform)
}
