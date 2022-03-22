rm(list=ls())
################################## BC progression ##################################
# load package
library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(eoffice)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggnewscale)
library(cowplot)
library(enrichplot)
library(DOSE)
library(GOSemSim)
library(ReactomePA)
library(RColorBrewer)

# set working directory
args <- commandArgs(T)
data_dir <- paste0(args[1], "Data/")
result_dir <- paste0(args[1], "Results/")


################################### Function module ##################################
# phenotype
pheno <- read.xlsx(paste0(data_dir, "combined/Clinicopathologic_v2.xlsx"), na.strings=c("","NA","Nan","#N/A"))
datasets_all <- list.files(paste0(result_dir, "1_data_clean/"))
datasets_all <- datasets_all[grep(".txt$", datasets_all)]
datasets <- datasets_all[!grepl("^a|GSE88|GSE89|PMID", datasets_all)]

# 1) get combined cell score matirx
# output 1: cell score matrix
ssgsea_matrix <- data.frame(matrix(nrow=9, ncol=0))
cell_z_matrix <- data.frame(matrix(nrow=9, ncol=0))

for(i in 1:length(datasets)){
	dataset_name <- gsub("_clean.txt", "", datasets[i])
	ssgsea_score <- read.xlsx(paste0(result_dir, "8_Cell_scores/8.3_ssGSEA_heatmap/", dataset_name, "_ssGSEA_score.xlsx"), rowNames=TRUE)
	cell_z <- read.xlsx(paste0(result_dir, "8_Cell_scores/8.4_cell_z-score_heatmap/", dataset_name, "_cell_z-score.xlsx"), rowNames=TRUE)

	if(i==1){
		ssgsea_matrix <- ssgsea_score
		cell_z_matrix <- cell_z
	}else{
		ssgsea_matrix <- merge(ssgsea_matrix, ssgsea_score, by="row.names")
		cell_z_matrix <- merge(cell_z_matrix, cell_z, by="row.names")
		rownames(ssgsea_matrix) <- ssgsea_matrix$Row.names
		ssgsea_matrix <- ssgsea_matrix[,-grep("Row.names", names(ssgsea_matrix))]
		rownames(cell_z_matrix) <- cell_z_matrix$Row.names
		cell_z_matrix <- cell_z_matrix[,-grep("Row.names", names(cell_z_matrix))]
	}
}

write.xlsx(ssgsea_matrix, paste0(result_dir, "8_Cell_scores/8.3_ssGSEA_matrix.xlsx"), overwrite=TRUE, rowNames=TRUE)
write.xlsx(cell_z_matrix, paste0(result_dir, "8_Cell_scores/8.4_cell_z-score_matrix.xlsx"), overwrite=TRUE, rowNames=TRUE)

# 2) Heatmap for immune cells
for(method in c("ssgsea", "cell_z")){
	matrix_name <- paste0(method, "_matrix")
	plot_data <- assign(matrix_name, get(matrix_name))

	# annotation for immune cells
	gene_survival_sig <- read.xlsx(paste0(result_dir, "6_Gene_Survival/6.2_Heatmap_data.xlsx"))
	cell_gene_stat <- table(gene_survival_sig$cell_type)
	row_ha = rowAnnotation(Gene_count=anno_barplot(as.numeric(cell_gene_stat), gp=gpar(fill="#6baed6", border=NA, lty="blank")))

	# annotation for samples
	pheno_all <- pheno[match(names(plot_data), pheno$Sample_name), ]
	number_cols <- c("Age_ori", "Progression_free_survival", "Disease_specific_survival", "OS")
	pheno_all[number_cols] <- lapply(pheno_all[number_cols], as.numeric)
	pheno_all$Staging <- factor(pheno_all$Staging, levels=c("T0", "Ta", "Ta-1", "T1", "Tis"))
	pheno_all$Grade98 <- factor(pheno_all$Grade98, levels=c("Low", "High"))
	pheno_all$Tumor_status <- factor(pheno_all$Tumor_status, levels=c("Primary", "Recurrent", "Progress"))

	column_color <- list(
		# pred=colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c(brewer.pal(5, "Purples"))), 
		Source_dataset=c("E-MTAB-1940"="#8c510a","E-MTAB-4321"="#bf812d","GSE120736"="#dfc27d","GSE128959"="#f6e8c3","GSE13507"="#f5f5f5","GSE3167"="#c7eae5","GSE32894"="#80cdc1","GSE48075"="#35978f","GSE83586"="#01665e"),
		Age=colorRamp2(c(0,25,50,75,100), c(brewer.pal(5,"Blues"))),
		Sex=c("Female"="#d7191c", "Male"="#6baed6"),
		Staging=c("T0"="#fef0d9", "Ta"="#fdcc8a", "Ta-1"="#fc8d59", "T1"="#e34a33", "Tis"="#b30000"),
		Grade73=c("G0"="#feebe2", "G1"="#fbb4b9", "G2"="#f768a1", "G3"="#ae017e", "Gx"="grey"),
		Grade98=c("Low"="#c994c7", "High"="#dd1c77"),
		Tumor_size=c("<=3cm"="#a6bddb", ">3cm"="#1c9099"),
		Tumor_status=c("Primary"="#fde0dd", "Recurrent"="#fa9fb5", "Progress"="#c51b8a"),
		CIS_in_disease_course=c("CIS-"="#fee0d2", "CIS+"="#de2d26"),
		Progression_beyond_T2_Progression=c("FALSE"="#6baed6", "TRUE"="#d7191c"),
		Cancer_specific_satus_DOD_event=c("FALSE"="#6baed6", "TRUE"="#d7191c"),
		Vital_status=c("FALSE"="#6baed6", "TRUE"="#d7191c"))
	column_ha = HeatmapAnnotation(
		# pred=pheno_all$pred, 
		Source_dataset=pheno_all$Source_dataset,
		Age=pheno_all$Age_ori,
		Sex=pheno_all$Sex,
		Staging=pheno_all$Staging,
		Grade73=pheno_all$Grade73,
		Grade98=pheno_all$Grade98,
		Tumor_size=pheno_all$Tumor_size,
		Tumor_status=pheno_all$Tumor_status,
		CIS_in_disease_course=pheno_all$CIS_in_disease_course,
		Progression_beyond_T2_Progression=pheno_all$Progression_beyond_T2_Progression,
		Cancer_specific_satus_DOD_event=pheno_all$Cancer_specific_satus_DOD_event,
		Vital_status=pheno_all$Vital_status,
		col=column_color,
		show_legend=FALSE,
		gap=unit(0.5,"mm"),
		simple_anno_size=unit(3,"mm"),
		annotation_name_gp=gpar(fontsize=8))
	column_ha_legend = HeatmapAnnotation(
		# pred=pheno_all$pred, 
		Source_dataset=pheno_all$Source_dataset,
		Age=pheno_all$Age_ori,
		Sex=pheno_all$Sex,
		Staging=pheno_all$Staging,
		Grade73=pheno_all$Grade73,
		Grade98=pheno_all$Grade98,
		Tumor_size=pheno_all$Tumor_size,
		Tumor_status=pheno_all$Tumor_status,
		CIS_in_disease_course=pheno_all$CIS_in_disease_course,
		Progression_beyond_T2_Progression=pheno_all$Progression_beyond_T2_Progression,
		Cancer_specific_satus_DOD_event=pheno_all$Cancer_specific_satus_DOD_event,
		Vital_status=pheno_all$Vital_status,
		col=column_color,
		show_legend=TRUE,
		gap=unit(0.5,"mm"),
		simple_anno_size=unit(3,"mm"),
		annotation_name_gp=gpar(fontsize=8))

	if(method=="ssgsea"){
		group <- as.factor(pheno_all$Progression_beyond_T2_Progression)
		names(group) <- pheno_all$Sample_name
		heatmap <- Heatmap(as.matrix(plot_data), 
			col=colorRamp2(c(0.6,0.4,0.2,0,-0.2), c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")), 
			heatmap_height=unit(8,"cm"), 
			top_annotation=column_ha, 
			right_annotation=row_ha, 
			show_column_dend=FALSE, 
			show_column_names=FALSE, 
			show_row_dend=FALSE, 
			row_dend_reorder=TRUE,
			column_split=group, 
			column_title=NULL)
		heatmap_legend <- Heatmap(as.matrix(plot_data), 
			col=colorRamp2(c(0.6,0.4,0.2,0,-0.2), c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")), 
			heatmap_height=unit(15,"cm"), 
			top_annotation=column_ha_legend, 
			right_annotation=row_ha, 
			show_column_dend=FALSE, 
			show_column_names=FALSE, 
			show_row_dend=FALSE)

		pdf(paste0(result_dir, "8_Cell_scores/8.3_ssGSEA_heatmap.pdf"))
		print(heatmap)
		dev.off()

		pdf(paste0(result_dir, "8_Cell_scores/8.3_ssGSEA_heatmap_legend.pdf"))
		print(heatmap_legend)
		dev.off()

		jpeg(paste0(result_dir, "8_Cell_scores/8.3_ssGSEA_heatmap.jpg"))
		print(heatmap)
		dev.off()

		jpeg(paste0(result_dir, "8_Cell_scores/8.3_ssGSEA_heatmap_legend.jpg"))
		print(heatmap_legend)
		dev.off()

		# for pptx text annotation (Figure version 1.3)
		legend_only <- packLegend(max_height = unit(14, "cm"), list = list(
			Legend(title="Cell_score", col_fun=colorRamp2(c(0.6,0.4,0.2,0,-0.2), c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba"))), 
			Legend(title="Source_dataset", at=c("E-MTAB-1940", "E-MTAB-4321", "GSE120736", "GSE128959", "GSE13507", "GSE3167", "GSE32894", "GSE48075", "GSE83586"), legend_gp=gpar(fill=c("#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#f5f5f5", "#c7eae5", "#80cdc1", "#35978f", "#01665e"))), 
			Legend(title="Age", col_fun=colorRamp2(c(0,25,50,75,100), c(brewer.pal(5,"Blues")))), 
			Legend(title="Sex", at=c("Female", "Male"), legend_gp=gpar(fill=c("#d7191c", "#6baed6"))), 
			Legend(title="Staging", at=c("T0", "Ta", "Ta-1", "T1", "Tis"), legend_gp=gpar(fill=c("#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"))), 
			Legend(title="Grade73", at=c("G0", "G1", "G2", "G3", "Gx"), legend_gp=gpar(fill=c("#feebe2", "#fbb4b9", "#f768a1", "#ae017e", "grey"))), 
			Legend(title="Grade98", at=c("Low", "High"), legend_gp=gpar(fill=c("#c994c7", "#dd1c77"))), 
			Legend(title="Tumor_size", at=c("<=3cm", ">3cm"), legend_gp=gpar(fill=c("#a6bddb", "#1c9099"))), 
			Legend(title="Tumor_status", at=c("Primary", "Recurrent", "Progress"), legend_gp=gpar(fill=c("#fde0dd", "#fa9fb5", "#c51b8a"))), 
			Legend(title="CIS_in_disease_course", at=c("CIS-", "CIS+"), legend_gp=gpar(fill=c("#fee0d2", "#de2d26"))), 
			Legend(title="Progression_beyond_T2_Progression", at=c("FALSE", "TRUE"), legend_gp=gpar(fill=c("#6baed6", "#d7191c"))), 
			Legend(title="Cancer_specific_satus_DOD_event", at=c("FALSE", "TRUE"), legend_gp=gpar(fill=c("#6baed6", "#d7191c"))), 
			Legend(title="Vital_status", at=c("FALSE", "TRUE"), legend_gp=gpar(fill=c("#6baed6", "#d7191c")))))

		topptx(draw(legend_only), paste0(result_dir, "8_Cell_scores/8.3_ssGSEA_heatmap_legend.pptx"))

		column_ha_text = HeatmapAnnotation(
			# pred=pheno_all$pred, 
			Source_dataset=pheno_all$Source_dataset[1:2],
			Age=pheno_all$Age_ori[1:2],
			Sex=pheno_all$Sex[1:2],
			Staging=pheno_all$Staging[1:2],
			Grade73=pheno_all$Grade73[1:2],
			Grade98=pheno_all$Grade98[1:2],
			Tumor_size=pheno_all$Tumor_size[1:2],
			Tumor_status=pheno_all$Tumor_status[1:2],
			CIS_in_disease_course=pheno_all$CIS_in_disease_course[1:2],
			Progression_beyond_T2_Progression=pheno_all$Progression_beyond_T2_Progression[1:2],
			Cancer_specific_satus_DOD_event=pheno_all$Cancer_specific_satus_DOD_event[1:2],
			Vital_status=pheno_all$Vital_status[1:2],
			col=column_color,
			show_legend=FALSE,
			gap=unit(0.5,"mm"),
			simple_anno_size=unit(3,"mm"),
			annotation_name_gp=gpar(fontsize=8))
		bar_text <- Heatmap(as.matrix(plot_data[,1:2]), 
			col=colorRamp2(c(0.6,0.4,0.2,0,-0.2), c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")), 
			heatmap_height=unit(8,"cm"), 
			top_annotation=column_ha_text, 
			right_annotation=row_ha, 
			show_column_dend=FALSE, 
			show_column_names=FALSE, 
			show_row_dend=FALSE, 
			row_dend_reorder=TRUE, 
			row_order=row_order(heatmap))

		topptx(bar_text, paste0(result_dir, "8_Cell_scores/8.3_ssGSEA_heatmap_bar_text.pptx"))
	}else{
		heatmap <- Heatmap(as.matrix(plot_data), col=colorRamp2(c(100,25,0,-25,-100), c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")), top_annotation=column_ha, right_annotation=row_ha, show_column_dend=FALSE, show_column_names=FALSE, show_row_dend=FALSE, heatmap_height=unit(15,"cm"))
		heatmap_legend <- Heatmap(as.matrix(plot_data), col=colorRamp2(c(100,25,0,-25,-100), c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")), top_annotation=column_ha_legend, right_annotation=row_ha, show_column_dend=FALSE, show_column_names=FALSE, show_row_dend=FALSE, heatmap_height=unit(15,"cm"))

		pdf(paste0(result_dir, "8_Cell_scores/8.4_cell_z-score_heatmap.pdf"))
		print(heatmap)
		dev.off()

		pdf(paste0(result_dir, "8_Cell_scores/8.4_cell_z-score_heatmap_legend.pdf"))
		print(heatmap_legend)
		dev.off()

		topptx(heatmap, paste0(result_dir, "8_Cell_scores/8.4_cell_z-score_heatmap.pptx"))

		topptx(heatmap_legend, paste0(result_dir, "8_Cell_scores/8.4_cell_z-score_heatmap_legend.pptx"))
	}
}
