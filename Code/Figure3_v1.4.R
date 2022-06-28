# modified according to reviewer1's suggestion: update with "6_Gene_Survival_v2.R"

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
library(ggtree)
library(aplot)
library(dplyr)
library(tidytree)

# set working directory
args <- commandArgs(T)
data_dir <- paste0(args[1], "Data/")
result_dir <- paste0(args[1], "Results/")
# dir.create(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2/"), recursive=TRUE)
dir.create(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/"), recursive=TRUE)
# dir.create(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.4_2.5/"), recursive=TRUE)

###################################  main ##################################
pheno <- read.xlsx(paste0(data_dir, "combined/Clinicopathologic_v2.xlsx"), na.strings=c("","NA","Nan","#N/A"))
gene_survival <- read.xlsx(paste0(result_dir, "6_Gene_Survival/6.1_Gene_Survival.xlsx"))
gene_survival_sig <- gene_survival[gene_survival$KM_Pvalue < 0.05 & gene_survival$Forest_Pvalue < 0.05 & (gene_survival$HR_High < 0.5 | gene_survival$HR_High > 2.5) & gene_survival$HR_High != "not-converge",]
write.xlsx(gene_survival_sig, paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.2_Heatmap_data.xlsx"), overwrite=TRUE)

for(survival_time in c("Progression_free_survival", "OS", "all")){

	################ 0) Pre-processing
	if(survival_time=="Progression_free_survival"){
		survival_type <- "Progression_beyond_T2_Progression"
	}else if(survival_time=="OS"){
		survival_type <- "Vital_status"
	}

	if(survival_time!="all"){
		heatmap_data <- gene_survival_sig[gene_survival_sig$survival_time==survival_time,]
		dataset_name <- unique(heatmap_data$dataset_name)
		candidate_genes <- unique(heatmap_data$gene_name)
	
		# expression data
		exp_data <- read.csv(paste0(result_dir, "1_data_clean/", dataset_name, "_clean.txt"), sep="\t", row.names=1)
		exp_matrix <- exp_data[match(candidate_genes, rownames(exp_data)),]
	
		# survival data
		survival_data <- pheno[pheno$Source_dataset==dataset_name,]
		survival_matrix <- survival_data[match(names(exp_matrix), survival_data$Sample_name),]
		survival_matrix <- survival_matrix[!is.na(survival_matrix[,survival_time]),]
	
		# re-order by sample
		survival_matrix <- survival_matrix[order(-survival_matrix[,survival_type], survival_matrix[,survival_time], decreasing=TRUE),]
		exp_matrix <- exp_matrix[,match(survival_matrix$Sample_name, names(exp_matrix))]
		exp_matrix <- as.matrix(exp_matrix)

		################ 1) Heatmap
		column_ha = HeatmapAnnotation(Status=survival_matrix[,survival_type], Survival=anno_barplot(as.numeric(survival_matrix[,survival_time]), gp=gpar(fill="#6baed6", border=NA, lty="blank")), col=list(Status = c("TRUE"="#d7191c", "FALSE"="#6baed6")))
		row_ha = rowAnnotation(Haze_ratio=anno_barplot(as.numeric(unique(heatmap_data$HR_High)), baseline=1, gp=gpar(fill="#6baed6", border=NA, lty="blank")))
		
		if(survival_time=="Progression_free_survival"){
			group <- as.factor(t(survival_matrix$Progression_beyond_T2_Progression))
			names(group) <- survival_matrix$Sample_name
			heatmap <- Heatmap(exp_matrix, 
				col=colorRamp2(c(10,7,4,2,0), c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")), 
				width=unit(9,"cm"), 
				height=unit(15,"cm"), 
				row_names_gp=gpar(fontsize=4), 
				top_annotation=column_ha, 
				right_annotation=row_ha, 
				show_column_dend=FALSE, 
				show_column_names=FALSE, 
				show_row_dend=FALSE, 
				row_dend_reorder=TRUE, 
				cluster_columns = FALSE, 
				# column_dend_reorder=TRUE, 
				# cluster_columns=cluster_within_group(exp_matrix, group), 
				column_split = group, 
				column_title = NULL)
			
			# for pptx bar plot and text annotation (Figure3A version 1.4)
			column_ha_text = HeatmapAnnotation(Status=survival_matrix[c(1,445),survival_type], Survival=anno_barplot(as.numeric(survival_matrix[c(1,445),survival_time]), gp=gpar(fill="#6baed6", border=NA, lty="blank")), col=list(Status = c("TRUE"="#d7191c", "FALSE"="#6baed6")))
			bar_text <- Heatmap(exp_matrix[,c(1,445)], 
				col=colorRamp2(c(10,7,4,2,0), c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")), 
				width=unit(9,"cm"), 
				height=unit(14,"cm"), 
				row_names_gp=gpar(fontsize=4), 
				# top_annotation=column_ha_text, 
				right_annotation=row_ha, 
				show_column_dend=FALSE, 
				show_column_names=FALSE, 
				show_row_dend=FALSE, 
				row_dend_reorder=TRUE, 
				row_order = row_order(heatmap), 
				cluster_columns = FALSE)

			topptx(bar_text, paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.2_", survival_time, "_DEG_heatmap_bar_text.pptx"))

			pdf(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.2_", survival_time, "_DEG_heatmap_bar_text.pdf"))
			print(bar_text)
			dev.off()
		}else{
			group <- as.factor(t(survival_matrix$Vital_status))
			names(group) <- survival_matrix$Sample_name
			heatmap <- Heatmap(exp_matrix, 
				col=colorRamp2(c(16,12,10,9,7), c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")), 
				width=unit(9,"cm"), 
				height=unit(7,"cm"), 
				row_names_gp=gpar(fontsize=6), 
				top_annotation=column_ha, 
				right_annotation=row_ha, 
				show_column_dend=FALSE, 
				show_column_names=FALSE, 
				show_row_dend=FALSE, 
				row_dend_reorder=TRUE, 
				cluster_columns = FALSE, 
				# column_dend_reorder=TRUE, 
				# cluster_columns=cluster_within_group(exp_matrix, group), 
				column_split = group, 
				column_title = NULL)
		}
	
		pdf(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.2_", survival_time, "_DEG_heatmap.pdf"))
		print(heatmap)
		dev.off()

		topptx(heatmap, paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.2_", survival_time, "_DEG_heatmap.pptx"))

		# ggsave(heatmap, paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.2_", survival_time, "_DEG_heatmap.svg"))
	
	}else{
		heatmap_data <- gene_survival_sig[gene_survival_sig$survival_time!="Disease_specific_survival",]
	}
	
	################ 2) KEGG, GO Enrichment: ALL
	risk_factor <- mapIds(org.Hs.eg.db, keys=unique(heatmap_data[heatmap_data$HR_High > 1,"gene_name"]), keytype="SYMBOL", column="ENTREZID")
	protect_factor <- mapIds(org.Hs.eg.db, keys=unique(heatmap_data[heatmap_data$HR_High < 1,"gene_name"]), keytype="SYMBOL", column="ENTREZID")
	top_genes <- unique(union(risk_factor, protect_factor))

	# Pathway enrichment: KEGG, Reactome, WikiPathway
	top_genes_kegg <- enrichKEGG(gene=top_genes, organism="hsa", pvalueCutoff=0.05)
	risk_kegg <- enrichKEGG(gene=risk_factor, organism="hsa", pvalueCutoff=0.05)
	risk_reactome <- enrichPathway(gene=risk_factor, pvalueCutoff=0.05, readable=TRUE)
	risk_wiki <- enrichWP(gene=risk_factor, organism="Homo sapiens")
	protect_kegg <- enrichKEGG(gene=protect_factor, organism="hsa", pvalueCutoff=0.05)
	protect_reactome <- enrichPathway(gene=protect_factor, pvalueCutoff=0.05, readable=TRUE)
	protect_wiki <- enrichWP(gene=protect_factor, organism="Homo sapiens")

	for(kegg in c("top_genes_kegg", "risk_kegg", "protect_kegg")){
		plot_data <- assign(kegg, get(kegg))
		if(nrow(summary(plot_data)) > 0){
			kegg_dot <- dotplot(plot_data)
			kegg_bar <- barplot(plot_data)

			topptx(kegg_dot, paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.3_", survival_time, "_", kegg, "_dot.pptx"))
			topptx(kegg_bar, paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.3_", survival_time, "_", kegg, "_bar.pptx"))

			jpeg(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.3_", survival_time, "_", kegg, "_dot.jpg"))
			print(kegg_dot)
			dev.off()

			jpeg(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.3_", survival_time, "_", kegg, "_bar.jpg"))
			print(kegg_bar)
			dev.off()

			pdf(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.3_", survival_time, "_", kegg, "_dot.pdf"))
			print(kegg_dot)
			dev.off()

			pdf(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.3_", survival_time, "_", kegg, "_bar.pdf"))
			print(kegg_bar)
			dev.off()

			# delete disease related pathways of top_genes_kegg (Figure3C version v1.4)
			if(kegg=="top_genes_kegg"){
				kegg_dot_select <- dotplot(plot_data, showCategory = c("Focal adhesion", "Cytokine−cytokine receptor interaction", "PI3K-Akt signaling pathway", "Chemokine signaling pathway", "ECM-receptor interaction", "Protein digestion and absorption"))

				pdf(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.3_", survival_time, "_", kegg, "_dot_select.pdf"))
				print(kegg_dot_select + theme(plot.margin=unit(c(6,2,6,2), "cm")))
				dev.off()

				topptx(kegg_dot_select + theme(plot.margin=unit(c(5,2,5,2), "cm")), paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.3_", survival_time, "_", kegg, "_dot_select.pptx"))
			}
		}
	}

	# GO enrichment: CC, MF, BP
	top_genes_go_cc <- enrichGO(gene=top_genes, OrgDb=org.Hs.eg.db, ont="CC", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE)
	top_genes_go_mf <- enrichGO(gene=top_genes, OrgDb=org.Hs.eg.db, ont="MF", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE)
	top_genes_go_bp <- enrichGO(gene=top_genes, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE)
	risk_go_cc <- enrichGO(gene=risk_factor, OrgDb=org.Hs.eg.db, ont="CC", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE)
	risk_go_mf <- enrichGO(gene=risk_factor, OrgDb=org.Hs.eg.db, ont="MF", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE)
	risk_go_bp <- enrichGO(gene=risk_factor, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE)
	protect_go_cc <- enrichGO(gene=protect_factor, OrgDb=org.Hs.eg.db, ont="CC", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE)
	protect_go_mf <- enrichGO(gene=protect_factor, OrgDb=org.Hs.eg.db, ont="MF", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE)
	protect_go_bp <- enrichGO(gene=protect_factor, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE)

	for(go in c("top_genes_go_cc", "top_genes_go_mf", "top_genes_go_bp")){
		plot_data <- assign(go, get(go))
		if(nrow(summary(plot_data)) > 1){
			go_dot <- dotplot(plot_data)
			go_tree <- try(treeplot(pairwise_termsim(plot_data)), silent=TRUE)

			topptx(go_dot, paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.4_", survival_time, "_", go, "_dot.pptx"))

			jpeg(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.4_", survival_time, "_", go, "_dot.jpg"))
			print(go_dot)
			dev.off()

			pdf(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.4_", survival_time, "_", go, "_dot.pdf"))
			print(go_dot)
			dev.off()

			if(class(go_tree)!="try-error"){
				jpeg(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.4_", survival_time, "_", go, "_tree.jpg"))
				print(go_tree)
				dev.off()

				pdf(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.4_", survival_time, "_", go, "_tree.pdf"))
				print(go_tree)
				dev.off()

				# for compressed tree plots (Figure3D,3E,3F version 1.4)
				# find node numbers using following code:
				# go_tree_node <- try(treeplot(pairwise_termsim(plot_data)) + geom_text2(aes(subset=!isTip, label=node)), silent=TRUE)

				if(go=="top_genes_go_cc"){
					topptx(scaleClade(go_tree, node=8, scale=0.3), paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.4_", survival_time, "_", go, "_tree.pptx"))
				}else if(go=="top_genes_go_mf"){
					topptx(scaleClade(go_tree, node=12, scale=0.4), paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.4_", survival_time, "_", go, "_tree.pptx"))
				}else{
					topptx(go_tree, paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.4_", survival_time, "_", go, "_tree.pptx"))
				}
			}
		}
	}

	################ 3) KEGG, GO Enrichment: Cell
	unique_cells <- unique(heatmap_data$cell_type)
	for(cell in unique_cells){
		cell_genes <- mapIds(org.Hs.eg.db, keys=heatmap_data[heatmap_data$cell_type==cell,"gene_name"], keytype="SYMBOL", column="ENTREZID")

		# Enrichment anlysis: KEGG, GO
		kegg <- enrichKEGG(gene=cell_genes, organism="hsa", pvalueCutoff=0.05)

		if(as.numeric(nrow(summary(kegg))) > 1){
			topptx(dotplot(kegg), paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.5_", survival_time, "_", cell, "_kegg.pptx"))

			jpeg(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.5_", survival_time, "_", cell, "_kegg.jpg"))
			print(dotplot(kegg))
			dev.off()

			pdf(paste0(result_dir, "6_Gene_Survival/6.1+_HR_0.5_2.5/6.5_", survival_time, "_", cell, "_kegg.pdf"))
			print(dotplot(kegg))
			dev.off()
		}
	}
}
