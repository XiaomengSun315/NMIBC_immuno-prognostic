rm(list=ls())
################################## BC progression ##################################
# load package
# install.packages("openxlsx", dependencies=TRUE)
# install.packages("ggsignif", dependencies=TRUE)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(eoffice)

# set working directory
args <- commandArgs(T)
data_dir <- paste0(args[1], "Data/")
result_dir <- paste0(args[1], "Results/")


################################### Function module ################################## 
cell_clin <- read.xlsx(paste0(result_dir, "3_Cell_Clin_Box/3.1_Compare.xlsx"))
cell_clin <- cell_clin[!cell_clin$Group_name %in% c("Sex", "Age"),]

# 1) refine results
uniq_comb <- unique(cell_clin[,c("Source_dataset", "Source_method", "Cell_name", "Group_name")])
uniq_comb$Pvalue <- ""
uniq_comb$Padjust <- ""
for(i in 1:nrow(uniq_comb)){
	hit_data <- cell_clin[cell_clin$Source_dataset==uniq_comb[i,"Source_dataset"] & cell_clin$Source_method == uniq_comb[i,"Source_method"] & cell_clin$Cell_name == uniq_comb[i,"Cell_name"] & cell_clin$Group_name == uniq_comb[i,"Group_name"],]
	uniq_comb[i,"Pvalue"] <- mean(hit_data$Pvalue)
	uniq_comb[i,"Padjust"] <- mean(hit_data$Padjust)
}

uniq_comb_sig <- na.omit(uniq_comb[uniq_comb$Padjust < 0.05,])
write.xlsx(uniq_comb_sig, paste0(result_dir, "3_Cell_Clin_Box/3.2_uniq_comb_sig.xlsx"), overwrite=TRUE)

# 2) statistics for each cell type
cell_type <- read.xlsx(paste0(result_dir, "3_Cell_Clin_Box/3.0_Cell_types.xlsx"))
cell_stat_vote <- cell_type
cell_stat_FDRp <- cell_type
row = 1
col = 3
for(row in 1:nrow(cell_type)){
	for(col in 3:ncol(cell_type)){
		if(!is.na(cell_type[row, col])){
			cell_sig_data <- na.omit(uniq_comb_sig[uniq_comb_sig$Source_method==names(cell_type)[col] & uniq_comb_sig$Cell_name==cell_type[row,col],])
			cell_stat_vote[row,col] <- as.numeric(nrow(cell_sig_data))
			cell_stat_FDRp[row,col] <- mean(as.numeric(cell_sig_data$Padjust))
		}else{
			cell_stat_vote[row,col] <- 0
			cell_stat_FDRp[row,col] <- ""
		}
	}
}

write.xlsx(cell_stat_vote, paste0(result_dir, "3_Cell_Clin_Box/3.2_Cell_stat_vote.xlsx"), overwrite=TRUE)
write.xlsx(cell_stat_FDRp, paste0(result_dir, "3_Cell_Clin_Box/3.2_Cell_stat_FDRp.xlsx"), overwrite=TRUE)

# 3) Figure 2A
plot_data <- data.frame(rep(names(cell_type[,c(3,5,6,7,8,9)]), times=rep(nrow(cell_type),6)), rep(cell_type$Cell_type_short_name, 6), as.numeric(apply(cell_stat_vote[,c(3,5,6,7,8,9)], 2, unlist)), as.numeric(apply(cell_stat_FDRp[,c(3,5,6,7,8,9)], 2, unlist)))
names(plot_data) <- c("Source_method", "Cell_name", "vote", "FDRp")
p <- ggplot(plot_data, aes(x=Source_method, y=Cell_name, size=vote, color=FDRp)) + geom_point(alpha=0.7)
pdf(paste0(result_dir, "3_Cell_Clin_Box/3.3_Cell_stat.pdf"))
print(p)
dev.off()


plot_data_sig <- na.omit(plot_data)
plot_data_sig <- plot_data_sig[plot_data_sig$Cell_name %in% c("CD4..T.cells", "CD8..T.cells", "B.cells", "Fibroblasts", "NK.cells", "Endothelial.cells", "DC", "Macrophages", "Neutrophils"),]
p_sig <- ggplot(plot_data_sig, aes(x=Cell_name, y=Source_method, size=vote, color=FDRp)) +
	geom_point(alpha=0.6, shape=16) +
	scale_size(range=c(1, 10)) +
	scale_colour_distiller(palette="Blues", limits=c(0, 0.05)) +
	theme_bw() +
	theme(panel.grid=element_blank(), 
		panel.border=element_blank(), 
		axis.text.x=element_text(angle=45, hjust=1), 
		axis.ticks.x=element_blank(), 
		axis.ticks.y=element_blank()) +
	coord_fixed(ratio=0.8)

pdf(paste0(result_dir, "3_Cell_Clin_Box/3.3_Cell_stat_sig.pdf"))
print(p_sig)
dev.off()

topptx(p_sig, filename=paste0(result_dir, "3_Cell_Clin_Box/3.3_Cell_stat_sig.pptx"))
