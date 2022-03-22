# modified according to reviewer1's suggestion: use RUVg normalization instead

rm(list=ls())
################################## BC progression ##################################
# install.packages("BiocManager")
# BiocManager::install(update=TRUE, ask=FALSE)
# BiocManager::install(version='3.13')
# BiocManager::install("impute", force=TRUE)
# BiocManager::install("RUVSeq")
library(openxlsx)
library(preprocessCore)
library(impute)
library(RUVSeq)
library(RColorBrewer)


# set working directory
data_dir <- "C:/0_xmsun/xmsun/Graduate/20210224_NMIBC/Data/"
result_dir <- "C:/0_xmsun/xmsun/Graduate/20210224_NMIBC/Results/"
dir.create(paste0(result_dir, "7_Housekeeping_normalization/"))


########## load data
# expression count matrix including: exp_data_all, exp_data_all_house, exp_data_enough, exp_data_all_house_knn, exp_data_all_house_knn_quant, exp_data_all_house_knn_quant_log
load(file=paste0(result_dir, "7_Housekeeping_normalization_v1_before_review/7.99_ALL_norm_process.Rdata"))

# phenotype matrix
pheno <- read.xlsx(paste0(data_dir, "combined/Clinicopathologic_v2.xlsx"), na.strings=c("","NA","Nan","#N/A"))
rownames(pheno) <- pheno$Sample_name

# housekeeping genes
housekeeping_genes <- read.xlsx(paste0(result_dir, "7_Housekeeping_normalization_v1_before_review/7.0_Housekeeping_genes.xlsx"))
my_housekeeping <- housekeeping_genes[housekeeping_genes$inBladder_BladderCancer=="yes", "Gene_name"]
GEP_housekeeping <- housekeeping_genes[housekeeping_genes$Source=="GEP", "Gene_name"]


########## data processing
# expression count matrix: genes in row, samples in column, numeric format
exp_matrix <- apply(exp_data_all, 2, function(x) {round(x, digits=0)})
rownames(exp_matrix) <- rownames(exp_data_all)
exp_matrix <- na.omit(exp_matrix)

# phenotype matirx: samples in row, phenotypes in column
pheno_matrix <- pheno[match(colnames(exp_matrix), rownames(pheno)),]
pheno_matrix$Source_dataset <- as.factor(pheno_matrix$Source_dataset)
pheno_matrix$Platform <- as.factor(pheno_matrix$Platform)
pheno_matrix$Staging <- factor(pheno_matrix$Staging, levels=c("T0", "Ta", "Ta-1", "T1", "Tis"))

########## PCA & RLE plot: Non-normalized expression matrix
set <- newSeqExpressionSet(as.matrix(exp_matrix), phenoData = pheno_matrix)
colors <- brewer.pal(9, "Set1")
pdf(paste0(result_dir, "7_Housekeeping_normalization/0_Non_norm.pdf"))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[pheno_matrix$Source_dataset])
plotPCA(set, col=colors[pheno_matrix$Source_dataset])
dev.off()

########## PCA & RLE plot: RUVg-normalized expression matrix, with seven pre-defined housekeeping genes
my_housekeeping_kept <- my_housekeeping[my_housekeeping %in% rownames(exp_matrix)]
set1 <- RUVg(set, my_housekeeping_kept, k=10)
colors <- brewer.pal(9, "Set1")

pdf(paste0(result_dir, "7_Housekeeping_normalization/1_RUVg_house.pdf"))
plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[pheno_matrix$Source_dataset])
plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[pheno_matrix$Platform])
plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[pheno_matrix$Staging])
plotPCA(set1, col=colors[pheno_matrix$Source_dataset])
plotPCA(set1, col=colors[pheno_matrix$Platform])
plotPCA(set1, col=colors[pheno_matrix$Staging])
dev.off()

########## PCA & RLE plot: RUVg-normalized expression matrix, with empirical control genes
design <- model.matrix(~ Staging, data=pData(set))
y <- DGEList(counts=counts(set), group=pData(set)$Staging)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]
set2 <- RUVg(set, empirical, k=10)
colors <- brewer.pal(9, "Set1")

pdf(paste0(result_dir, "7_Housekeeping_normalization/2_RUVg_empirical.pdf"))
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[pheno_matrix$Source_dataset])
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[pheno_matrix$Platform])
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[pheno_matrix$Staging])
plotPCA(set2, col=colors[pheno_matrix$Source_dataset])
plotPCA(set2, col=colors[pheno_matrix$Platform])
plotPCA(set2, col=colors[pheno_matrix$Staging])
dev.off()


# output count matrix
write.xlsx(counts(set), file=paste0(result_dir, "7_Housekeeping_normalization/7.1_exp.xlsx"))
write.xlsx(normCounts(set1), file=paste0(result_dir, "7_Housekeeping_normalization/7.2_RUVg_house.xlsx"))
write.xlsx(normCounts(set2), file=paste0(result_dir, "7_Housekeeping_normalization/7.2_RUVg_empirical.xlsx"))
save(set, set1, set2, file=paste0(result_dir, "7_Housekeeping_normalization/7.99_RUVg_norm_process.Rdata"))
