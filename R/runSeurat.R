library(Seurat)
source("R/tximport/tximport.R")
args <- commandArgs(trailingOnly = TRUE)
files <- file.path(args[1])
iter = args[2]
setwd(args[3])
cluster_resolution = args[4]
file.exists(files)
txi <- tximport(files, type="alevin", em_iter = iter)
sample <- CreateSeuratObject(counts = txi$counts , min.cells = 3, min.features = 0, project = "10X_sample")

print("Sample before filter:")
sample
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
sample <- subset(sample, subset = nFeature_RNA > 200 & percent.mt < 20)
print("Sample after filter:")
sample

sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)

sample <- ScaleData(sample)
sample <- RunPCA(sample, features = VariableFeatures(object = sample), verbose = FALSE)
sample <- FindNeighbors(sample, dims = 1:10)
sample <- FindClusters(sample, resolution = cluster_resolution)
sample <- RunUMAP(sample, dims = 1:10)

cluster_filename <- paste(iter, 'iters_cluster.txt', sep = '')
write.table(Idents(sample), file = cluster_filename)

UMAP_filename <- paste(iter, 'iters_UMAP.pdf', sep = '')
pdf(UMAP_filename)
DimPlot(sample, reduction = "umap")
dev.off()

pca.mat <- Loadings(object = sample, reduction = "pca")
cluster_num <- length(unique(Idents(sample)))
sig_all <- matrix(, nrow = 50, ncol = 0)
for(i in c(0:(cluster_num - 1))) {
	sub_sample <- subset(sample, ident = i)
	scale.data <- GetAssayData(object = sub_sample, slot = "scale.data")
	sig_exp <- rowMeans(scale.data)
	gene_name <- rownames(scale.data)
	sig_pca <- matrix(, nrow = 50, ncol = 1)
	colnames(sig_pca) <- c(i)
	for (j in c(1:50)) {
		tmp_sum <- 0
		pc_label <- paste("PC", j, sep = '_')
		for (w in c(1:length(gene_name))) {
			tmp_sum <- tmp_sum + (sig_exp[gene_name[w]] * pca.mat[gene_name[w], pc_label])
		}
		sig_pca[j,1] <- tmp_sum
	}

	if (i == 0) {
		sig_all <- sig_pca
	} else {
		sig_all <- cbind(sig_all, sig_pca)
	}
}

sig_filename <- paste(iter, 'Sig_PCA.txt', sep = '')
write.table(sig_all, file = sig_filename)

# Used fixed selected dimension
selected_dim <- 15

PCA <- Embeddings(object = sample[["pca"]])[, 1:selected_dim]
PCA_filename <- paste(iter, 'PCA_coordinate.txt', sep = '')
write.table(PCA, file = PCA_filename, sep = ' ', row.names = T, col.names = F, quote = F)




