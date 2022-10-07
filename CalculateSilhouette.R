rm(list=ls())
library(dplyr)
library(Seurat)
library(patchwork)

setwd('D:/test')
pbmc.data <- Read10X(data.dir = 'filtered_gene_bc_matrices/hg19/')
pbmc <- CreateSeuratObject(counts = pbmc.data,project = 'pbmc3k',min.cells = 3,min.features = 200)
pbmc

pbmc[['percent.mt']] <- PercentageFeatureSet(pbmc,pattern = '^MT-')
pbmc <- subset(pbmc,subset = nFeature_RNA > 200 & nFeature_RNA <2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc,selection.method = 'vst',nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc),10)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.genes)

pbmc <- RunPCA(pbmc,features = VariableFeatures(object = pbmc))

pbmc <- JackStraw(pbmc,num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc,dims = 1:20)

CalculateSilhouette<- function(object, dims = 1:50){
  if (length(dims) > ncol(object@reductions$pca@cell.embeddings)) {
    stop("please specify PCA dims smaller than calculated")
  }
  cell_distance<- dist(object@reductions$pca@cell.embeddings[, dims])
  # or as.integer
  cell_cluster<- as.numeric(as.character(Idents(object)))
  silhouette_score<- cluster::silhouette(cell_cluster, cell_distance)
  silhouette_score<- tibble::tibble(cluster = silhouette_score[,1],
                                    width = silhouette_score[,3],
                                    cell = colnames(object)) %>%
    dplyr::mutate(cluster = as.factor(cluster))
  return(silhouette_score)
}
pbmc <- FindNeighbors(pbmc,dims = 1:20)
pbmc <- FindClusters(pbmc,resolution = 0.3)
mean(CalculateSilhouette(pbmc,dims = 1:7)$width)
