# sample code
# analysing PBMC data
####################################################################

# load require libraries and codes

library(Seurat)
library(Matrix)
# library(SingleR)  # optional
# library(celldex)  # optional
library(igraph)
library(scales)
library(ggplot2)
library(gplots)
library(minet)

# BootCellNet codes
source("bootcellnet.R")

#=============================
# loading and processing the PBMC data

# data retrieved from 10x
# for example...
# wget https://cf.10xgenomics.com/samples/cell-exp/3.1.0/manual_5k_pbmc_NGSC3_ch5/manual_5k_pbmc_NGSC3_ch5_filtered_feature_bc_matrix.tar.gz
# gunzip manual_5k_pbmc_NGSC3_ch5_filtered_feature_bc_matrix.tar.gz

data <- Read10X("filtered_feature_bc_matrix/")
so <- CreateSeuratObject(data)

ngenes <- dim(so)[1] # number of genes in the data. We use this later.

# filtering
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)

# normalization etc
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
so <- ScaleData(so, features=rownames(so))

so <- RunPCA(so, features=VariableFeatures(object=so), dim=1:20)
# DimPlot(so, reduction = "pca") + NoLegend() # optional
so <- RunUMAP(so, dims=1:20)
# DimPlot(so, reduction = "umap") + NoLegend() # optional
so <- FindNeighbors(so, dims=1:10)

# optional: usual clustering workflow
#
# so <- FindNeighbors(so, dims=1:10)
# so <- FindClusters(so, resolution=.5, algorithm=1)
# 
# DimPlot(so, reduction = "umap") + NoLegend()
# 
# saveRDS(so, "pbmc_test.rds")

####################################################################
####################################################################
# BootCellNet procedure starts from here

#====================================
# k-NN smoothing
# here, we used k=10
data <- so@assays$RNA$data
umap.coord <- so@reductions$umap@cell.embeddings

data.kns <- knnsmo(data=data, umap=umap.coord, 
                   K=10, file="pbmc_kns.mtx.gz")

#---------------------
# optional 
# if the smoothing takes time...
# knnsmopbmc_kns.mtx.gz.dir(data=data, umap=umap.coord, K=10, file="pbmc_kns.mtx.gz")
# Please cat the header + the matrix market file by running
# cat header_pbmc_kns.mtx.gz pbmc_kns.mtx.gz >wh_pbmc_kns.mtx.gz
# 
# # if you stopped the workflow, start over from here
# data.kns <- as(readMM("wh_pbmc_kns.mtx.gz"), "dgCMatrix")

rownames(data.kns) <- rownames(data)
colnames(data.kns) <- colnames(data)

#====================================
# resampling

#-------------------
# preparation
# read in the annotation table for PBMC data, which contains GeneIDs
annot <- read.table("pbmc_annot_geneid.txt", header=T, fill=T, na.strings=NA, stringsAsFactors = F)

# read in the table of transcription factros, 
# used in Joung et al, "A transcription factor atlas of directed differentiation", Cell, 2023
goi <- read.table("MORF_genelist_geneids.txt", header=T, stringsAsFactors = F)
annot.goi <- annot[annot[,2]%in%unlist(goi),1]

#-------------------
# resampling and aggregation

N = 3000  # number of cells to be resampled
M = 100   # number of resampling
n.vargenes = 2000 # number of variable expression genes to be used in calculation

bootcellnet(data.kns, N=N, M=M, nfeatures=n.vargenes, goi=annot.goi, directory = "pbmc_bcn/")

# directory "pbmc_bcn" is generated, and 100 files should be there.

# START OVER FROM HERE (if necessary)
agr.cells <- aggregateBCN(directory="pbmc_bcn/", R=0)

#====================================
# calc p-values
# resampling and shuffling

Mhat = 1000  # number of resampling for support calculation
K = 100      # number of shuffling

bootcellnet(data.kns, N=N, M=Mhat, suffix=NULL, goi=annot.goi, directory = "pbmc_bcn_shuffle/")

agr.cells.shuffle <- aggregateBCN(directory="pbmc_bcn_shuffle/", R=K)

node.pvals <- computeNodePvals(agr.cells$node.count, M*n.vargenes/ngenes)
edge.pvals <- computeEdgePvals(agr.cells$edge.count, M*agr.cells.shuffle$edge.count/(Mhat*K))

#====================================
# clean up graphs based on the p-values
# here, p.adjust is performed to obtain q-values
# we will incorporate this as a part of cleanGRN function in future
ar.graph <- cleanGRN(agr.cells$aracne, p.adjust(node.pvals, "BY"), p.adjust(edge.pvals, "BY"))

# check the graphs if necessary
# we will provide this as a function in future
pdf("pbmc-GRNs.pdf")

# aggregated graphs: contains all the nodes and edges
plot(ar.graph$graph.all$graph, layout=ar.graph$graph.all$layout,
     vertex.color="white", vertex.cex=0.5, vertex.size=10, vertex.label.cex=0.5,
     edge.width=ar.graph$graph.all$weight)

# clean graph
plot(ar.graph$graph.clean$graph, layout=ar.graph$graph.clean$layout,
     vertex.color="white", vertex.cex=0.5, vertex.size=10, vertex.label.cex=0.5,
     edge.width=ar.graph$graph.clean$weight)
dev.off()


#====================================
# compute MDS
ar.graph <- MDSGRN(ar.graph)

pdf("pbmc-GRN_MDS_graphs.pdf")
node.color <- rep("white", vcount(ar.graph$graph.clean$graph))
node.color[names(V(ar.graph$graph.clean$graph)) %in% ar.graph$graph.clean$MDS] <- "red"

plot(ar.graph$graph.clean$graph, layout=ar.graph$graph.clean$layout,
     vertex.color=node.color, vertex.cex=0.5, vertex.size=10, vertex.label.cex=0.5,
     vertex.label.color="black",
     edge.width=ar.graph$graph.clean$weight, 
     vertex.shape=ar.graph$graph.clean$node_shape_MDS, 
     main="MDS")
dev.off()


#====================================
# now ready for clustering...

#-------------------
# first, determine the resolution to obtain cell clumps
# what resolution is suitable?
clust.tab <- list()
for(r in c(1,2,3,4,6,8,10,15,20,50)){
  so2 <- FindClusters(so, graph.name="RNA_snn", resolution = r)
  clust.tab[[paste0("res.",r)]] <- table(Idents(so2))
}
rm(so2)

# check the contents of clust.tab
# and you will find that it starts to have cluster with less than 10 cells after resolution > 10
# thus use resolution=10
so <- FindClusters(so, graph.name="RNA_snn", resolution = 10)

#-------------------
# clustering by MDS
Idents(so) <- "RNA_snn_res.10"

data.mds <- avgGS(ar.graph$graph.clean$MDS, so)
data.mds.z <- t(apply(data.mds2, 1, function(x) {(return(x-mean(x))/sd(x))}))

data.mds.genecl <- hclust(as.dist(t(1-cor(t(data.mds.z)))), method="ward.D2")
data.mds.cellcl <- hclust(as.dist(1-cor(data.mds.z)), method="ward.D2")

pdf("pbmc-MDS_hmp_noclust.pdf")
heatmap.2(data.mds.z, col=colorRampPalette(c("blue","white","red"))(50), 
          scale="row", Rowv=as.dendrogram(data.mds.genecl), Colv=as.dendrogram(data.mds.cellcl), 
          cexRow=.5, cexCol=.5, trace="none")
dev.off()

# check the heatmap by eyes...
# determine cluster number
# here we set it 13
n.clust <- 13

pdf("pbmc-MDS_hmp.pdf")
heatmap.2(data.mds.z, col=colorRampPalette(c("blue","white","red"))(50), 
          scale="row", Rowv=as.dendrogram(data.mds.genecl), Colv=as.dendrogram(data.mds.cellcl), 
          cexRow=.5, cexCol=.5, trace="none",
          ColSideColors = hue_pal()(n.clust)[cutree(data.mds.cellcl,k=n.clust)])
dev.off()

# map the clusters to Seurat object
Idents(so) <- "RNA_snn_res.10"
so <- newclst(so, cutree(data.mds.cellcl,k=n.clust))


pdf("pbmc_MDS_clst.pdf")
Idents(so) <- "RNA_snn_res.10"
DimPlot(so, reduction="umap", label = T)+NoLegend()
Idents(so) <- "MDS_clst"
DimPlot(so, reduction="umap", label = T)
DimPlot(so, reduction="umap", label = T)+NoLegend()
DimPlot(so, reduction="umap", label = F)+NoLegend()

for(i in 1:n.clust){
  # DimPlot(so, reduction="umap", label = T, cells = colnames(so)[so@active.ident==i])
  dimcol <- rep("#EEEEEE10", n.clust); dimcol[i] <- hue_pal()(n.clust)[i]
  p <- DimPlot(so, reduction="umap", label = T, cols = dimcol)
  print(p)
}

dev.off()


#====================================
# now you can use the clustering result to any analysis
# for example

Idents(so) <- "MDS_clst"
mds.clust.markers <- FindAllMarkers(so, only.pos=TRUE, min.pct=.1, logfc.threshold=1)
write.table(mds.clust.markers, file="pbmc_MDS_clust_markers.txt", quote=F, sep="\t", row.names=F)


