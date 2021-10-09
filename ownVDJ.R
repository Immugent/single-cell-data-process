
rm(list =ls())
library(Seurat)


setwd("")
sample_old <- Read10X(data.dir = "")
#colnames(pbmc.data) <- paste('BC11_TUMOR1', colnames(pbmc.data), sep = '_')

sample_old <- CreateSeuratObject(counts = sample_old, project = "Old", min.cells = 3, min.features = 200)
sample_old


sample_old[["percent.mt"]] <- PercentageFeatureSet(sample_old, pattern = "^mt-")
VlnPlot(sample_old, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),col=3, ncol = 3)

sample_old <- subset(sample_old, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & percent.mt < 10)
#sample_old
#table(sample_old$orig.ident)

sample_old <- NormalizeData(sample_old, normalization.method = "LogNormalize", scale.factor = 10000)
sample_old
sample_old<- FindVariableFeatures(sample_old, selection.method = "vst", nfeatures = 3000)

all.genes <- rownames(sample_old)
sample_old <- ScaleData(sample_old, features = all.genes)
sample_old <- ScaleData(sample_old, vars.to.regress = "percent.mt")#

sample_old <- RunPCA(sample_old, features = VariableFeatures(object = sample_old))
ElbowPlot(sample_old, ndims = 50)

sample_old <- FindNeighbors(sample_old, dims = 1:50)

sample_old <- FindClusters(sample_old, resolution = 1)
sample_old <- RunUMAP(sample_old, dims = 1:50)
sample_old <- RunTSNE(sample_old, dims = 1:50)

sample_vdj<-sample_old
tcr <- read.csv("")
tcr[1:4,]
tcr <- tcr[!duplicated(tcr$barcode), ] 
tcr[1:4,]
names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id" 
tcr[1:4,]
table(tcr$v_gene)
clono <- read.csv("")   
head(clono)
tcr <- merge(tcr, clono)
tcr[1:4,]
rownames(tcr) <- tcr[,2]
clono_seurat <- AddMetaData(object=sample_vdj, metadata=tcr)   
table(rownames(clono_seurat@meta.data)  %in% rownames(tcr))
#clono_seurat<- subset(clono_seurat,cells =  rownames(tcr))

setwd("")
sample_1 <- readRDS("")

head(clono_seurat@meta.data)
sample_1<-clono_seurat
table(clono_seurat$chain)


DimPlot(sample_1, reduction="umap", group.by="chain",label=FALSE,pt.size=1.5)
DimPlot(sample_1, reduction="tsne", group.by="chain",label=FALSE,pt.size=1.5)
FeaturePlot(sample_1,features="frequency", reduction="tsne",pt.size=1.5)
FeaturePlot(sample_1, features = c("Cd3e", "Cd3d", "Cd8a","Cd4"), reduction = "tsne",cols = c("gray", "red"),pt.size=1)

setwd("")
saveRDS(sample_1, file = "")
save(sample_1,file="")

                                                  
rm(list = ls())
setwd("")
sample_young <- Read10X(data.dir = "")
#colnames(pbmc.data) <- paste('BC11_TUMOR1', colnames(pbmc.data), sep = '_')

sample_young <- CreateSeuratObject(counts = sample_young, project = "young", min.cells = 3, min.features = 200)
sample_young

sample_young[["percent.mt"]] <- PercentageFeatureSet(sample_young, pattern = "^mt-")

VlnPlot(sample_young, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),col=3, ncol = 3)

sample_young <- subset(sample_young, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & percent.mt < 10)
sample_young

sample_young <- NormalizeData(sample_young, normalization.method = "LogNormalize", scale.factor = 10000)
sample_young

sample_young<- FindVariableFeatures(sample_young, selection.method = "vst", nfeatures = 3000)

all.genes <- rownames(sample_young)
sample_young <- ScaleData(sample_young, features = all.genes)
sample_young <- ScaleData(sample_young, vars.to.regress = "percent.mt")

sample_young <- RunPCA(sample_young, features = VariableFeatures(object = sample_young))
ElbowPlot(sample_young, ndims = 50)
sample_young <- FindNeighbors(sample_young, dims = 1:50)

sample_young <- FindClusters(sample_young, resolution = 1)
sample_young <- RunUMAP(sample_young, dims = 1:50)
sample_young <- RunTSNE(sample_young, dims = 1:50)

sample_vdj<-sample_young
tcr <- read.csv("")
tcr[1:4,]
tcr <- tcr[!duplicated(tcr$barcode), ]
names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
table(tcr$v_gene)
clono <- read.csv("")
head(clono)
tcr <- merge(tcr, clono)
rownames(tcr) <- tcr[,2]
clono_seurat <- AddMetaData(object=sample_vdj, metadata=tcr)   
table(rownames(clono_seurat@meta.data)  %in% rownames(tcr))

setwd("")
sample_2 <- readRDS("")
head(clono_seurat@meta.data)
sample_2<-clono_seurat
table(clono_seurat$chain)
table(sample_2$clonotype_id)
DimPlot(sample_2, reduction="umap", group.by="chain",label=FALSE,pt.size=1.5)
DimPlot(sample_2, reduction="tsne", group.by="clonotype_id",label=FALSE,pt.size=1.5)
DimPlot(sample_2, reduction="tsne", group.by="chain",label=FALSE,pt.size=1.5)
FeaturePlot(sample_2,features="frequency", reduction="tsne",pt.size=1.5)
FeaturePlot(sample_2, features = c("Cd3e", "Cd3d", "Cd8a","Cd4"), reduction = "tsne",cols = c("gray", "red"),pt.size=1)
RidgePlot(sample_2,features=rownames(mk),group.by = "chain")

setwd("")
saveRDS(sample_2, file = "")
save(sample_2,file="")

                            
library(scRepertoire)
rm(list=ls())
setwd("")
csv1 <- read.csv("filtered_contig_annotations.csv", stringsAsFactors = F)
csv2 <- read.csv("filtered_contig_annotations.csv", stringsAsFactors = F)
contig_list1 <- list(csv1,csv2)


combined <- combineTCR(contig_list1, samples = 
                         c("TUY", "TUO"), ID = c("T", "T"), cells ="T-AB")
head(combined[[1]])


p1 <- quantContig(combined, cloneCall="gene+nt", scale = F)
p2 <- quantContig(combined, cloneCall="gene+nt", scale = T)
plotc = p1 + p2
plotc

quantContig(combined, cloneCall="gene+nt", scale = T, exportTable = T)

p1 <- lengthContig(combined, cloneCall="aa", chains = "combined") 
p2 <- lengthContig(combined, cloneCall="aa", chains = "single")
plotc = p1 + p2
plotc
#ggsave('VDJ/lengthContig.png', plotc, width = 8, height = 4)


p1 = compareClonotypes(combined, numbers = 10, samples = c("TUY_T", "TUO_T"), 
                       cloneCall="aa", graph = "alluvial")
p1
#ggsave('VDJ/compareClonotypes.png', p1, width = 8, height = 4)

clonalHomeostasis(combined, cloneCall = "aa")
clonalProportion(combined, cloneCall = "aa") 

clonalDiversity(combined, cloneCall = "aa", group = "samples")


seurat <- readRDS("merge\\data\\merged_tutorial.rds")
names(seurat@meta.data$sample_type)
table(seurat@meta.data$sample_type)
dim(seurat)
seurat@meta.data$sample_type

DimPlot(seurat, label = T) #+ NoLegend()

names(seurat@meta.data)

combined$TUY_T$sample
seurat <- combineExpression(combined, seurat, cloneCall="aa", groupBy = "sample")
names(seurat@meta.data)

seurat$cloneType <- factor(seurat$cloneType, levels = c("Hyperexpanded (100 < X <= 500)", 
                                                        "Large (20 < X <= 100)", "Medium (5 < X <= 20)", 
                                                        "Small (1 < X <= 5)", "Single (0 < X <= 1)", NA))
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
DimPlot(seurat, group.by = "cloneType") + 
  scale_color_manual(values = colorblind_vector(5), na.value="grey")


seurat <- highlightClonotypes(seurat, cloneCall= "aa", 
                              sequence = c("CAVNGGSQGNLIF_CSAEREDTDTQYF", "NA_CATSATLRVVAEKLFF"))
DimPlot(seurat, group.by = "highlight")

alluvialClonotypes(seurat, cloneCall = "aa", 
                   y.axes = c("Patient", "cluster", "Type"), 
                   color = "CAVNGGSQGNLIF_CSAEREDTDTQYF") + 
  scale_fill_manual(values = c("grey", colorblind_vector(1)))

alluvialClonotypes(seurat, cloneCall = "aa", 
                   y.axes = c("Patient", "cluster", "Type"), 
                   color = "cluster")


combined2 <- expression2List(seurat, group = "cluster")
p1 = clonalDiversity(combined2, cloneCall = "aa") + ggtitle("clonalDiversity")
p2 = clonalHomeostasis(combined2, cloneCall = "aa") + ggtitle("clonalHomeostasis")
p3 = clonalProportion(combined2, cloneCall = "aa") + ggtitle("clonalProportion")
p4 = clonalOverlap(combined2, cloneCall="aa", method="overlap")  + ggtitle("clonalOverlap")  
plotc = (p2|p3)/(p1|p4)
ggsave('VDJ/clonal.png', plotc, width = 10, height = 8)





