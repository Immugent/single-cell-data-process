#20200930  

setwd("")
##载入seurat包
library(dplyr)
library(Seurat)
rm(list=ls())
##读入pbmc数据
pbmc.data <- Read10X(data.dir = "")
##创建Seurat对象与数据过滤
data <- CreateSeuratObject(counts = pbmc.data, project = "Young", min.cells = 3, min.features = 200)
data
##读入pbmc1数据
pbmc1.data <- Read10X(data.dir = "")
##创建Seurat对象与数据过滤
data1 <- CreateSeuratObject(counts = pbmc1.data, project = "Old", min.cells = 3, min.features = 200)
data1
merged<-merge(data,data1)
table(merged$orig.ident)

##计算每个细胞的线粒体基因转录本数的百分比（%）,使用[[ ]] 操作符存放到metadata中
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^mt-")

##展示基因及线粒体百分比
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(merged, features = "percent.mt",y.max=10)

##展示基因及线粒体百分比
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(merged, features = "percent.mt",y.max=10)

plot1 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

str(merged)
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & percent.mt < 10)#
#apply(merged@assays$RNA@counts,2,mean)
median(merged$nFeature_RNA)
median(pbmc1$nFeature_RNA)

table(merged$orig.ident)
##表达量数据标准化,LogNormalize的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000 )
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
#merged <- NormalizeData(merged) 或者用默认的
merged
##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 3000)

##提取表达量变变化最高的10个基因；
top10 <- head(VariableFeatures(merged), 10)
top10

plot1 <- VariableFeaturePlot(merged)
plot2 <- LabelPoints(plot = plot1, points = top10)
CombinePlots(plots = list(plot1, plot2))

##PCA分析：
#PCA分析数据准备,使用ScaleData()进行数据归一化；默认只是标准化高变基因（2000个）,速度更快,不影响PCA和分群,但影响热图的绘制。
#merged <- ScaleData(merged,vars.to.regress ="percent.mt")

##而对所有基因进行标准化的方法如下：
all.genes <- rownames(merged)
merged <- ScaleData(merged, features = all.genes)
#用高变基因标准化代码如下，默认参数
merged <- ScaleData(merged, vars.to.regress = "percent.mt")#耗时

##线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；
merged <- RunPCA(merged, features = VariableFeatures(object = merged))

##检查PCA分群结果, 这里只展示前5个PC,每个PC只显示5个基因；
print(merged[["pca"]], dims = 1:5, nfeatures = 5)

##展示主成分基因分值
VizDimLoadings(merged, dims = 1:2, reduction = "pca")

##绘制pca散点图
DimPlot(merged, reduction = "pca",pt.size = 1.5)

##画第1个或15个主成分的热图；
DimHeatmap(merged, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(merged, dims = 1:15, cells = 500, balanced = TRUE)

#方法2：肘部图（碎石图）,基于每个主成分对方差解释率的排名。
ElbowPlot(merged, ndims = 50)
#dev.off()
##分群个数这里选择30,建议尝试选择多个主成分个数做下游分析,对整体影响不大；在选择此参数时,建议选择偏高的数字（为了获取更多的稀有分群,“宁滥勿缺”）；有些亚群很罕见,如果没有先验知识,很难将这种大小的数据集与背景噪声区分开来。

##非线性降维（UMAP/tSNE)基于PCA空间中的欧氏距离计算nearest neighbor graph,优化任意两个细胞间的距离权重（输入上一步得到的PC维数）。
merged <- FindNeighbors(merged, dims = 1:50)

##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
merged <- FindClusters(merged, resolution = 1.2)
##使用Idents（）函数可查看不同细胞的分群；
head(Idents(merged), 5)
table(merged$seurat_clusters)
##Seurat提供了几种非线性降维的方法进行数据可视化（在低维空间把相似的细胞聚在一起）,比如UMAP和t-SNE,运行UMAP需要先安装'umap-learn'包,这里不做介绍。
merged <- RunTSNE(merged, dims = 1:50)#耗时
merged1 <- RunUMAP(merged, dims = 1:50)
##用DimPlot()函数绘制散点图,reduction = "tsne",指定绘制类型；如果不指定,默认先从搜索 umap,然后 tsne, 再然后 pca；也可以直接使用这3个函数PCAPlot()、TSNEPlot()、UMAPPlot()； cols,pt.size分别调整分组颜色和点的大小；

DimPlot(merged,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(merged1,reduction = "umap",label = TRUE,pt.size = 1.5)

plot1<-DimPlot(subset(merged, subset = orig.ident=='Old'),reduction = "tsne",label = TRUE,pt.size = 1.5)
plot2<-DimPlot(subset(merged, subset = orig.ident=='Young'),reduction = "tsne",label = TRUE,pt.size = 1.5)
CombinePlots(plots = list(plot1, plot2))

DimPlot(merged,reduction = "tsne",label = TRUE,group.by="RNA_snn_res.1.2",label.size = 5,pt.size = 1.5)

##细胞周期归类
merged<- CellCycleScoring(object = merged, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
head(x = merged@meta.data)
DimPlot(merged,reduction = "tsne",label = TRUE,group.by="Phase",label.size = 0,pt.size = 1.5)
merged$RNA_snn_res.1.2

#对分群进行重新命名
new.cluster.ids <- c("Macrophages", "Macrophages", "CD3+CD4-CD8-T", "cDC", "Monocytes", "CD8+CD4-T", "Granulocytes", "Macrophages", "Granulocytes",
                     "NK", "CD4+CD8-T", "cDC", "CD4+CD8-T", "CD3+CD4-CD8-T","B cells", "pDC", "cDC", "CD8+CD4-T", "CD8+CD4-T",
                     "Macrophages", "Macrophages", "Macrophages","Monocytes", "B cells", "Basophils")
names(new.cluster.ids) <- levels(merged)
merged <- RenameIdents(merged, new.cluster.ids)
DimPlot(merged, reduction = "tsne", label = TRUE, pt.size = 1.5) + NoLegend()

##比较2组数据中特定cluster的差异基因
merged@meta.data$sample_type <- paste(merged@meta.data$orig.ident, merged@meta.data$RNA_snn_res.1.2, sep = "_")#增加了一列，可以更多方向cluster的比较
merged@meta.data$celltype <- paste(merged@meta.data$orig.ident, Idents(merged), sep = "_")#增加了一列，可以更多方向cluster的比较
merged@meta.data$origin <- paste(merged@meta.data$RNA_snn_res.1.2, Idents(merged), sep = "_")#增加了一列，可以更多方向cluster的比较
table(merged$celltype)

a<- table(merged@meta.data$origin)
write.csv(a,"origin.csv")
b <- table(merged@meta.data$sample_type)
write.csv(b,"sample_type.csv")
c <- table(merged@meta.data$celltype)
write.csv(c,"celltype.csv")
#验证
# B cell NK
FeaturePlot(merged, features = c("Ptprc","Cd3e","Cd4", "Cd8a", "Cd79a", "Cd79b",  "Nkg7", "Itgam", "Csf1r","Itgax","Flt3","Csf3r"), reduction = "tsne",cols = c("gray", "red"))

FeaturePlot(merged, features = c("Cd19", "Cd38","Ms4a1", "Cd37", "Pax5", "Ncam1", "Cd7", "Klf2","Nkg7","Fosb"),reduction = "umap",cols = c("gray", "red"))
# DC Mφ Mo
FeaturePlot(merged, features = c("Itgax", "Cd80", "Nrp1","Cd83", "Cd86", "Cd14", "F13a1", "Itgam", "Fcer1a", "Siglech","Ly6c1","Ly6g"),reduction = "umap",cols = c("gray", "red"))
# gdT  NKT  Neu
FeaturePlot(merged, features = c("Trdc", "Sox13", "Myb", "Klrc1", "Jun", "Prdm16", "Klrc3", "Klrc2", "Zbtb16","Tbx21","Nkg7","Klrb1c"),reduction = "umap",cols = c("gray", "red"))



VlnPlot(merged, features = c("Cd3e", "Cd3d","Cd4", "Cd8a","Itgam", "Itgax", "Cd79b","Csf3r","Csf1r"),pt.size = 0)




a<- table(merged@meta.data$origin)
write.csv(a,"origin.csv")
b <- table(merged@meta.data$sample_type)
write.csv(b,"sample_type.csv")
c <- table(merged@meta.data$celltype)
write.csv(c,"celltype.csv")
##存储结果
setwd("")
saveRDS(merged, file = "")
save(merged,file="") 
           
                          ########亚群分析############
library(dplyr)
library(Seurat)
rm(list=ls())
setwd("")
merged <- readRDS("merged_tutorial.rds")
merged.markers <- read.table("merged_allmarker.txt")

merged <- readRDS("zhong.rds")
table(Idents(merged))
table(merged$celltype)
merged
merged<-subset(merged, idents = c("Monocytes", "Macrophages", "CD3+CD4-CD8-T", "cDC",  "CD8+CD4-T","Granulocytes",
                                  "NK",  "CD4+CD8-T", "B cells", "pDC"))#c(0,1)
merged
new.cluster.ids <- c("Monocytes", "Macrophages", "CD3_T", "DC_1",  "CD8_T", "Granulocytes","Granulocytes",
                     "NK_cells", "CD4_T", "CD4_T", "B cells", "DC_2", "DC_2")
names(new.cluster.ids) <- levels(merged)
merged <- RenameIdents(merged, new.cluster.ids)
DimPlot(merged, reduction = "tsne", label = TRUE, pt.size = 1.5,label.size = 0) + NoLegend()
CombinePlots(plots = list(plot1, plot2))

plot1<-DimPlot(subset(merged, subset = orig.ident=='Old'),reduction = "tsne",label = TRUE,pt.size = 2,label.size = 0)+ NoLegend()
plot2<-DimPlot(subset(merged, subset = orig.ident=='Young'),reduction = "tsne",label = TRUE,pt.size = 2,label.size = 0)+ NoLegend()
CombinePlots(plots = list(plot1, plot2))


##寻找cluster1的marker
cluster1.markers <- FindMarkers(merged, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)#耗时

##寻找每一cluster的marker
merged.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
merged.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)#耗时

##存储marker
write.csv(merged.markers,file="marker.csv")
##绘制分cluster的热图，最好用全部基因画
top20 <- merged.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(merged, features = top20$gene,label = 0)

top10 <- merged.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(merged, features = top10$gene) + NoLegend()

#寻找具体的某两个群的差异基因
VlnPlot(merged, features = c("Cd4", "Cd8a"),split.by = "orig.ident",pt.size = 0)
VlnPlot(merged, features = c("Cd4", "Cd8a"),pt.size = 0)

table(merged$celltype)
marker<-FindMarkers(merged, group.by="celltype",ident.1 = "Old_Granulocytes", ident.2 = "Young_Granulocytes", logfc.threshold = 0, min.pct = 0.1)
head(marker)
View(marker)

write.csv(marker,"GN_oldyoung.csv")
#20201228 NK cell
merged<-subset(merged, idents = c("NK_cells"))
merged
table(merged$celltype)
VlnPlot(merged, features = c("Stat3", "Gzmb","Tnfrsf9", "Nfkbia","Klrc1", "Nkg7", "Nfkb1","Klrk1", "Prf1","Klra7","Klra4","Klrd1"),pt.size = 0,group.by = "celltype",ncol = 4)
DotPlot(merged, features = c( "Wfdc21", "Cd24a","Fcer1g","G0s2",  "Cd14","Clec4d", "Ccl4", "Ccl3","Il1rn","Bri3"), cols = c("blue", "red"),group.by = "celltype")
#20201228 B cell
merged<-subset(merged, idents = c("B cells"))
merged
table(merged$celltype)
VlnPlot(merged, features = c("Ighm", "Igkc","Cd44", "Tgfb1","Ly6a", "Jund", "Il2rg","Cd19", "Cd79a","Cd79b","Cd37","Cd74"),pt.size = 0,group.by = "celltype",ncol = 4)
DotPlot(merged, features = c( "Wfdc21", "Cd24a","Fcer1g","G0s2",  "Cd14","Clec4d", "Ccl4", "Ccl3","Il1rn","Bri3"), cols = c("blue", "red"),group.by = "celltype")



DotPlot(merged, features = c("cd74", "Itgax", "Cd3e", "Cd4", "Cd8a", "Ctla4","Pdcd1", "Foxp3", "Nkg7", "Itgam", "Cd3d"), cols = c("blue", "red"))
#CD8T STAT
DotPlot(merged, features = c("Tcf7", "Lef1", "Sell", "Il7r", "Prf1", "Fasl", "Gzmb", "Itgam", "Ctla4","Pdcd1","Havcr2","Tigit","Lag3"), cols = c("blue", "red"))
#B 细胞
DotPlot(merged, features = c("Prf1", "Ptprc", "Spib", "Sdc1", "Cd74", "Cd19", "Cd79a", "Cd79b"), cols = c("blue", "red"))
#淋系
DotPlot(merged, features = c("Cd3e", "Cd3d","Cd8a","Cd4", "Nkg7","Klra4", "Prf1","Cd19", "Cd79a", "Cd79b"), cols = c("blue", "red"))
#髓系
DotPlot(merged, features = c("Itgam", "Itgax","Cd14", "Fcgr3","Csf1r","Csf3r", "Ly6c2", "Mafb","Flt3","Kit"), cols = c("blue", "red"))
DotPlot(merged, features = c("Itgam", "Itgax","Cd14", "Fcgr3","Csf1r","Csf3r", "Ly6c2", "Mafb","Flt3","Kit","Cd3d","Cd8a","Cd4", "Nkg7","Klra4", "Prf1","Cd19", "Cd79a", "Cd79b"), cols = c("blue", "red"))


table(merged@meta.data$seurat_clusters)
                        
                                 ##########  T 细胞分群########
setwd("")
merged <- readRDS("sub_Tcell.rds")
sub_Tcell.markers <- read.table("sub_Tcell_allmarker.txt")

sub_Tcell<-subset(merged, idents = c('2',"5","10","12","13","17","18"))#c(0,1)
table(sub_Tcell@meta.data$orig.ident)
table(sub_Tcell@meta.data$RNA_snn_res.1.2)

new.cluster.ids <- c("CD8+T", "Macrophages_1", "Monocytes_1", "DC_1", "Macrophages_2", "Naive CD8+T", "Neutrophils", "Monocytes_2", "Granulocytes",
                     "NK", "CD4+ T", "DC_2", "Treg", "Actived B", "DC_3", "DC_4", "Exhausted CD8+T", "cytotoxic CD8+T",
                     "Macrophages_3", "Monocytes_3", "Macrophages_4", "B cells", "Basophils")
names(new.cluster.ids) <- levels(merged)
merged <- RenameIdents(merged, new.cluster.ids)
DimPlot(merged, reduction = "tsne", label = TRUE, pt.size = 1.5) + NoLegend()
write.table(Idents(merged),"cell-type.xls",sep="\t",quote = F)

##存储细胞归类结果
saveRDS(merged, file = "merged_final.rds")
setwd("")
merged <- readRDS("merged_final.rds")


##########对取出的cluster重新进行分群
sub_Tcell <- FindVariableFeatures(sub_Tcell, selection.method = "vst", nfeatures = 3000)

all.genes <- rownames(sub_Tcell)
sub_Tcell <- ScaleData(sub_Tcell, features = all.genes)
sub_Tcell <- ScaleData(sub_Tcell, vars.to.regress = "percent.mt")

sub_Tcell <- RunPCA(sub_Tcell, features = VariableFeatures(object = sub_Tcell))
ElbowPlot(sub_Tcell, ndims = 50)
sub_Tcell <- FindNeighbors(sub_Tcell, dims = 1:30)

sub_Tcell <- FindClusters(sub_Tcell, resolution = 0.52)
table(sub_Tcell$seurat_clusters)
table(sub_Tcell$RNA_snn_res.0.5)
sub_Tcell <- RunTSNE(sub_Tcell, dims = 1:30)
sub_Tcell <- RunUMAP(sub_Tcell, dims = 1:30)
DimPlot(sub_Tcell, reduction = "tsne",label = TRUE, pt.size = 1.5)
DimPlot(sub_Tcell,reduction = "umap",label = TRUE,pt.size = 1.5)



VlnPlot(sub_Tcell, features = c("Cd4","Foxp3", "Cd8a"),split.by = "seurat_clusters",pt.size = 0)
VlnPlot(sub_Tcell, features = c("Cd4", "Foxp3", "Cd8a"),split.by = "orig.ident",pt.size = 0)
VlnPlot(sub_Tcell, features = c("Cd3e", "Cd3d"),split.by = "RNA_snn_res.0.5",pt.size = 0)
VlnPlot(sub_Tcell, features = c("Ctla4","Pdcd1"),split.by = "RNA_snn_res.0.5",pt.size = 0)
VlnPlot(sub_Tcell, features = c("Cxcr6","Tcf7"),split.by = "RNA_snn_res.0.5",pt.size = 0)
VlnPlot(sub_Tcell, features = c("Prf1","Gzmb"),split.by = "RNA_snn_res.0.5",pt.size = 0)

##比较2组数据中特定cluster的差异基因
sub_Tcell@meta.data$sample_type <- paste(sub_Tcell@meta.data$orig.ident, sub_Tcell@meta.data$RNA_snn_res.0.52, sep = "_")#增加了一列，可以更多方向cluster的比较
sub_Tcell@meta.data$seurat_clusters
table(sub_Tcell@meta.data$sample_type)

marker<-FindMarkers(merged, group.by="sample_type",ident.1 = "Old_1", ident.2 = "Young_1", logfc.threshold = 0, min.pct = 0.1)
head(marker)
View(marker)
write.table(marker,"sample_typemarker.txt",sep="\t")

##寻找每一cluster的marker
sub_Tcell.markers <- FindAllMarkers(sub_Tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sub_Tcell.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)#耗时

##绘制分cluster的热图
top20 <- sub_Tcell.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(sub_Tcell, features = top20$gene) + NoLegend()

top10 <- sub_Tcell.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(sub_Tcell, features = top10$gene) + NoLegend()


##存储marker
write.table(sub_Tcell.markers,file="sub_Tcell_allmarker.txt")
write.csv(sub_Tcell.markers,file="sub_Tcell_allmarker.csv")

setwd("")
saveRDS(sub_Tcell, file = "sub_Tcell.rds")
save(sub_Tcell,file="sub_Tcell.Robj") 

#write.csv(merged@meta.data$sample_type,"sample_type.csv",quote=T)

marker<-FindMarkers(merged, group.by="sample_type",ident.1 = "Old_1", ident.2 = "Young_1", logfc.threshold = 0, min.pct = 0.1)
head(marker)
View(marker)
write.table(marker,"sample_typemarker.txt",sep="\t")
VlnPlot(merged, features = c("Cd4", "Cd8a"),split.by = "orig.ident")


DotPlot(merged, features = c("Cd74", "Itgax", "Cd3e", "Cd4", "Cd8a", "Ctla4","Pdcd1", "Foxp3", "Nkg7", "Itgam", "Cd3d"), cols = c("blue", "red"))
#CD8T STAT
DotPlot(sub_Tcell, features = c("Tcf7", "S1pr1","Lef1", "Sell", "Il7r","Xcl1","Cxcr3","Il2ra","Slamf6","Id3","Ccr7" ,"Prf1", "Fasl", "Gzmb", "Gzmk","Id2","Tox","Batf","Entpd1", "Ctla4","Pdcd1","Havcr2","Tigit","Lag3"), cols = c("blue", "red"))

DotPlot(sub_Tcell, features = c("Tcf7", "Ccl5", "Trbc2", "Cd69", "Bcl6", "Il7r", "Cd44", "Eomes", "Ccr7","Prdm1","Klrg1","Tbx21"), cols = c("blue", "red"))
DotPlot(sub_Tcell, features = c("Prdm1", "Il2ra", "Btla", "Il2rb", "Tnfrsf8", "Il7r", "Cxcr3", "Sell", "Ccr7","Cxcr6","Icos"), cols = c("blue", "red"))
DotPlot(sub_Tcell, features = c("Klrd1", "Ifngr1", "Cxcr3", "Batf", "Il17e", "Gata3", "Irf4", "Nrp1", "Izumo1r","Trbc2","Il2ra","Tcf7"), cols = c("blue", "red"))

VlnPlot(sub_Tcell, features = c("Cd4", "Cd8a"),split.by = "RNA_snn_res.1.2",pt.size = 0)
VlnPlot(sub_Tcell, features = c("Cd3e", "Cd3d"),split.by = "RNA_snn_res.1.2",pt.size = 0)
VlnPlot(sub_Tcell, features = c("Ctla4","Pdcd1"),split.by = "RNA_snn_res.1.2",pt.size = 0)
VlnPlot(sub_Tcell, features = c("Cxcr6","Tcf7"),split.by = "RNA_snn_res.1.2",pt.size = 0)
VlnPlot(sub_Tcell, features = c("Prf1","Gzmb"),split.by = "RNA_snn_res.1.2",pt.size = 0)

dim(sub_Tcell)
dev.off()

VlnPlot(sub_Tcell, features = c("Cd3e", "Cd3d","Cd4", "Cd8a"),split.by = "orig.ident",pt.size = 0)
VlnPlot(sub_Tcell, features = c("Cd3e", "Cd3d"),split.by = "orig.ident",pt.size = 0)


                                  ##==鉴定细胞类型==##
table(sub_Tcell$seurat_clusters)
library(patchwork)
library(SingleR)
library(tidyverse)
dir.create("CellType_T")
refdata <- ImmGenData()
testdata <- GetAssayData(sub_Tcell, slot="data")
table(sub_Tcell@meta.data$seurat_clusters)
table(sub_Tcell@meta.data$RNA_snn_res.1.2)
clusters <- sub_Tcell@meta.data$seurat_clusters
#使用ImmGen参考数据库鉴定
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"CellType_T/celltype_ImmGen.csv",row.names = F)
sub_Tcell@meta.data$celltype_ImmGen = "NA"
for(i in 1:nrow(celltype)){
  sub_Tcell@meta.data[which(sub_Tcell@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_ImmGen'] <- celltype$celltype[i]}
p1 = DimPlot(sub_Tcell, group.by="celltype_ImmGen", repel=T, label=T, label.size=5, reduction='tsne')
p1
p2 = DimPlot(merged, group.by="celltype_ImmGen", repel=T, label=T, label.size=5, reduction='umap')

p3 = p1+p2+ plot_layout(guides = 'collect')
ggsave("CellType/tSNE_celltype_ImmGen.png", p1, width=7 ,height=6)
ggsave("CellType/UMAP_celltype_ImmGen.png", p2, width=7 ,height=6)
ggsave("CellType/celltype_ImmGen.png", p3, width=10 ,height=5)
#使用MouseRNAseqData参考数据库鉴定
refdata <- MouseRNAseqData()
# load('~/database/SingleR_ref/ref_DICE_1561s.RData')
# refdata <- ref_DICE
testdata <- GetAssayData(merged, slot="data")
clusters <- merged@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"CellType/celltype_RNAseq.csv",row.names = F)
merged@meta.data$celltype_RNAseq = "NA"
for(i in 1:nrow(celltype)){
  merged@meta.data[which(merged@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_RNAseq'] <- celltype$celltype[i]}
p4 = DimPlot(merged, group.by="celltype_RNAseq", repel=T, label=T, label.size=5, reduction='tsne')
p4
p5 = DimPlot(scRNA, group.by="celltype_RNAseq", repel=T, label=T, label.size=5, reduction='umap')
p6 = p3+p4+ plot_layout(guides = 'collect')
ggsave("CellType/tSNE_celltype_RNAseq.png", p4, width=7 ,height=6)
ggsave("CellType/UMAP_celltype_RNAseq.png", p5, width=7 ,height=6)
ggsave("CellType/celltype_DICE.RNAseq", p6, width=10 ,height=5)
#对比两种数据库鉴定的结果
p8 = p1+p4
p8
ggsave("CellType/ImmGen_RNAseq.png", p8, width=12 ,height=5)
table(merged@meta.data$celltype_ImmGen)
##保存数据
saveRDS(merged,'merged.rds')
#DimPlot(sub_pbmc, reduction = "tsne",label = TRUE, pt.size = 1.5)

DotPlot(sub_merged, features = c("Ms4a1", "Itgax", "Cd3e", "Cd4", "Cd8a", "Ctla4","Pdcd1", "Foxp3", "Nkg7", "Itgam"), cols = c("blue", "red"))
FeaturePlot(sub_merged, features = c("Cd3e", "Cd3d", "Trbc2","Tnfrsf9", "Foxp3", "Ctla4", "Zap70", "Bcl11b", "Trbc1", "Il7r", "Cd8a","Cd4"), reduction = "umap",cols = c("gray", "red"))
##每个cluster按照marker规定细胞类型
new.cluster.ids <- c("CD8+T", "Macrophages_1", "Monocytes_1", "DC_1", "Macrophages_2", "Naive CD8+T", "Neutrophils", "Monocytes_2", "Granulocytes",
                     "NK", "CD4+ T", "DC_2", "Treg", "Actived B", "DC_3", "DC_4", "Exhausted CD8+T", "cytotoxic CD8+T",
                     "Macrophages_3", "Monocytes_3", "Macrophages_4", "B cells", "Basophils")
names(new.cluster.ids) <- levels(merged)
merged <- RenameIdents(merged, new.cluster.ids)
DimPlot(merged, reduction = "tsne", label = TRUE, pt.size = 1.5) + NoLegend()
write.table(Idents(merged),"cell-type.xls",sep="\t",quote = F)

##存储细胞归类结果
saveRDS(merged, file = "merged_final.rds")
setwd("")
merged <- readRDS("merged_final.rds")


















