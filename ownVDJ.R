#20200915  
rm(list =ls())
library(Seurat)


                                                   ######老年组#######
setwd("")
sample_old <- Read10X(data.dir = "")
#colnames(pbmc.data) <- paste('BC11_TUMOR1', colnames(pbmc.data), sep = '_')
##创建Seurat对象与数据过滤
sample_old <- CreateSeuratObject(counts = sample_old, project = "Old", min.cells = 3, min.features = 200)
sample_old

##计算每个细胞的线粒体基因转录本数的百分比（%）,使用[[ ]] 操作符存放到metadata中
sample_old[["percent.mt"]] <- PercentageFeatureSet(sample_old, pattern = "^mt-")
##展示基因及线粒体百分比
VlnPlot(sample_old, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),col=3, ncol = 3)

sample_old <- subset(sample_old, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & percent.mt < 10)
#sample_old
#table(sample_old$orig.ident)
##表达量数据标准化,LogNormalize的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000 )
sample_old <- NormalizeData(sample_old, normalization.method = "LogNormalize", scale.factor = 10000)
#merged <- NormalizeData(merged) 或者用默认的
sample_old
##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
sample_old<- FindVariableFeatures(sample_old, selection.method = "vst", nfeatures = 3000)
##而对所有基因进行标准化的方法如下：
all.genes <- rownames(sample_old)
sample_old <- ScaleData(sample_old, features = all.genes)
sample_old <- ScaleData(sample_old, vars.to.regress = "percent.mt")#耗时
##线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；
sample_old <- RunPCA(sample_old, features = VariableFeatures(object = sample_old))
ElbowPlot(sample_old, ndims = 50)
##非线性降维（UMAP/tSNE)基于PCA空间中的欧氏距离计算nearest neighbor graph,优化任意两个细胞间的距离权重（输入上一步得到的PC维数）。
sample_old <- FindNeighbors(sample_old, dims = 1:50)

##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
sample_old <- FindClusters(sample_old, resolution = 1)
sample_old <- RunUMAP(sample_old, dims = 1:50)
sample_old <- RunTSNE(sample_old, dims = 1:50)#耗时

sample_vdj<-sample_old
tcr <- read.csv("")
tcr[1:4,]
tcr <- tcr[!duplicated(tcr$barcode), ] #挑出productive=ture
tcr[1:4,]
names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"  #改一下名字raw_clonotype_id
tcr[1:4,]
table(tcr$v_gene)
clono <- read.csv("")   #读入克隆型信息
head(clono)
tcr <- merge(tcr, clono)#合并两种信息
tcr[1:4,]
rownames(tcr) <- tcr[,2]
clono_seurat <- AddMetaData(object=sample_vdj, metadata=tcr)   #VDJ整合到seurat之中
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

                                                  ######青年组#######
rm(list = ls())
setwd("")
sample_young <- Read10X(data.dir = "")
#colnames(pbmc.data) <- paste('BC11_TUMOR1', colnames(pbmc.data), sep = '_')
##创建Seurat对象与数据过滤
sample_young <- CreateSeuratObject(counts = sample_young, project = "young", min.cells = 3, min.features = 200)
sample_young

##计算每个细胞的线粒体基因转录本数的百分比（%）,使用[[ ]] 操作符存放到metadata中
sample_young[["percent.mt"]] <- PercentageFeatureSet(sample_young, pattern = "^mt-")
##展示基因及线粒体百分比
VlnPlot(sample_young, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),col=3, ncol = 3)
##过滤细胞：保留gene数大于200小于5000的细胞；目的是去掉空GEMs和1个GEMs包含2个以上细胞的数据；而保留线粒体基因的转录本数低于10%的细胞,为了过滤掉死细胞等低质量的细胞数据。
sample_young <- subset(sample_young, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & percent.mt < 10)
sample_young
#sample_old
#table(sample_old$orig.ident)
##表达量数据标准化,LogNormalize的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000 )
sample_young <- NormalizeData(sample_young, normalization.method = "LogNormalize", scale.factor = 10000)
#merged <- NormalizeData(merged) 或者用默认的
sample_young
##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
sample_young<- FindVariableFeatures(sample_young, selection.method = "vst", nfeatures = 3000)
##而对所有基因进行标准化的方法如下：
all.genes <- rownames(sample_young)
sample_young <- ScaleData(sample_young, features = all.genes)
sample_young <- ScaleData(sample_young, vars.to.regress = "percent.mt")#耗时
##线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；
sample_young <- RunPCA(sample_young, features = VariableFeatures(object = sample_young))
ElbowPlot(sample_young, ndims = 50)
##非线性降维（UMAP/tSNE)基于PCA空间中的欧氏距离计算nearest neighbor graph,优化任意两个细胞间的距离权重（输入上一步得到的PC维数）。
sample_young <- FindNeighbors(sample_young, dims = 1:50)

##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
sample_young <- FindClusters(sample_young, resolution = 1)
sample_young <- RunUMAP(sample_young, dims = 1:50)
sample_young <- RunTSNE(sample_young, dims = 1:50)#耗时

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
clono_seurat <- AddMetaData(object=sample_vdj, metadata=tcr)   #VDJ整合到seurat之中
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

                             #########两组间比较####
library(scRepertoire)
rm(list=ls())
setwd("")
csv1 <- read.csv("filtered_contig_annotations.csv", stringsAsFactors = F)
csv2 <- read.csv("filtered_contig_annotations.csv", stringsAsFactors = F)
contig_list1 <- list(csv1,csv2)


combined <- combineTCR(contig_list1, samples = 
                         c("TUY", "TUO"), ID = c("T", "T"), cells ="T-AB")
head(combined[[1]])

##展示每个样本的克隆型数量
#?quantContig
p1 <- quantContig(combined, cloneCall="gene+nt", scale = F)
p2 <- quantContig(combined, cloneCall="gene+nt", scale = T)
plotc = p1 + p2
plotc
#设置exportTable = T，则输出表格而非图形
quantContig(combined, cloneCall="gene+nt", scale = T, exportTable = T)



##CDR3序列的长度分布，“aa”代表统计氨基酸序列长度
p1 <- lengthContig(combined, cloneCall="aa", chains = "combined") 
p2 <- lengthContig(combined, cloneCall="aa", chains = "single")
plotc = p1 + p2
plotc
#ggsave('VDJ/lengthContig.png', plotc, width = 8, height = 4)


##对比两个样本的克隆型
p1 = compareClonotypes(combined, numbers = 10, samples = c("TUY_T", "TUO_T"), 
                       cloneCall="aa", graph = "alluvial")
p1
#ggsave('VDJ/compareClonotypes.png', p1, width = 8, height = 4)


##Clonal Space Homeostasis，这个不知道怎么翻译
clonalHomeostasis(combined, cloneCall = "aa")
#此分析将克隆型按其相对丰度分为rare, small, medium, large, hyperexpanded
#5大类，并统计各个类别的占比

##克隆型分类占比，与上一个分析相似，只是分类方法调整了
clonalProportion(combined, cloneCall = "aa") 

##Overlap Analysis，分析样本相似性
#clonalOverlap(combined, cloneCall = "aa", method = "morisita")
##克隆多样性指数
clonalDiversity(combined, cloneCall = "aa", group = "samples")

##VDJ与scRNA整合分析
seurat <- readRDS("merge\\data\\merged_tutorial.rds")
names(seurat@meta.data$sample_type)
table(seurat@meta.data$sample_type)
dim(seurat)
seurat@meta.data$sample_type
#查看UMAP图
DimPlot(seurat, label = T) #+ NoLegend()


#查看原始seurat的meta.data有哪些内容
names(seurat@meta.data)
# [1] "nCount_RNA"  "nFeature_RNA"  "integrated_snn_res.0.5"  "seurat_clusters"
# [5] "Patient"     "Type"          "RawBarcode"  
##将VDJ数据添加到seurat对象的meta.data中
combined$TUY_T$sample
seurat <- combineExpression(combined, seurat, cloneCall="aa", groupBy = "sample")
#查看添加VDJ数据后seurat的meta.data有哪些内容
names(seurat@meta.data)

##UMAP图展示cloneType的分布
seurat$cloneType <- factor(seurat$cloneType, levels = c("Hyperexpanded (100 < X <= 500)", 
                                                        "Large (20 < X <= 100)", "Medium (5 < X <= 20)", 
                                                        "Small (1 < X <= 5)", "Single (0 < X <= 1)", NA))
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
DimPlot(seurat, group.by = "cloneType") + 
  scale_color_manual(values = colorblind_vector(5), na.value="grey")


##UMAP图展示特定克隆型的分布
seurat <- highlightClonotypes(seurat, cloneCall= "aa", 
                              sequence = c("CAVNGGSQGNLIF_CSAEREDTDTQYF", "NA_CATSATLRVVAEKLFF"))
DimPlot(seurat, group.by = "highlight")


##桑基图展示特定克隆型的来源
alluvialClonotypes(seurat, cloneCall = "aa", 
                   y.axes = c("Patient", "cluster", "Type"), 
                   color = "CAVNGGSQGNLIF_CSAEREDTDTQYF") + 
  scale_fill_manual(values = c("grey", colorblind_vector(1)))


##桑基图展示克隆型在病例-cluster-组织类型之间的关系
alluvialClonotypes(seurat, cloneCall = "aa", 
                   y.axes = c("Patient", "cluster", "Type"), 
                   color = "cluster")

##按cluster分析克隆型
combined2 <- expression2List(seurat, group = "cluster")
p1 = clonalDiversity(combined2, cloneCall = "aa") + ggtitle("clonalDiversity")
p2 = clonalHomeostasis(combined2, cloneCall = "aa") + ggtitle("clonalHomeostasis")
p3 = clonalProportion(combined2, cloneCall = "aa") + ggtitle("clonalProportion")
p4 = clonalOverlap(combined2, cloneCall="aa", method="overlap")  + ggtitle("clonalOverlap")  
plotc = (p2|p3)/(p1|p4)
ggsave('VDJ/clonal.png', plotc, width = 10, height = 8)





