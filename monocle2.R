##安装monocle
rm(list=ls())
##载入monocle�?
library(Seurat)
library(dplyr)
library(Matrix)
library(monocle)


##读入pbmc数据
setwd("")
load("data\\sub_Tcell.Robj")
table(sub_Tcell$origin)
table(sub_Tcell@meta.data$RNA_snn_res.0.52)
##取需要做轨迹�?2个亚�?
table(Idents(sub_Tcell))
sub_Tcell <- subset(sub_Tcell, idents=c("Cytotoxic_CD8 T", "Activated_CD8 T", "Exhausted_CD8 T"))
sub_Tcell
markers <- FindAllMarkers(sub_Tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top100 <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
top100
save(sub_Tcell,file="data\\sub_CD8.Robj") 
table(sub_Tcell$orig.ident)

##往monocle加载Seurat对象
seurat_object <- sub_Tcell
data <- as(as.matrix(seurat_object@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat_object@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)#构建需要输入的三个文件：data，meta.data,fd


#Construct monocle cds
cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size());
cds <- detectGenes(cds,  min_expr = 3)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))
print(head(pData(cds)))
save(cds,file="data\\cds_normal.Robj")


ordering_genes<-as.matrix(top100)

cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree',auto_param_selection = F,num_paths=20) 
cds <- orderCells(cds)
cds <- orderCells(cds,reverse=T)
save(cds,file="data\\orderCells.Robj")
write.table(pData(cds),file="data\\my_pseudotime.txt")

table(cds$origin)
head(pData(cds))
#不同指标的轨迹图
plot_cell_trajectory(cds, color_by = "orig.ident")
plot_cell_trajectory(cds, color_by = "RNA_snn_res.0.52")
plot_cell_trajectory(cds, color_by = "Pseudotime",show_branch_points=F)
plot_cell_trajectory(cds, color_by = "origin",show_branch_points=F)
plot_cell_trajectory(cds, color_by = "celltype")
plot_cell_trajectory(cds, color_by = "origin",show_branch_points=F)

#轨迹图分面显�?
p1 <- plot_cell_trajectory(cds, color_by = "orig.ident",show_branch_points=F) + facet_wrap(~orig.ident, nrow = 1)
p2 <- plot_cell_trajectory(cds, color_by = "origin",show_branch_points=F) + facet_wrap(~origin, nrow = 1)
plotc <- p1/p2
plotc
p3 <- plot_cell_trajectory(cds, color_by = "origin",show_branch_points=F) + facet_wrap(~celltype, nrow = 1)
p3

