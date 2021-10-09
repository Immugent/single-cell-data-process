
rm(list=ls())
setwd('D:\\mouse_GSEA\\Cytotoxic_Tcell')

library(Seurat)
library(tidyverse)
dir.create("GSEA")   
dir.create("GSEA/input")
dir.create("GSEA/output")
scRNA <- readRDS("Cytotoxic_Tcell.rds")

tmp <- scRNA@meta.data
table(tmp$celltype)
tmp <- subset(tmp, subset = (tmp$celltype=='Old_Cytotoxic_CD8 T'|tmp$celltype=='Young_Cytotoxic_CD8 T'))
sub.cells <- rownames(tmp)
scRNAsub <- subset(scRNA, cells=sub.cells)
table(scRNAsub$cell_ontology_class)
scRNAsub$celltype <- gsub('Old_Cytotoxic_CD8 T','Old_Cytotoxic_CD8',scRNAsub$celltype)
scRNAsub$celltype <- gsub('Young_Cytotoxic_CD8 T','Young_Cytotoxic_CD8',scRNAsub$celltype)

expr <- GetAssayData(scRNAsub, slot = 'data')
head(expr)
expr <- data.frame(NAME=rownames(expr), Description=rep('na', nrow(expr)), expr, stringsAsFactors=F)
write('#1.2', "GSEA/input/expr.gct", ncolumns=1)
write(c(nrow(expr),(ncol(expr)-2)), "GSEA/input/expr.gct", ncolumns=2, append=T, sep='\t')
write.table(expr, "GSEA/input/expr.gct", row.names=F, sep='\t', append=T, quote=F)
line.1 <- c((ncol(expr)-2), 2, 1)
tmp <- table(as.character(scRNAsub@meta.data$celltype))
line.2 <- c("#", names(tmp))
line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2],tmp[2]))
write(line.1, 'GSEA/input/group.cls', ncolumns=length(line.1), append=T, sep='\t')
write(line.2, 'GSEA/input/group.cls', ncolumns=length(line.2), append=T, sep='\t')
write(line.3, 'GSEA/input/group.cls', ncolumns=length(line.3), append=T, sep='\t')



