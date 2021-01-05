#20200103 GSVA

library(GSVA)
library(Seurat)
rm(list=ls())
setwd("D:\\mouse_GSEA\\NK")
#saveRDS(sub_Treg, file = "sub_Treg.rds")
sub_T <- readRDS("subNK.rds")

#sub_Treg@assays$RNA@counts未标准化的数据
expr <- as.data.frame(sub_T@assays$RNA@counts)
head(expr)
expr=as.matrix(expr)

#构建参考数据集
library(msigdbr)

#小鼠所有的基因集
m_df = msigdbr(species = "Mus musculus")
##2.3 查看基因集类别：
a <- m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
View(a)
##2.5 检索鼠类C2 (curated) CGP (chemical and genetic perturbations)基因集：
m_df = msigdbr(species = "Mus musculus", category = "C2" )
#GSVA
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)
kegg2 <- gsva(expr, msigdbr_list, kcdf="Poisson",method = "gsva",parallel.sz=1)
head(kegg2)
setwd("D:\\mouse_GSEA\\NK\\01 diff")
write.csv(kegg2,"C2.csv")

#差异分析
library(limma)
rt <- kegg2
head(rt)
#rt <- read.csv("")
logFCcutoff=0.32  
adjPvalueCutoff=0.05
## 1.构建比较矩阵
type=c( rep("con",238),rep("treat",190) )
## 1.构建比较矩阵
design=model.matrix(~ type)
## 比较矩阵命名
colnames(design)=c("con", "treat")
##2.线性模型拟合
fit=lmFit(rt, design)
##3.贝叶斯检验
fit=eBayes(fit)
#全部差异分析结果
all=topTable(fit, coef="con", number=Inf,adjust.method="holm")
all=rbind(id=colnames(all),all)
write.table(all,file="C2all.txt",sep="\t",quote=F,col.names=F)
#有统计学意义的结果
diff <- topTable(fit, coef="con", number=Inf,
                 p.value=adjPvalueCutoff, adjust="holm", lfc=logFCcutoff)
diffName=row.names(diff)
diff=rbind(id=colnames(diff),diff)
write.csv(diff,file="C2diff.csv",quote=F,col.names=F)

hmExp=rt[diffName,]
hmExp
hmExp=rbind(id=colnames(hmExp),hmExp)

write.table(hmExp,file="heatmap.txt",sep="\t",quote=F,col.names=T,)