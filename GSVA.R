
library(GSVA)
library(Seurat)
rm(list=ls())
setwd("D:\\mouse_GSEA\\NK")
#saveRDS(sub_Treg, file = "sub_Treg.rds")
sub_T <- readRDS("subNK.rds")

#sub_Treg@assays$RNA@counts
expr <- as.data.frame(sub_T@assays$RNA@counts)
head(expr)
expr=as.matrix(expr)

library(msigdbr)

m_df = msigdbr(species = "Mus musculus")
a <- m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
View(a)
##C2 (curated) CGP (chemical and genetic perturbations)
m_df = msigdbr(species = "Mus musculus", category = "C2" )
#GSVA
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)
kegg2 <- gsva(expr, msigdbr_list, kcdf="Poisson",method = "gsva",parallel.sz=1)
head(kegg2)
setwd("D:\\mouse_GSEA\\NK\\01 diff")
write.csv(kegg2,"C2.csv")


library(limma)
rt <- kegg2
head(rt)
#rt <- read.csv("")
logFCcutoff=0.32  
adjPvalueCutoff=0.05
type=c( rep("con",238),rep("treat",190) )
design=model.matrix(~ type)
colnames(design)=c("con", "treat")
fit=lmFit(rt, design)
fit=eBayes(fit)
all=topTable(fit, coef="con", number=Inf,adjust.method="holm")
all=rbind(id=colnames(all),all)
write.table(all,file="C2all.txt",sep="\t",quote=F,col.names=F)
diff <- topTable(fit, coef="con", number=Inf,
                 p.value=adjPvalueCutoff, adjust="holm", lfc=logFCcutoff)
diffName=row.names(diff)
diff=rbind(id=colnames(diff),diff)
write.csv(diff,file="C2diff.csv",quote=F,col.names=F)

hmExp=rt[diffName,]
hmExp
hmExp=rbind(id=colnames(hmExp),hmExp)

write.table(hmExp,file="heatmap.txt",sep="\t",quote=F,col.names=T,)
