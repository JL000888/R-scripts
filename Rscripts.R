options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")


######## Expression profile detection################
library(affyPLM)
Data<-ReadAffy()
Pset<-fitPLM (Data)
image(Data[,1])
image(Pset,type="weights",which=1,main="Weights")
image(Pset,type="resids",which=1,main="Residuals")
image(Pset,type="sign.resids",which=1,main="Residuals.sign")


library(affyPLM)
library(RColorBrewer)
Data<-ReadAffy()
Pset<-fitPLM (Data)
colors<-brewer.pal(12,"Set3")
Mbox(Pset,col=colors,main="RLE",las=3) 
Mbox(Pset,ylim=c(-1,1),col=colors,main="RLE",las=3)

library(affyPLM)
library(RColorBrewer)
Pset<-fitPLM (Data)
colors<-brewer.pal(12,"Set3")
boxplot(Pset,col=colors,main="NUSE",las=3)


data.deg<-AffyRNAdeg(Data)
plotAffyRNAdeg(data.deg,col=colors)
legend("topleft",sampleNames(Data),col=colors,lwd=1,inset=0.05,cex=0.2)

library(affyPLM)
library(affy)
library(RColorBrewer)
colors<-brewer.pal(12,"Set3")
Data<-ReadAffy()
Pset<-fitPLM (Data)
#Pre RMA 
Mbox(Pset,ylim=c(-1,1),col=colors,main="RLE",las=3)
###############
sampleNames(Data)
N=length(Data)
#Post-RMA
eset.rma<-rma(Data)
####after nor#######
Mbox(eset.rma,ylim=c(-1,1),col=colors,main="RLE",las=3)
normal_exprs<-exprs(eset.rma)
probeid<-rownames(normal_exprs)
normal_exprs<-cbind(probeid,normal_exprs)
write.table(normal_exprs,file="All_probeid.txt",sep='\t',quote=F,row.names=F)

############Annotation
probe_exp<-read.csv("All_probeid.txt",header=T,sep="\t",row.names=1)
probe_exp<- as.matrix(probe_exp)
probeid_geneid<-read.csv("GPL6244.annotation.txt",header=T,sep="\t")
probe_name<-rownames(probe_exp)
loc<-match(probeid_geneid[,1],probe_name)
probe_exp<-probe_exp[loc,]
raw_geneid<-as.numeric(as.matrix(probeid_geneid[,3]))
index<-which(!is.na(raw_geneid))
geneid<-raw_geneid[index]
exp_matrix<-probe_exp[index,]
geneidfactor<-factor(geneid)
gene_exp_matrix<-apply(exp_matrix,2,function(x) tapply(x,geneidfactor,mean))
rownames(gene_exp_matrix)<-levels(geneidfactor)
geneid<-rownames(gene_exp_matrix)
gene_exp_matrix2<-cbind(geneid,gene_exp_matrix)
write.table(gene_exp_matrix2,file="Merge_probe.txt",sep='\t',quote=F,row.names=F)
####gene id  To gene symbol
loc<-match(rownames(gene_exp_matrix),probeid_geneid[,3])
rownames(gene_exp_matrix)=probeid_geneid[loc,2]
genesymbol<-rownames(gene_exp_matrix)
gene_exp_matrix3<-cbind(genesymbol,gene_exp_matrix)
write.table(gene_exp_matrix3,file="Merge_genesybmol.txt",sep='\t',quote=F,row.names=F)

############# NA value procession
library(impute)
gene_exp_matrix<-read.table("Merge_genesybmol.txt",header=T,sep="\t")
gene_exp_matrix <-gene_exp_matrix[!duplicated(gene_exp_matrix[,1]),] 
rownames(gene_exp_matrix) <- gene_exp_matrix[,1]
gene_exp_matrix1 <- gene_exp_matrix[,-1]
gene_exp_matrix2<-as.matrix(gene_exp_matrix1)
write.table(gene_exp_matrix2,file="Merge_genesybmol_genesybmol.txt",sep='\t',quote=T,row.names=T)
######KNN 
imputed_gene_exp<-impute.knn(gene_exp_matrix2,k=10,rowmax=0.5,colmax=0.8,maxp=3000,rng.seed=362436069)

GeneExp<-imputed_gene_exp$data
genesymbol<-rownames(GeneExp)
GeneExp<-cbind(genesymbol,GeneExp)
write.table(GeneExp,file="genesybmol_KNN0.txt",sep='\t',quote=T,row.names=T)

#######################################################


library(sva)
## Loading required package: mgcv
## Loading required package: nlme
## This is mgcv 1.8-33. For overview type 'help("mgcv-package")'.
## Loading required package: genefilter
## Loading required package: BiocParallel
library(cluster)
library(oompaBase)

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

# PCA plot
source("batchPCA.R")

# GEO microarray
GSE53408.expr <- read.csv("GSE53408_genesybmol.csv", row.names = 1)
GSE117261.expr <- read.csv("GSE117261genesybmol_KNN0.csv", row.names = 1)
GSE24988.expr <- read.csv("GSE24988_genesybmol_KNN0.csv", row.names = 1)
range(GSE53408.expr) 
range(GSE117261.expr) 
range(GSE24988.expr) 

blue <- "#2874C5"
yellow <- "#EABF00"
green <- "#008B8A"

comgene <- intersect(intersect(rownames(GSE24988.expr), rownames(GSE53408.expr)), rownames(GSE117261.expr))
combined.expr <- cbind.data.frame(GSE24988.expr[comgene,],
                                  GSE53408.expr[comgene,],
                                  GSE117261.expr[comgene,])

#PCA plot
batchPCA(indata = t(scale(t(combined.expr))),
         batch = rep(c("GSE24988","GSE53408","GSE117261"), times = c(ncol(GSE24988.expr),ncol(GSE53408.expr),ncol(GSE117261.expr))),
         fig.dir = ".",
         PCA.fig.title = "Raw PCA for combined expression profile",
         cols = c(blue, yellow, green),
         showID = F,
         cex = 1.0,
         showLegend = T) 
## Warning in if (class(groupi) == "matrix") {: the condition has length > 1 and
## only the first element will be used

## Warning in if (class(groupi) == "matrix") {: the condition has length > 1 and
## only the first element will be used

## Warning in if (class(groupi) == "matrix") {: the condition has length > 1 and
## only the first element will be used
range(combined.expr)


#######combat
batch <- data.frame(batch = rep(c("GSE24988","GSE53408","GSE117261"), times = c(ncol(GSE24988.expr),ncol(GSE53408.expr),ncol(GSE117261.expr))))
modcombat = model.matrix(~1, data = batch)
combined.expr.combat <- as.data.frame(ComBat(dat=as.matrix(combined.expr), batch=batch$batch, mod=modcombat))
## Found3batches
## Adjusting for0covariate(s) or covariate level(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
write.csv(combined.expr.combat[1:3,], "output_combined_expr.csv", quote = F)
#write.csv(combined.expr.combat, "output_combined_expr.csv", quote = F)

#PCA plot
batchPCA(indata = t(scale(t(combined.expr.combat))),
         batch = rep(c("GSE24988","GSE53408","GSE117261"), times = c(ncol(GSE24988.expr),ncol(GSE53408.expr),ncol(GSE117261.expr))),
         fig.dir = ".",
         PCA.fig.title = "Combat PCA for combined expression profile",
         cols = c(blue, yellow, green),
         showID = F,
         cex = 1.0,
         showLegend = T) 
## Warning in if (class(groupi) == "matrix") {: the condition has length > 1 and
## only the first element will be used

## Warning in if (class(groupi) == "matrix") {: the condition has length > 1 and
## only the first element will be used

## Warning in if (class(groupi) == "matrix") {: the condition has length > 1 and
## only the first element will be used
range(combined.expr.combat)

#############
library(affy)
library(limma)
##import phenotype data
phenoData = read.AnnotatedDataFrame('SVA_GroupControl.txt')
pheno = pData(phenoData)
#View(pheno)
rt=read.table("All_Exp_Group Contl.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

##sva--combat
library(sva)
library(pamr)
batch = pheno[,c('batch')]
batch
modcombat = model.matrix(~1,data=pheno)
combat_edata = ComBat(dat=data,batch=batch,
                      mod=modcombat,par.prior=T,
                      prior.plot=T)

write.table(combat_edata,file="SVA_Expdata_Control.txt",sep="\t")

#################################################################
normal_exprs<-read.table("SVA_Expdata_Control.txt",header=T,sep="\t")
tumor_exprs<-read.table("SVA_Expdata_PAH.txt",header=T,sep="\t")

probe_exprs<-merge(normal_exprs,tumor_exprs,by="ID")
write.table(probe_exprs,file="Exp_ALL.txt",sep='\t',quote=F,row.names=F)
#####################
##################### 58  Control   and  164 PAH 
#####################
#######PCA For (Pre_combat)
dat<-read.table("Exp_GroupPre.txt",header=TRUE,sep="\t",row.names = 1)
traits = read.table('SVA_GroupPre.txt',sep="\t",header= T)

library(RColorBrewer)
require(graphics)
heat.colors(6, alpha = 1)  
terrain.colors(6, alpha = 1) 
topo.colors(6, alpha = 1)  
cm.colors(6, alpha = 1) 


color = factor(traits$group,
               labels = c("#4C00FFFF", "#00E5FFFF", "#00FF4DFF", "#E6FF00FF" ,"#FFFF00FF", "#FFE0B3FF"),
               levels = c("GSE117261_PAH","GSE117261_Control",
                          "GSE53408_PAH","GSE53408_Control",
                          "GSE24988_PAH","GSE24988_Control"))

pca <- princomp(dat)
library(scatterplot3d)
scatterplot3d(pca$loadings[,1:3],main='PCA',color=color,type='p',
              highlight.3d=F,angle=60,grid=T,box=T,scale.y=1,
              cex.symbols=0.8,pch=16,col.grid='lightblue')
legend("topright",paste(c("GSE117261_PAH","GSE117261_Control",
                          "GSE53408_PAH","GSE53408_Control",
                          "GSE24988_PAH","GSE24988_Control"))
       ,fill=c("#4C00FFFF", "#00E5FFFF", "#00FF4DFF", "#E6FF00FF" ,"#FFFF00FF", "#FFE0B3FF"),box.col="grey")

library(scatterplot3d)
pdf('pre_Combat_pca_360_fine.pdf',onefile=TRUE,width=8,height=8)
diffangle <- function(ang){
  scatterplot3d(pca$loadings[,1:3],main='PCA',color=color,type='p',
                highlight.3d=F,angle=ang,grid=T,box=T,scale.y=1,
                cex.symbols=1.2,pch=16,col.grid='lightblue')
  legend("topright",paste(c("GSE117261_PAH","GSE117261_Control",
                            "GSE53408_PAH","GSE53408_Control",
                            "GSE24988_PAH","GSE24988_Control")),
         fill=c("#4C00FFFF", "#00E5FFFF", "#00FF4DFF", "#E6FF00FF" ,"#FFFF00FF", "#FFE0B3FF"),box.col="grey")
}
sapply(seq(-360,360,5),diffangle)
dev.off()
#####################
#####################
##PCA For (Post_combat)
dat<-read.table("Exp_ALL.txt",header=TRUE,sep="\t",row.names = 1)
traits = read.table('SVA_GroupPre.txt',sep="\t",header= T)

library(RColorBrewer)
require(graphics)
heat.colors(6, alpha = 1)  
terrain.colors(6, alpha = 1) 
topo.colors(6, alpha = 1)  
cm.colors(6, alpha = 1) 


color = factor(traits$group,
               labels = c("#4C00FFFF", "#00E5FFFF", "#00FF4DFF", "#E6FF00FF" ,"#FFFF00FF", "#FFE0B3FF"),
               levels = c("GSE117261_PAH","GSE117261_Control",
                          "GSE53408_PAH","GSE53408_Control",
                          "GSE24988_PAH","GSE24988_Control"))

pca <- princomp(dat)
library(scatterplot3d)
scatterplot3d(pca$loadings[,1:3],main='PCA',color=color,type='p',
              highlight.3d=F,angle=60,grid=T,box=T,scale.y=1,
              cex.symbols=0.8,pch=16,col.grid='lightblue')
legend("topright",paste(c("GSE117261_PAH","GSE117261_Control",
                          "GSE53408_PAH","GSE53408_Control",
                          "GSE24988_PAH","GSE24988_Control"))
       ,fill=c("#4C00FFFF", "#00E5FFFF", "#00FF4DFF", "#E6FF00FF" ,"#FFFF00FF", "#FFE0B3FF"),box.col="grey")

library(scatterplot3d)
pdf('post_Combat_pca_360_fine1.pdf',onefile=TRUE,width=8,height=8)
diffangle <- function(ang){
  scatterplot3d(pca$loadings[,1:3],main='PCA',color=color,type='p',
                highlight.3d=F,angle=ang,grid=T,box=T,scale.y=1,
                cex.symbols=1.2,pch=16,col.grid='lightblue')
  legend("topright",paste(c("GSE117261_PAH","GSE117261_Control",
                            "GSE53408_PAH","GSE53408_Control",
                            "GSE24988_PAH","GSE24988_Control")),
         fill=c("#4C00FFFF", "#00E5FFFF", "#00FF4DFF", "#E6FF00FF" ,"#FFFF00FF", "#FFE0B3FF"),box.col="grey")
}
sapply(seq(-360,360,5),diffangle)
dev.off()

#############################DEGs Detection ###############
#############################

library(limma)
rt<-read.table("Exp_ALL.txt",header=T,sep="\t",row.names=1)
#differential
class<-c(rep("control",58),rep("Is",164))   ## 58  Control   and  164 PAH 
design<-model.matrix(~factor(class))
colnames(design)<-c("control","Is")
fit<-lmFit(rt,design)
fit2<-eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=200000)
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)
############################################

#write table
diffLab<-allDiff[with(allDiff, ((logFC>1.0 |logFC<(-1.0)) & adj.P.Val<0.05)),]
write.table(diffLab,file="diffEXp_1.0.xls",sep="\t",quote=F)

diffExpLevel<-rt[rownames(diffLab),]
qvalue=allDiff[rownames(diffLab),]$adj.P.Val
diffExpQvalue=cbind(qvalue,diffExpLevel)
write.table(diffExpQvalue,file="diffExpLeve_1.0.xls",sep="\t",quote=F)

Upgene = allDiff[(allDiff$adj.P.Val < 0.05 & (allDiff$logFC>1.0)),]
write.table(Upgene, "Upgene1.0.xls",sep="\t",quote=F)
Downgene = allDiff[(allDiff$adj.P.Val < 0.05 & (allDiff$logFC<(-1.0))),]
write.table(Downgene, "Downgene1.0.xls",sep="\t",quote=F)


##########################GO Analysis###################
##########################

#install.packages('GOplot')
library(clusterProfiler)
library(AnnotationHub)
library(AnnotationDbi)
library(GOplot)
library(ggplot2)

library(clustree)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
#########################
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(cowplot)
library(scCATCH)
library(RColorBrewer)
require("SingleCellExperiment", quietly = T)
require("scater", quietly = T)
library(Seurat)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(tidyverse)
library(SeuratData)

#rm(list=ls())
options(stringsAsFactors = F) 
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
# library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)  #BiocManager::install("pathview",ask = F,update = F)


log2FC_cutoff = log2(10)
pvalue_cutoff = 0.05
padj_cutoff = 0.05

#

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
gsym.fc <- read.csv("limmaTab_GO.csv", as.is = T)
dim(gsym.fc)

#keytypes(org.Mm.eg.db)
gsym.id <- bitr(gsym.fc$SYMBOL, #
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = "org.Hs.eg.db") #Mm  Hs  Rn   org.Rn.eg.db
#####org.Hs.eg.db     org.Mm.eg.db
head(gsym.id)
#####ENTREZID_____foldchange
idvec <- gsym.id$ENTREZID
names(idvec) <- gsym.id$SYMBOL
gsym.fc$ENTREZID <- idvec[gsym.fc$SYMBOL]
head(gsym.fc)
###save
write.csv(gsym.fc[,c(3,2)], "DEGs_GO_input.csv", quote = F, row.names = F)

#####GO
ego_ALL <- enrichGO(gene          = gsym.id$SYMBOL,
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)
write.csv(ego_ALL, "DEGs_ego_ALL—GO.csv", row.names = T)

ego_CC <- enrichGO(gene          = gsym.id$SYMBOL,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
write.csv(ego_CC, "DEGs_ego_CC—GO.csv", row.names = T)

ego_MF <- enrichGO(gene          = gsym.id$SYMBOL,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
write.csv(ego_MF, "DEGs_ego_MF—GO.csv", row.names = T)

ego_BP <- enrichGO(gene          = gsym.id$SYMBOL,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) 
write.csv(ego_BP, "DEGs_ego_BP—GO.csv", row.names = T)


#ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
#ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
#ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)
p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("Biological process")
p_CC <- barplot(ego_CC,showCategory = 10) + ggtitle("Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("Molecular function")
plotc <- p_BP|p_CC|p_MF

ggsave("DEGs_Go—barplot.pdf", plotc, width=18 ,height=6)

p_BP1 <- dotplot(ego_BP,showCategory = 10) + ggtitle("Biological process")
p_CC1 <- dotplot(ego_CC,showCategory = 10) + ggtitle("Cellular component")
p_MF1 <- dotplot(ego_MF,showCategory = 10) + ggtitle("Molecular function")
plotc1 <- p_BP1|p_CC1|p_MF1
plotc1
ggsave("DEGs_Go—dotplot.pdf", plotc1, width=18 ,height=6)


id.fc <- read.csv("DEGs_GO_input.csv", as.is = T)
head(id.fc)
###https://www.genome.jp/kegg/catalog/org_list.html 
#### KEGG GO ####
kegg_enrich_results <- enrichKEGG(gene  = id.fc$ENTREZID,
                                  organism  = "hsa", 
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.2)
kegg_enrich_results <- DOSE::setReadable(kegg_enrich_results, 
                                         OrgDb="org.Hs.eg.db", 
                                         keyType='ENTREZID')#ENTREZID to gene Symbol
write.csv(kegg_enrich_results@result,'DEGs_KEGG_gene_up_enrichresults.csv') 
save(kegg_enrich_results, file ='DEGs_KEGG_gene_up_enrichresults.Rdata')

go_enrich_results <- enrichGO(gene = id.fc$ENTREZID,
                              OrgDb = "org.Hs.eg.db",
                              ont   = "ALL"  ,     #One of "BP", "MF"  "CC"  "ALL" 
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.2,
                              readable      = TRUE)
write.csv(go_enrich_results@result, 'GO_gene_BP_enrichresults.csv') 
save(go_enrich_results, file ='GO_gene_enrichresults.Rdata')

browseKEGG(kegg_enrich_results, 'map04137') 




kegg_enrich_results@result$Description[1:25] 
i=1 
select_pathway <- kegg_enrich_results@result$ID[i] 

select_pathway
pathview(gene.data     = id.fc$ENTREZID,
         pathway.id    = select_pathway,
         species       = 'hsa' ,     
         kegg.native   = F,# 
         new.signature = T, 
         limit         = list(gene=2, cpd=1) 
)

#ggsave("pathview_hsa04151.pdf", pp1, width=12 ,height=12)


pathview(gene.data     = id.fc$ENTREZID,
         pathway.id    = select_pathway,
         species       = 'hsa' ,     
         kegg.native   = F,
         new.signature = T, 
         limit= list(gene=2.5, cpd=1) 
)

### goplot : Gene Ontology is organized as a directed acyclic graph.
#gop <- goplot(kegg_enrich_results, showCategory = 10)

#ggsave(gop, filename = "KEGGplot.pdf", width=10, height=10)

p_kegg <- barplot(kegg_enrich_results,showCategory = 10) + ggtitle("KEGG Pathway")
ggsave("DEGs_p_kegg—barplot.pdf", p_kegg, width=6 ,height=6)

p_kegg1 <- dotplot(kegg_enrich_results,showCategory = 10) + ggtitle("KEGG Pathway")
ggsave("DEGs_p_kegg—dotplot.pdf", p_kegg1, width=6 ,height=6)
### cnetplot: Gene-Concept Network

genelist <- as.numeric(gsym.id[,2]) 
names(genelist) <- row.names(gsym.id)

cnetp1 <- cnetplot(go_enrich_results,  foldChange = genelist,
                   showCategory = 6,
                   colorEdge = T,
                   node_label = 'all',
                   color_category ='steelblue')
cnetp2 <- cnetplot(go_enrich_results,  foldChange = genelist,
                   showCategory = 6,
                   node_label = 'gene',
                   circular = T, 
                   colorEdge = T)
ggsave(cnetp1,filename ='DEGs_Gocnetplot.pdf', width =12,height =10)
ggsave(cnetp2,filename = 'DEGs_Gocnetplot_cir.pdf', width =15,height =10)


cnetp1 <- cnetplot(kegg_enrich_results,  foldChange = genelist,
                   showCategory = 6,
                   colorEdge = T,
                   node_label = 'all',
                   color_category ='steelblue')
cnetp2 <- cnetplot(kegg_enrich_results,  foldChange = genelist,
                   showCategory = 6,
                   node_label = 'gene',
                   circular = T, 
                   colorEdge = T)
ggsave(cnetp1,filename ='DEGs_Keggcnetplot.pdf', width =18,height =18)
ggsave(cnetp2,filename = 'DEGs_Keggcnetplot_cir.pdf', width =18,height =18)         

## upsetplot  # emphasizes the gene overlapping among different gene sets
upsetp <- upsetplot(go_enrich_results, n = 10)+
  theme(plot.margin=unit(c(1,1,1,1),'lines')) 
ggsave(upsetp, filename = 'DEGs_upsetplot.pdf', width=15, height=10)

## upsetplot  # emphasizes the gene overlapping among different gene sets
upsetp <- upsetplot(kegg_enrich_results, n = 10)+
  theme(plot.margin=unit(c(1,1,1,1),'lines')) 
ggsave(upsetp, filename = 'DEGs_upsetplot.pdf', width=15, height=10)

## Tree plot  #
pt <- pairwise_termsim(kegg_enrich_results)
treep <- treeplot(pt,
                  showCategory = 30)
ggsave(treep, filename = 'DEGs_treeplot.pdf', width=18, height=10)


############################GSEA######################################
############################
# install.packages("msigdbr")
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(aplot)

msigdbr_species()


# choose C2 BP gene sets
genesets <-  msigdbr(species = "Homo sapiens",
                     category = "C2")

# check
head(genesets,3)
# # A tibble: 3 x 18
# gs_cat gs_subcat gs_name    gene_symbol entrez_gene ensembl_gene human_gene_symb~ h
# #   gs_exact_source <chr>, gs_url <chr>, gs_description <chr>, taxon_id <int>,
# #   ortholog_sources <chr>, num_ortholog_sources <dbl>

# TERM2GENE
gmt <- genesets %>% dplyr::select(gs_name,gene_symbol)


# TERM2GENE
gmt <- genesets %>% dplyr::select(gs_name,gene_symbol)

# load gene
gene <- read.table('limmaTab_GSEA.txt',header = T) %>%
  arrange(desc(log2FoldChange))

# check
head(gene,3)
#       gene_name log2FoldChange
# 1         Ecscr       6.015527
# 2       Gm32341       5.962540
# 3 B130034C11Rik       5.841789

geneList <- gene$log2FoldChange
names(geneList) <- gene$gene_name

# check
head(geneList,3)
#    Ecscr       Gm32341 B130034C11Rik
# 6.015527      5.962540      5.841789

#GSEA 
# GO enrich
gseaRes <- GSEA(geneList = geneList,
                TERM2GENE = gmt,
                minGSSize    = 1,
                maxGSSize    = 500,
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                verbose      = FALSE)

# to dataframe
data_ga <- data.frame(gseaRes) %>%
  filter(pvalue < 0.05)
save(gseaRes, file ='C2 data_gsea.Rdata')

write.table(gseaRes@result,file="gse_C2 all1.txt",sep='\t',quote=T,row.names=T)
write.table(data_ga,file="gse_C2 sig1.txt",sep='\t',quote=T,row.names=T)

#########################C5

genesets <-  msigdbr(species = "Homo sapiens",
                     category = "C5"#,subcategory = "BP"
)
# check
head(genesets,3)
# # A tibble: 3 x 18
# gs_cat gs_subcat gs_name    gene_symbol entrez_gene ensembl_gene human_gene_symb~ h
# #   gs_exact_source <chr>, gs_url <chr>, gs_description <chr>, taxon_id <int>,
# #   ortholog_sources <chr>, num_ortholog_sources <dbl>

# TERM2GENE
gmt <- genesets %>% dplyr::select(gs_name,gene_symbol)


# TERM2GENE
gmt <- genesets %>% dplyr::select(gs_name,gene_symbol)

# load gene
gene <- read.table('limmaTab_GSEA.txt',header = T) %>%
  arrange(desc(log2FoldChange))

# check
head(gene,3)
#       gene_name log2FoldChange
# 1         Ecscr       6.015527
# 2       Gm32341       5.962540
# 3 B130034C11Rik       5.841789

geneList <- gene$log2FoldChange
names(geneList) <- gene$gene_name

# check
head(geneList,3)
#    Ecscr       Gm32341 B130034C11Rik
# 6.015527      5.962540      5.841789

#GSEA 
# GO enrich
gseaRes <- GSEA(geneList = geneList,
                TERM2GENE = gmt,
                minGSSize    = 1,
                maxGSSize    = 500,
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                verbose      = FALSE)

# to dataframe
data_ga <- data.frame(gseaRes) %>%
  filter(pvalue < 0.05)
save(gseaRes, file ='C5 data_gsea.Rdata')

write.table(gseaRes@result,file="gse_C5 all.txt",sep='\t',quote=T,row.names=T)
write.table(data_ga,file="gse_C5 sig.txt",sep='\t',quote=T,row.names=T)



#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################

# gseaplot2
gseaplot2(gseaRes, geneSetID = 'GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN',
          title = gseaRes$Description['GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN'])

gseaScores <- getFromNamespace("gseaScores", "DOSE")

# define function
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}

# get data
gsdata <- gsInfo(gseaRes, geneSetID = 'GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN')


gsdata1 <- gsdata %>%
  mutate("gene_name" = gene$gene_name) %>%
  filter(position == 1)

# check
head(gsdata1,3)
#     x runningScore position        ymin       ymax geneList                                   Description gene_name
# 1   4   0.05996027        1 -0.02241342 0.02241342 5.803995 GOBP_NUCLEOSIDE_DIPHOSPHATE_METABOLIC_PROCESS     Hkdc1
# 2  29   0.11168581        1 -0.02241342 0.02241342 5.081191 GOBP_NUCLEOSIDE_DIPHOSPHATE_METABOLIC_PROCESS    Entpd3
# 3 122   0.15419642        1 -0.02241342 0.02241342 4.426756 GOBP_NUCLEOSIDE_DIPHOSPHATE_METABOLIC_PROCESS      Dlg2


# colnames
colnames(gsdata1)
# [1] "x"            "runningScore" "position"     "ymin"         "ymax"         "geneList"
# [7] "Description"  "gene_name"



##############curve plot

# plot
pcurve <- ggplot(gsdata,aes(x = x,y = runningScore,
                            color = runningScore)) +
  geom_hline(yintercept = 0,size = 1,color = 'black',
             lty = 'dashed') +
  geom_line() +
  # geom_segment(data = gsdata1,aes(xend = x,yend = 0)) +
  theme_bw(base_size = 14) +
  scale_color_gradient(low = '#76BA99',high = '#EB4747') +
  scale_x_continuous(expand = c(0,0)) +
  # scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = 'none',
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        plot.margin = margin(t = .2,r = .2, b = 0,l = .2,unit = "cm")) +
  ylab('Running Enrichment Score')

pcurve

#egment plot

pseg <- ggplot(gsdata,aes(x = x,y = runningScore)) +
  geom_segment(data = gsdata1,
               aes(x = x,xend = x,y = 0,yend = 1),
               color = 'black',show.legend = F) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw(base_size = 14) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = margin(t = 0,r = .2,b = .2,l = .2,unit = "cm")) +
  xlab('Rank in Ordered Dataset')

pseg

#segment and heatmap

v <- seq(1, sum(gsdata$position), length.out = 9)
inv <- findInterval(rev(cumsum(gsdata$position)), v)
if (min(inv) == 0) inv <- inv + 1

col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))

# ymin <- min(p2$data$ymin)
# yy <- max(p2$data$ymax - p2$data$ymin) * .3
ymin <- 0
yy <- 0.3
xmin <- which(!duplicated(inv))
xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
d <- data.frame(ymin = ymin, ymax = yy,
                xmin = xmin,
                xmax = xmax,
                col = col[unique(inv)])

pseg_ht <- pseg + geom_rect(
  aes_(xmin = ~xmin,xmax = ~xmax,
       ymin = ~ymin,ymax = ~ymax,
       fill = ~I(col)),
  data = d,
  alpha = 0.8,
  inherit.aes = FALSE)

pseg_ht



#gene rank plot

# add gene rank
pseg_ht1 <- pseg_ht + xlab('') +
  theme(axis.title.x = element_blank(),
        plot.margin = margin(t = -.1,r = .2,b = 0,l = .2,unit = "cm"))

prank <- ggplot(gsdata,aes(x = x,y = geneList)) +
  geom_col(width = 1,fill = 'grey80',color = NA) +
  # geom_col(aes(fill = geneList),
  #          width = 1,color = NA,show.legend = F) +
  # scale_fill_gradient2(low = col[1],mid = 'white',high = col[length(col)],midpoint = 0) +
  geom_hline(yintercept = 0,size = 0.8,color = 'black',
             lty = 'dashed') +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        plot.margin = margin(t = -.1,r = .2,b = .2,l = .2,unit = "cm")) +
  coord_cartesian(expand = 0) +
  ylab('Ranked List') +
  xlab('Rank in Ordered Dataset')

prank


#拼图

# combine
pall <- aplot::plot_list(gglist = list(pcurve,pseg_ht1,prank),
                         ncol = 1, heights = c(0.5,0.2,0.3))

pall

# add gene name
geneLabel <- gsdata1 %>% arrange(desc(runningScore)) %>%
  head(20)

plabel <- pcurve +
  geom_segment(data = geneLabel,aes(xend = x,yend = 0),
               color = 'red') +
  geom_text_repel(data = geneLabel,
                  aes(label = gene_name),
                  force = 20,
                  max.overlaps = 50,
                  # nudge_y = 0.2,
                  size = 4,
                  fontface = 'italic')

aplot::plot_list(gglist = list(plabel,pseg_ht1,prank),
                 ncol = 1, heights = c(0.5,0.2,0.3))


panother <- ggplot(gsdata,aes(x = x,y = runningScore,color = runningScore)) +
  geom_hline(yintercept = 0,size = 0.8,color = 'black',
             lty = 'dashed') +
  geom_point() +
  geom_line() +
  geom_segment(data = gsdata1,aes(xend = x,yend = 0)) +
  theme_bw(base_size = 14) +
  # scale_color_gradient(low = '#336699',high = '#993399') +
  scale_color_gradient2(low = '#336699',mid = 'white',high = '#993399',midpoint = 0.2) +
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.position = 'none',
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        plot.margin = margin(t = .2,r = .2, b = 0,l = .2,unit = "cm")) +
  ylab('Running Enrichment Score')

panother

# add gene name
panother_label <-
  panother +
  geom_text_repel(data = geneLabel,
                  aes(label = gene_name),
                  force = 20,
                  max.overlaps = 50,
                  # nudge_y = 0.15,
                  size = 4,
                  # color = 'black',
                  fontface = 'italic')

panother_label


# new color
color <- rev(colorRampPalette(c("#336699","white", "#993399"))(10))

ht <- ggplot(gsdata,aes(x = x,y = runningScore)) +
  geom_rect(aes_(xmin = ~xmin,xmax = ~xmax,
                 ymin = ~ymin,ymax = ~ymax,
                 fill = ~I(color)),
            data = d,
            alpha = 0.8,
            inherit.aes = FALSE) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(t = 0,r = .2, b = .2,l = .2,unit = "cm"))

# combine
aplot::plot_list(gglist = list(panother_label,ht),
                 ncol = 1, heights = c(0.9,0.1))


###################################################################################
###################################################################################
###################################################################################

library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
# options(mc.cores = parallel::detectCores())

msigdbr_show_species()
##  [1] "Bos taurus"               "Caenorhabditis elegans"  
##  [3] "Canis lupus familiaris"   "Danio rerio"             
##  [5] "Drosophila melanogaster"  "Gallus gallus"           
##  [7] "Homo sapiens"             "Mus musculus"            
##  [9] "Rattus norvegicus"        "Saccharomyces cerevisiae"
## [11] "Sus scrofa"
#http://software.broadinstitute.org/gsea/downloads.jsp
h <- msigdbr(species = "Homo sapiens", 
             category = "C2")### C2 C5  C7
######################################################
h <- select(h, gs_name, gene_symbol) %>% 
  as.data.frame %>% 
  split(., .$gs_name) %>% 
  lapply(., function(x)(x$gene_symbol))

gs <- lapply(h, unique)

count <- table(unlist(gs))
keep <- names(which(table(unlist(gs)) < 20000))
gs <- lapply(gs, function(x) intersect(keep, x))
gs <- gs[lapply(gs, length) > 0]
head(gs)

save(gs, file = "C2.gs.RData")
###C2#################GSVA############################

(load("C2.gs.RData")) #
## [1] "gs"

gsym.expr <- read.csv("Exp_ALL.csv",row.names = 1)
head(gsym.expr)
##############GSVA Score########################
gsva_es <- gsva(as.matrix(gsym.expr), gs)

head(gsva_es)
################################ͨ·Score#####################
write.csv(gsva_es, "gsva_output_C2.gs.csv", quote = F)
#############################################################
####@@@@@@@@@@@@  TernaryCluster  

options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")


Packages <- c("pheatmap","gplots")
for(pk in Packages){
  if(!require(pk,character.only = T,quietly = T)) install.packages(pk)
  suppressMessages(library(pk,character.only = T))
}
## 
## Attaching package: 'gplots'
## The following object is masked from 'package:stats':
## 
##     lowess
if(!require("GSVA",character.only = T,quietly = T)){
  if(!require("BiocManager",character.only = T, quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("GSVA")
} 
suppressMessages(library("GSVA",character.only = T))
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

tumornumber <- 58 
normalnumber <- 164 
signaturename <- "ELECTRON_TRANSPORT" 

halfwidth <- 0.025

ClusterNumber <- 4

ClusterColor = c("#B31E22","#529B40","#020105","#383D8E")


expr <- read.csv("Exp_ALL.csv",check.names = F,row.names = 1)
expr[1:3,1:3]
##        LIHC-A114-01A LIHC-A4NB-01A LIHC-A5RG-01A
## TSPAN6         41.24         28.55         59.39
## TNMD            0.03          0.00          0.02
## DPM1           17.56         21.47         28.95
signature <- read.csv("ELECTRON_TRANSPORT.csv")[,1]
signature
##  [1] "SMAD3"  "SPTBN1" "SMAD4"  "SMAD5"  "SMAD7"  "TGFBR1" "TGFBR2"
##  [8] "ACVR2A" "EP300"  "CREBBP" "RUNX2"  "CTCF"   "DLG1"   "BMPR1A"
## [15] "BMPR2"  "TGFB1"  "TGFB2"  "TGFB3"

siglist <- list()
siglist[[signaturename]] <- signature
signature.gsva <- gsva(as.matrix(expr[,1:tumornumber]),siglist,method = "gsva")
## Estimating GSVA scores for 1 gene sets.
## Computing observed enrichment scores
## Estimating ECDFs with Gaussian kernels
## Using parallel with 4 cores

rowids <- intersect(rownames(expr), signature)
logdata <- log2(expr[rowids, ] + 0.5)
tumordata <- logdata[, 1:tumornumber] 
var <- apply(tumordata, 1, sd, na.rm=T) 
selvar <- var[var>0] 
tumordata <- tumordata[names(selvar), ]
normaldata <- logdata[names(selvar), (tumornumber+1):(tumornumber+normalnumber)] 

halfwidth <- 0.025
normaldown <- apply(normaldata, 1, function(x) quantile(x, probs=halfwidth, na.rm=T) ) 
normalup <- apply(normaldata, 1, function(x) quantile(x, probs=1-halfwidth, na.rm=T) )


for (k in 1:nrow(tumordata)) {
  rowk <- as.numeric(tumordata[k, ])
  out <- rep(0, times=ncol(tumordata)) 
  out[rowk>normalup[k]] <- 1 
  out[rowk<normaldown[k]] <- -1 
  tumordata[k, ] <- out
}


outdata <- tumordata
outdata[outdata==1] <- "UP" 
outdata[outdata==-1] <- "DOWN" 
outdata[outdata==0] <- "NOCHANGE" 
write.csv(outdata,"ELECTRON_TRANSPORT_Ternary.csv",row.names = T,col.names = NA)


hcg <- hclust(dist(tumordata), "ward.D")
hcs <- hclust(dist(t(tumordata)), "ward.D")
k <- ClusterNumber 
group <- paste0("C",cutree(hcs,k))
names(group) <- colnames(tumordata)


annCol <- data.frame(Cluster=group,score=signature.gsva[signaturename,])

p <- kruskal.test(annCol$score ~ as.factor(annCol$Cluster))$p.value
colnames(annCol)[2] <- paste0(signaturename,"_score P = ",format(p,digits = 3)) 


annColors <- list(Cluster = c(C1 = ClusterColor[1],
                              C2 = ClusterColor[2],
                              C3 = ClusterColor[3],
                              C4 = ClusterColor[4]),
                  bluered(64))
names(annColors)[2] <- colnames(annCol)[2]

pct <- paste0(round(rowSums(tumordata == 1)/ncol(tumordata),2) * 100,"%") 
rownames(tumordata) <- paste0(rownames(tumordata)," (",pct,")") 


pheatmap(tumordata,
         cluster_rows = hcg,
         cluster_cols = hcs,
         color = c("#283285","#DEDEDE","#651430"),
         annotation_col = annCol,
         annotation_colors = annColors,
         fontsize_col = 2, 
         fontsize_row = 8, 
         fontsize = 6, 
         cutree_cols = k, 
         legend_breaks = c(-1,0,1), 
         legend_labels = c("mRNA\ndown-regulation","mRNA\nno-change","mRNA\nup-regulation"), 
         filename = "Heatmap_TernaryCluster_ELECTRON_TRANSPORT.pdf",width = 8,height = 8)


########################## GSVA score 
###########################################################################################
normal_exprs<-read.table("AllSampleGroup.txt",header=T,sep="\t")
tumor_exprs<-read.table("Mito_PathwayScore.txt",header=T,sep="\t")

probe_exprs<-merge(normal_exprs,tumor_exprs,by="ID")
write.table(probe_exprs,file="Mito_PathwayScore_Group.txt",sep='\t',quote=F,row.names=F)


knitr::opts_chunk$set(echo = TRUE)
set.seed(987654321)
library(knitr)
library(kableExtra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(janitor)
library(ggalluvial)
library(xgboost)
#install.packages('rBayesianOptimization')
library(rBayesianOptimization)
#install.packages('Amelia')
library(Amelia)
library(patchwork)
#install.packages('tidyquant')
library(SHAPforxgboost)
library(tidyquant)
library(tidyverse)
library(caret)
library(PRROC)
#remotes::install_github("AppliedDataSciencePartners/xgboostExplainer")
library(xgboostExplainer)
library(viridis)
library(xlsx)

protocol_fill_color = "grey25"

theme_bluewhite <- function (base_size = 11, base_family = "serif") {
  theme_bw() %+replace% 
    theme(
      text = element_text(family = "serif"),
      panel.grid.major  = element_line(color = "white"),
      panel.background = element_rect(fill = "grey97"),
      panel.border = element_rect(color = "darkred", fill = NA, size = 1), ##05014a
      axis.line = element_line(color = "grey97"),
      axis.ticks = element_line(color = "grey25"),
      axis.title = element_text(size = 10),
      axis.text = element_text(color = "grey25", size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 15, hjust = 0.5),
      strip.background = element_rect(fill = '#05014a'),
      strip.text = element_text(size = 10, colour = 'white'), # changes the facet wrap text size and color
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}


metadata <- read.csv("ELECTRON_TRANSPORT_EXPModel.csv")
#metadata<-heart
patientFeatures <- metadata  %>% 
  dplyr::select(-c(X)) %>% 
  mutate(
    gender = gender - 1
  ) %>% 
  mutate_if(is.integer, as.numeric) %>% 
  mutate(
    gender = factor(gender)
  )

patientFeaturesWithMissing <- patientFeatures

patientFeatures <- patientFeatures %>% 
  discard(~sum(is.na(.x))/length(.x)* 100 >=50)


#################################
#####Train test##############
#################################

smp_size <- floor(0.75 * nrow(patientFeatures))
train_ind <- sample(seq_len(nrow(patientFeatures)), size = smp_size)

train <- patientFeatures[train_ind, ]
test <- patientFeatures[-train_ind, ]

X_train <- train %>% 
  mutate(gender = as.numeric(gender)) %>% 
  dplyr::select(-c(PATIENT_ID, Admission.date, Discharge.date, outcome)) %>% 
  as.matrix()

Y_train <- train %>% 
  dplyr::select(c(outcome)) %>% 
  as.matrix()

X_test <- test %>% 
  mutate(gender = as.numeric(gender)) %>% 
  dplyr::select(-c(PATIENT_ID, Admission.date, Discharge.date, outcome)) %>% 
  as.matrix()

Y_test <- test %>% 
  dplyr::select(c(outcome)) %>% 
  as.matrix()


MultiMLDataSet <- bind_cols(data.frame(Y_train), data.frame(X_train))
#####################################################################
#########################Logistic regression####################

set.seed(123)
Logistic_Model <- glm(outcome ~ ., data = na.omit(MultiMLDataSet), family = "binomial")

Logistic_Predictions_value = predict(Logistic_Model, data.frame(X_test), type = "response")
Logistic_Predictions = Logistic_Predictions_value
Logistic_Predictions_binomial=ifelse(Logistic_Predictions>0.5,1,0)


library(caret)
con_Logistic <- confusionMatrix(factor(Logistic_Predictions_binomial), factor(Y_test))

MCC_Logistic <- list(
  TP <- con_Logistic$table[[1]],
  FP <- con_Logistic$table[[2]],
  FN <- con_Logistic$table[[3]],
  TN <- con_Logistic$table[[4]]
) %>%  # # MCC <- ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP+FN) * (TN+FP) * (TN+FN))
  pmap_dbl(., ~ ((..1 * ..4) - (..2 * ..3))/sqrt((..1 + ..2) * (..1 + ..3) * (..4 + ..2) * (..4 + ..3)))


Logistic_pred_outcome <- cbind(as.numeric(Logistic_Predictions), as.numeric(Y_test)) %>% 
  data.frame() %>% 
  setNames(c("predictions", "outcome"))


### PR curve
Logistic_fg <- Logistic_pred_outcome %>% 
  filter(outcome == 1) %>%
  drop_na() %>% 
  pull(predictions)

Logistic_bg <- Logistic_pred_outcome %>% 
  filter(outcome == 0) %>%
  drop_na() %>% 
  pull(predictions)

Logistic_pr <- PRROC::pr.curve(
  scores.class0 = Logistic_fg, scores.class1 = Logistic_bg, curve = TRUE)

Logistic_pr_curve <- Logistic_pr$curve %>% 
  data.frame() %>% 
  mutate(X4 = "Logistic Regression")

library(pROC)
Logistic_roc <- roc(as.numeric(Y_test),as.numeric(Logistic_Predictions))


Logistic_pr_Plot <- 
  ggplot(Logistic_pr_curve,aes(x = X1, y = X2, color = X4)) +
  geom_line() +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = c("Logistic Regression")) +
  labs(title= paste0("Precision-Recall Curves, AUC=",round(Logistic_pr$auc.integral,3)), 
       y = "Precision",
       x = "Recall") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0, 1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )

Logistic_pr_Plot

######################### ROC Curves ##################################

Logistic_roc_Plot <- 
  ggroc(Logistic_roc,legacy.axes = TRUE, linetype = 1) +
  geom_abline(show.legend = TRUE, alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = c("Logistic Regression")) +
  labs(title=  paste0("Precision-Recall Curves, AUC=",round(Logistic_roc$auc,3)), 
       y = "Sensitivity",
       x = "1-Specificity") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "left",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )
Logistic_roc_Plot

##########################################################################################
(Logistic_roc_Plot + Logistic_pr_Plot) + plot_annotation(title = "Algorithm Comparison - Model Performance", tag_levels = "A", theme = theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")))

ggsave("./Logistic_rocPRPlot.pdf",width = 10,height = 7)

write.table(MultiMLDataSet,file="MultiMLDataSet.txt",quote=T,sep="\t",col.names=T)
write.table(Logistic_Model$coefficients,file="Logistic_Model_coefficients.txt",quote=T,sep="\t",col.names=T)


#####################################################################
########################Naive Bayes###########################


MultiMLDataSet <- bind_cols(data.frame(Y_train), data.frame(X_train))

### 2.train
set.seed(123)
Naive_Bayes_Model <- e1071::naiveBayes(as.factor(outcome) ~ ., data = MultiMLDataSet)
Naive_Bayes_Model

### 3. test
NB_Predictions = predict(Naive_Bayes_Model, X_test, type = "class")


library(caret)

con_NB <- confusionMatrix(NB_Predictions, factor(Y_test))

MCC_NB <- list(
  TP <- con_NB$table[[1]],
  FP <- con_NB$table[[2]],
  FN <- con_NB$table[[3]],
  TN <- con_NB$table[[4]]
) %>%  # # MCC <- ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP+FN) * (TN+FP) * (TN+FN))
  pmap_dbl(., ~ ((..1 * ..4) - (..2 * ..3))/sqrt((..1 + ..2) * (..1 + ..3) * (..4 + ..2) * (..4 + ..3)))

NB_pred_outcome <- cbind(as.numeric(NB_Predictions), as.numeric(Y_test)) %>% 
  data.frame() %>% 
  setNames(c("predictions", "outcome"))

NB_fg <- NB_pred_outcome %>% 
  filter(outcome == 1) %>%
  pull(predictions)

NB_bg <- NB_pred_outcome %>% 
  filter(outcome == 0) %>%
  pull(predictions)

NB_pr <- PRROC::pr.curve(
  scores.class0 = NB_fg, scores.class1 = NB_bg, curve = TRUE)
plot(NB_pr)


NB_pr_curve <- NB_pr$curve %>% 
  data.frame() %>% 
  mutate(X4 = "Naive Bayes")

library(pROC)
NB_roc <- roc(as.numeric(Y_test),as.numeric(NB_Predictions))
NB_pr_plot <- ggplot(NB_pr_curve,aes(x = X1, y = X2, color = X4)) +
  geom_line() +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = c("Naïve Bayes")) +
  labs(title= paste0("Precision-Recall Curves, AUC=",round(NB_pr$auc.integral,3)),  
       y = "Precision",
       x = "Recall") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0, 1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )
NB_pr_plot



NB_roc_Plot <- 
  ggroc(NB_roc,legacy.axes = TRUE, linetype = 1) +
  geom_abline(show.legend = TRUE, alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = c("Naïve Bayes")) +
  labs(title= paste0("Precision-Recall Curves, AUC=",round(NB_roc$auc,3)), 
       y = "Sensitivity",
       x = "1-Specificity") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "left",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )
NB_roc_Plot
###########################################################################

(NB_roc_Plot + NB_pr_plot) + plot_annotation(title = "Algorithm Comparison - Model Performance", tag_levels = "A", theme = theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")))

ggsave("./NB_rocPRPlot.pdf")




#####################################################################
#####################################################################
#####################################################################
###################### decision tree ######################### 

library(rpart)

classTree <- rpart(factor(outcome) ~ ., data = MultiMLDataSet, method = "class")
print(classTree)
plotcp(classTree)
summary(classTree)

classTreePredictions <- predict(classTree, data.frame(X_test), type = "class")


con_classTree <- confusionMatrix(classTreePredictions, factor(Y_test))

MCC_classTree <- list(
  TP <- con_classTree$table[[1]],
  FP <- con_classTree$table[[2]],
  FN <- con_classTree$table[[3]],
  TN <- con_classTree$table[[4]]
) %>%  # # MCC <- ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP+FN) * (TN+FP) * (TN+FN))
  pmap_dbl(., ~ ((..1 * ..4) - (..2 * ..3))/sqrt((..1 + ..2) * (..1 + ..3) * (..4 + ..2) * (..4 + ..3)))

classTree_pred_outcome <- cbind(as.numeric(classTreePredictions), as.numeric(Y_test)) %>% 
  data.frame() %>% 
  setNames(c("predictions", "outcome"))

classTree_fg <- classTree_pred_outcome %>% 
  filter(outcome == 1) %>%
  pull(predictions)

classTree_bg <- classTree_pred_outcome %>% 
  filter(outcome == 0) %>%
  pull(predictions)

classTree_pr <- PRROC::pr.curve(
  scores.class0 = classTree_fg, scores.class1 = classTree_bg, curve = TRUE)


classTree_pr_curve <- classTree_pr$curve %>% 
  data.frame() %>% 
  mutate(X4 = "Classification Tree")

classTree_prPlot <-
  ggplot(classTree_pr_curve,aes(x = X1, y = X2, color = X4)) +
  geom_line() +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = c("classTree")) +
  labs(title= paste0("Precision-Recall Curves, AUC=",round(classTree_pr$auc.integral,3)),  
       y = "Precision",
       x = "Recall") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0, 1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )

library(pROC)
classTree_roc <- roc(as.numeric(Y_test),as.numeric(classTreePredictions))



classTree_rocPlot <- 
  ggroc(classTree_roc,legacy.axes = TRUE, linetype = 1) +
  geom_abline(show.legend = TRUE, alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = "classTree") +
  labs(title= paste0("ROC Curves, AUC=",round(classTree_roc$auc,3)),  
       y = "Sensitivity",
       x = "1-Specificity") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "left",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )

###########################################################################

(classTree_rocPlot + classTree_prPlot) + plot_annotation(title = "Algorithm Comparison - Model Performance", tag_levels = "A", theme = theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")))

ggsave("classTree_rocPRPlot.pdf")

#install.packages('rpart.plot')
library(rpart.plot)
dev.off()
rpart.plot(classTree, # middle graph
           extra = 104, # show fitted class, probs, percentages
           box.palette = "GnBu", # color scheme
           branch.lty = 3, # dotted branch lines
           shadow.col = "white", # shadows under the node boxes
           nn = TRUE,
           yesno = 2, type = 1,
           #branch = .9,
           clip.right.labs = TRUE
) # display the node number



#####################################################################
###########################random forest##################################


#library(randomForest) 
Random_Forest_Model <- randomForest::randomForest(factor(outcome) ~., data = MultiMLDataSet, 
                                                  na.action = na.omit, ntree = 500, importance = TRUE)

RF_Predictions <- predict(Random_Forest_Model, X_test)


con_RF <- confusionMatrix(RF_Predictions, factor(Y_test))
MCC_RF <- list(
  TP <- con_RF$table[[1]],
  FP <- con_RF$table[[2]],
  FN <- con_RF$table[[3]],
  TN <- con_RF$table[[4]]
) %>%  # # MCC <- ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP+FN) * (TN+FP) * (TN+FN))
  pmap_dbl(., ~ ((..1 * ..4) - (..2 * ..3))/sqrt((..1 + ..2) * (..1 + ..3) * (..4 + ..2) * (..4 + ..3)))

RF_pred_outcome <- cbind(as.numeric(RF_Predictions), as.numeric(Y_test)) %>% 
  data.frame() %>% 
  setNames(c("predictions", "outcome"))

RF_fg <- RF_pred_outcome %>% 
  filter(outcome == 1) %>%
  drop_na() %>% 
  pull(predictions)

RF_bg <- RF_pred_outcome %>% 
  filter(outcome == 0) %>%
  drop_na() %>% 
  pull(predictions)

RF_pr <- PRROC::pr.curve(
  scores.class0 = RF_fg, scores.class1 = RF_bg, curve = TRUE)


RF_pr_curve <- RF_pr$curve %>% 
  data.frame() %>% 
  mutate(X4 = "Random forest")

RF_prPlot <-
  ggplot(RF_pr_curve,aes(x = X1, y = X2, color = X4)) +
  geom_line() +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = c("Random forest")) +
  labs(title= paste0("Precision-Recall Curves,AUC=",RF_pr$auc.integral), 
       y = "Precision",
       x = "Recall") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0, 1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )

library(pROC)
RF_roc <- roc(as.numeric(Y_test),as.numeric(RF_Predictions))

RF_rocPlot <- 
  ggroc(RF_roc,legacy.axes = TRUE, linetype = 1) +
  geom_abline(show.legend = TRUE, alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = "Random forest") +
  labs(title=paste0("ROC Curves,AUC=",RF_roc$auc), 
       y = "Sensitivity",
       x = "1-Specificity") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "left",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )

###########################################################################

(RF_rocPlot + RF_prPlot) + plot_annotation(title = "Algorithm Comparison - Model Performance", tag_levels = "A", theme = theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")))

ggsave("classTree_rocPRPlot.pdf")
write.table(Random_Forest_Model$importance,file="importance_RF.txt",quote=T,sep="\t",col.names=T)

write.table(Random_Forest_Model$importanceSD,file="importanceSD_RF.txt",quote=T,sep="\t",col.names=T)

#####################################################################
#####################################################################
#####################################################################
#####################################################################
##############################adaboost
library(JOUSBoost)
adaBoost_Y_train <- ifelse(Y_train == 0, -1, 1) # adaBoost expects a 0 and 1 prediction
adaBoost_Model <- adaboost(X_train, adaBoost_Y_train, tree_depth = 8, n_rounds = 1000, verbose = T, control = NULL)
adaBoost_Predictions <- predict(adaBoost_Model, X_test)
adaBoost_Predictions <- ifelse(adaBoost_Predictions == -1, 0, 1) # convert back to 0 and 1 predictions


con_adaBoost <- confusionMatrix(factor(adaBoost_Predictions), factor(Y_test))

MCC_adaBoost <- list(
  TP <- con_adaBoost$table[[1]],
  FP <- con_adaBoost$table[[2]],
  FN <- con_adaBoost$table[[3]],
  TN <- con_adaBoost$table[[4]]
) %>%  # # MCC <- ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP+FN) * (TN+FP) * (TN+FN))
  pmap_dbl(., ~ ((..1 * ..4) - (..2 * ..3))/sqrt((..1 + ..2) * (..1 + ..3) * (..4 + ..2) * (..4 + ..3)))

adaBoost_pred_outcome <- cbind(as.numeric(adaBoost_Predictions), as.numeric(Y_test)) %>% 
  data.frame() %>% 
  setNames(c("predictions", "outcome"))

adaBoost_fg <- adaBoost_pred_outcome %>% 
  filter(outcome == 1) %>%
  pull(predictions)

adaBoost_bg <- adaBoost_pred_outcome %>% 
  filter(outcome == 0) %>%
  pull(predictions)

adaBoost_pr <- PRROC::pr.curve(
  scores.class0 = adaBoost_fg, scores.class1 = adaBoost_bg, curve = TRUE)

adaBoost_pr_curve <- adaBoost_pr$curve %>% 
  data.frame() %>% 
  mutate(X4 = "adaBoost")


adaBoost_prPlot <-
  ggplot(adaBoost_pr_curve,aes(x = X1, y = X2, color = X4)) +
  geom_line() +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = c("adaBoost")) +
  labs(title= paste0("Precision-Recall Curves",adaBoost_pr$auc.integral), 
       y = "Precision",
       x = "Recall") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0, 1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )

library(pROC)
adaBoost_roc <- roc(as.numeric(Y_test),as.numeric(adaBoost_Predictions))



adaBoost_rocPlot <- 
  ggroc(adaBoost_roc,legacy.axes = TRUE, linetype = 1) +
  geom_abline(show.legend = TRUE, alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = "Random forest") +
  labs(title= paste0("ROC Curves",adaBoost_roc$auc), 
       y = "Sensitivity",
       x = "1-Specificity") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "left",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )

###########################################################################

(adaBoost_rocPlot + adaBoost_prPlot) + plot_annotation(title = "Algorithm Comparison - Model Performance", tag_levels = "A", theme = theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")))

ggsave("adaBoost_rocPRPlot.pdf")


############################################################################## ################################GBM########################### ########################################################################## library(caret)

fitControl<- trainControl(method = "repeatedcv",number = 10, repeats = 10)

gbm.grid <- expand.grid(interaction.depth = 6:8,
                        n.trees = c(100,200,300,400,500),
                        shrinkage =c(0.001,0.01,0.1),
                        n.minobsinnode=c(10,20,30))
set.seed(123)
Y_train=factor(Y_train)

fit.gbm<-train(X_train,Y_train,method='gbm',trControl=fitControl,tuneGrid = gbm.grid)
fit.gbm

results_df=fit.gbm$results
trellis.par.set(caretTheme())

plot(fit.gbm)

fit.gbm$bestTune

set.seed(123) 
fit.gbm.final<-train(X_train,Y_train,method='gbm',trControl=fitControl,tuneGrid = fit.gbm$bestTune)
result.final=fit.gbm.final$results
GBM_predictions<-predict.train(object=fit.gbm.final,X_test)

library(caret)
con_GBM<- confusionMatrix(factor(GBM_predictions), factor(Y_test))
MCC_GBM <- list(
  TP <- con_GBM$table[[1]],
  FP <- con_GBM$table[[2]],
  FN <- con_GBM$table[[3]],
  TN <- con_GBM$table[[4]]
)%>% ##MCC<-((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)* (TN+FN))
  pmap_dbl(., ~ ((..1 * ..4) - (..2 * ..3))/sqrt((..1 + ..2) * (..1 + ..3) *
                                                   (..4 + ..2) * (..4 + ..3)))

GBM_pred_outcome <- cbind(as.numeric(GBM_predictions), as.numeric(Y_test)) %>%
  data.frame() %>%
  setNames(c("predictions", "outcome"))
################### PR curve
GBM_fg <- GBM_pred_outcome %>%
  filter(outcome == 1) %>%
  drop_na() %>%
  pull(predictions)

GBM_bg <- GBM_pred_outcome %>%
  filter(outcome == 0) %>%
  drop_na() %>%
  pull(predictions)

GBM_pr <- PRROC::pr.curve(
  scores.class0 = GBM_fg, scores.class1 = GBM_bg, curve = TRUE)


GBM_pr_curve <- GBM_pr$curve %>%
  data.frame() %>%
  mutate(X4 = "GBM Regression")
library(pROC)
GBM_roc <- roc(as.numeric(Y_test),as.numeric(GBM_predictions))

GBM_pr_Plot <-
  ggplot(GBM_pr_curve,aes(x = X1, y = X2, color = X4)) +
  geom_line() +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = c("GBM Regression")) +
  labs(title= paste0("Precision-Recall Curves,
AUC=",round(GBM_pr$auc.integral,3)),
       y = "Precision",
       x = "Recall") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0,
                                                                              1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

GBM_pr_Plot
######################### ROC Curves ##################################
GBM_roc_Plot <-
  ggroc(GBM_roc,legacy.axes = TRUE, linetype = 1) +
  geom_abline(show.legend = TRUE, alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Model",
                        
                        labels = c("GBM Regression")) +
  labs(title=  paste0("ROC Curves, AUC=",round(GBM_roc$auc,3)),
       y = "Sensitivity",
       x = "1-Specificity") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "left",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )
GBM_roc_Plot
###########################################################################
(GBM_roc_Plot + GBM_pr_Plot) + plot_annotation(title = "Algorithm Comparison -
Model Performance", tag_levels = "A", theme = theme(plot.title =
                                                      element_text(hjust = 0.5), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")))
ggsave("./GBM_rocPRPlot.pdf",width = 9,height = 7)


######################### #########################
######################### #########################
######################### XGBoost #########################
dtrain <- xgb.DMatrix(data = X_train, label = Y_train)
cv_folds <- KFold(Y_train, nfolds = 10, stratified = TRUE, seed = 0)

xgb_cv_bayes <- function(eta, max.depth, min_child_weight, subsample) {
  cv <- xgb.cv(
    params = list(
      booster = "gbtree",
      eta = eta,
      max_depth = max.depth,
      min_child_weight = min_child_weight,
      subsample = subsample,
      colsample_bytree = 0.6,
      lambda = 1,
      alpha = 0,
      objective = "binary:logistic",
      eval_metric = "auc"),
    data = dtrain, 
    nround = 10,
    folds = cv_folds,
    prediction = TRUE,
    showsd = TRUE,
    early.stop.round = 5,
    maximize = TRUE,
    verbose = 0
  ) 
  list(Score = cv$evaluation_log[, max(test_auc_mean)], Pred = cv$pred
  )}


OPT_Res <- BayesianOptimization(
  xgb_cv_bayes,
  bounds = list(
    eta = c(0.01L, 0.05L, 0.1L, 0.3L),
    max.depth = c(6L, 8L, 12L),
    min_child_weight = c(1L, 10L),
    subsample = c(0.5, 0.8, 1)
  ),
  init_grid_dt = NULL,
  init_points = 10,

  n_iter = 10,
  acq = "ucb",
  kappa = 2.576,
  eps = 0.0,
  verbose = TRUE
)

params <- list(
  "eta" = unname(OPT_Res$Best_Par["eta"]),
  "max_depth" = unname(OPT_Res$Best_Par["max.depth"]),
  "colsample_bytree" = 1,
  "min_child_weight" = unname(OPT_Res$Best_Par["min_child_weight"]),
  "subsample"= unname(OPT_Res$Best_Par["subsample"]),
  "objective"="binary:logistic",
  "gamma" = 1,
  "lambda" = 1,
  "alpha" = 0,
  "max_delta_step" = 0,
  "colsample_bylevel" = 1,
  "eval_metric"= "auc",
  "set.seed" = 176
)

params <- list(
  "eta" = 0.05,
  "max_depth" = 6,
  "colsample_bytree" = 1,
  "min_child_weight" = 1,
  "subsample"= 0.73,
  "objective"="binary:logistic",
  "gamma" = 1,
  "lambda" = 1,
  "alpha" = 0,
  "max_delta_step" = 0,
  "colsample_bylevel" = 1,
  "eval_metric"= "auc",
  "set.seed" = 176
)
watchlist <- list("train" = dtrain)
nround = 20
xgb.model <- xgb.train(params, dtrain, nround, watchlist)
dtest <- xgb.DMatrix(data = X_test)
predictions <- predict(object = xgb.model, newdata = dtest, type = 'prob') 

XGBPredictions <- ifelse(predictions > 0.5, 1, 0)
results <- cbind(test, predictions) %>%
  mutate(
    pred_status = case_when
      predictions <= 0.50 ~ 0
    ),
    correct = case_when(
      outcome == pred_status ~ "Correct",
      TRUE ~ "Incorrect"
    ),
    outcome_text = case_when(
      outcome == 1 ~ "Perished",
      outcome == 0 ~ "Survived"
    ) ) %>%
  mutate_at(
    vars(outcome, pred_status, correct, outcome_text), funs(factor)
  ) %>%
  add_count(outcome, pred_status)

con_XGB <- confusionMatrix(factor(XGBPredictions), factor(Y_test))
MCC_XGB <- list(
  TP <- con_XGB$table[[1]],
  FP <- con_XGB$table[[2]],
  FN <- con_XGB$table[[3]],
  TN <- con_XGB$table[[4]]
)%>% ##MCC<-((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)* (TN+FN))
  pmap_dbl(., ~ ((..1 * ..4) - (..2 * ..3))/sqrt((..1 + ..2) * (..1 + ..3) *
                                                   (..4 + ..2) * (..4 + ..3)))

XGB_pred_outcome <- cbind(as.numeric(predictions), as.numeric(Y_test)) %>%
  data.frame() %>%
  setNames(c("predictions", "outcome"))
XGB_fg <- XGB_pred_outcome %>%
  filter(outcome == 1) %>%
  pull(predictions)
XGB_bg <- XGB_pred_outcome %>%
  filter(outcome == 0) %>%
  pull(predictions)
XGB_pr <- PRROC::pr.curve(
  scores.class0 = XGB_fg, scores.class1 = XGB_bg, curve = TRUE)

XGB_pr_curve <- XGB_pr$curve %>%
  data.frame() %>%
  mutate(X4 = "XGBoost")

library(pROC)

XGB_roc <- roc( as.numeric(Y_test),as.numeric(predictions))
pROC::plot.roc(as.numeric(Y_test),as.numeric(predictions),print.thres=T)
XGB_prPlot <-
  ggplot(XGB_pr_curve,aes(x = X1, y = X2, color = X4)) +
  geom_line() +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = 'XGBoost') +
  labs(title= paste0("Precision-Recall Curves,
AUC=",round(XGB_pr$auc.integral,3)),
       y = "Precision",
       x = "Recall") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0,
                                                                              1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

XGB_rocPlot<-
  ggroc(XGB_roc,legacy.axes = TRUE, linetype = 1) +
  geom_abline(show.legend = TRUE, alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = 'XGboost') +
  labs(title=  paste0("ROC Curves, AUC=",round(XGB_roc$auc,3)),
       y = "Sensitivity",
       x = "1-Specificity") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "left",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )
###########################################################################
(XGB_rocPlot + XGB_prPlot) + plot_annotation(title = "Algorithm Comparison -
Model Performance", tag_levels = "A", theme = theme(plot.title =
                                                      element_text(hjust = 0.5), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")))
ggsave("XGBPRPlot.pdf")


################################
#################################
#################light-GBM
#################################
#################################
#################################
library(lightgbm)
lgbm_train <- lgb.Dataset(data = X_train, label = Y_train)
lgbm_test <- lgb.Dataset(data = X_test, label = Y_test)

grid_search <- expand.grid(
  Depth = 8,
  L1 = 0:5,
  L2 = 0:5,
  MinHessianLeaf = 0:2,
  featureFraction = c(0.8L, 1L)
)
grid_search
model <- list()
perf <- numeric(nrow(grid_search)) 
valids <- list(test = lgbm_test) 
library(data.table)


for (i in 1:nrow(grid_search)) {
  model[[i]] <- lgb.train(
    list(
      objective = "binary",
      metric = "l2",
      lambda_l1 = grid_search[i, "L1"],
      lambda_l2 = grid_search[i, "L2"],
      max_depth = grid_search[i, "Depth"],
      min_sum_hessian_in_leaf = grid_search[i, "MinHessianLeaf"],
      feature_fraction = grid_search[i, "featureFraction"]
    ),
    lgbm_train,
    2,
    valids,
    min_data = 1,
    learning_rate = 1,
    early_stopping_rounds = 5)
  perf[i] <- min(rbindlist(model[[i]]$record_evals$test$l2))
}
cat("Model ", which.min(perf), " is lowest loss: ", min(perf), sep = "")
print(grid_search[which.min(perf), ])



params_lightGBM <- list(
  objective = "binary",
  metric = "auc",
  #默认0.0001 min_sum_hessian_in_leaf = 0,
  
  feature_fraction = 0.8, #默认0
  lambda_l1 = 1,
  #默认0
  lambda_l2 = 0,
  is_unbalance = FALSE
)
lightGBM.model <- lgb.train( data = lgbm_train,
                             params = params_lightGBM, #默认0.1
                             learning_rate = 0.1,
                             nrounds = 60 )
lightGBM_Predictions <- predict(object = lightGBM.model, data = X_test, rawscore = FALSE)
lightGBM_Predictions_bi <- ifelse(lightGBM_Predictions > 0.5, 1, 0)
con_LGBM <- confusionMatrix(factor(lightGBM_Predictions_bi), factor(Y_test))
############## Light GBM 评价 ##################### 
############## 
############## 
MCC_LGBM <- list(
  TP <- con_LGBM$table[[1]],
  FP <- con_LGBM$table[[2]],
  FN <- con_LGBM$table[[3]],
  TN <- con_LGBM$table[[4]]
)%>% ##MCC<-((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)* (TN+FN))
  pmap_dbl(., ~ ((..1 * ..4) - (..2 * ..3))/sqrt((..1 + ..2) * (..1 + ..3) *
                                                   (..4 + ..2) * (..4 + ..3)))

LGBM_pred_outcome <- cbind(as.numeric(lightGBM_Predictions), as.numeric(Y_test))%>%
  data.frame() %>%
  setNames(c("predictions", "outcome"))
LGBM_fg <- LGBM_pred_outcome %>%
  filter(outcome == 1) %>%
  pull(predictions)
LGBM_bg <- LGBM_pred_outcome %>%
  filter(outcome == 0) %>%
  pull(predictions)
LGBM_pr <- PRROC::pr.curve(
  scores.class0 = LGBM_fg, scores.class1 = LGBM_bg, curve = TRUE)


LGBM_pr_curve <- LGBM_pr$curve %>%
  data.frame() %>%
  
  mutate(X4 = "Light GBM")
LGBM_prPlot <-
  ggplot(LGBM_pr_curve,aes(x = X1, y = X2, color = X4)) +
  geom_line() +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = c("Light GBM")) +
  labs(title= paste0("Precision-Recall Curves,
AUC=",round(LGBM_pr$auc.integral,3)),
       y = "Precision",
       x = "Recall") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0,
                                                                              1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

library(pROC)
LGBM_roc <- roc(as.numeric(Y_test),as.numeric(lightGBM_Predictions))
LGBM_rocPlot <-
  ggroc(LGBM_roc,legacy.axes = TRUE, linetype = 1) +
  geom_abline(show.legend = TRUE, alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = "Light GBM") +
  labs(title=  paste0("ROC Curves, AUC=",round(LGBM_roc$auc,3)),
       y = "Sensitivity",
       x = "1-Specificity") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  theme_bluewhite() +
  theme(
    aspect.ratio = 1,
    legend.position = "left",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

###########################################################################
(LGBM_rocPlot + LGBM_prPlot) + plot_annotation(title = "Algorithm Comparison -
Model Performance", tag_levels = "A", theme = theme(plot.title =
                                                      element_text(hjust = 0.5), 
                                                    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")))
ggsave("LGBM_rocPRPlot.pdf")


##########################################################################
##########################################################################
#############################SHAP Model###################################
##########################################################################
library(data.table)
library(ggplot2)
library(xgboost)
library(SHAPforxgboost)


params = list(objective = 'reg:squarederror',
              #num_boost_round  = 500,
              learning_rate = 0.05,
              max_depth = 8,
              #gamma = 0.009,
              subsample = 0.8,
              colsample_bytree = 0.8
)

dim(X_train)
dtrain = xgb.DMatrix(X_train, label = Y_train)
model = xgb.train(params,
                  dtrain,
                  verbose = 50,
                  nrounds = 500)
imp = xgb.importance(model = model)
xgb.plot.importance(imp, 
                    rel_to_first = T, 
                    xlab = 'Relative importance',
                    top_n = 10)

#Now we use `xgb.plot.shap` to show feature importance magnitude.
shap_values = shap.values(xgb_model = model, X_train = X_train)

shap_long = shap.prep(xgb_model = model, X_train = X_train)
# is the same as: using given shap_contrib
shap_long = shap.prep(shap_contrib = shap_values$shap_score, X_train = X_train)

shap.plot.summary.wrap1(model, X = X_train, top_n = 15) 
shap.plot.summary.wrap1(model, X = X_train) 
shap.plot.summary.wrap1(model, X = X_train, top_n = 4) 
#ow in order to show detailed relation between features and target output, we create some SHAP Dependence Plots. In R package `xgboost`, it has a bulit-in function to plot SHAP dependence for most important features. 

xgb.plot.shap(X_train, model = model, features = 'COX5A')

########`xgb.plot.shap` could also shows multiple features:
########
xgb.plot.shap(X_train, model = model, top_n = 4, n_col = 2)

####Interaction effects
#The SHAP interaction values take time since it calculates all the combinations.

# prepare the data using either: 
# (this step is slow since it calculates all the combinations of features.)
shap_int <- shap.prep.interaction(xgb_mod = model, X_train = X_train)
# or:
shap_int <- predict(model, X_train, predinteraction = TRUE) # (the same)
# **SHAP interaction effect plot **
# if `data_int` is supplied, the same function will plot the interaction effect:
g3 <- shap.plot.dependence(data_long = shap_long,
                           data_int = shap_int,
                           x= "COX5A", y = "COX4I2", 
                           color_feature = "COX4I2")
g4 <- shap.plot.dependence(data_long = shap_long,
                           data_int = shap_int,
                           x= "COX5A", y = "COX5B", 
                           color_feature = "COX5B")
gridExtra::grid.arrange(g3, g4, ncol=2)

#Plot SHAP value against feature value, without color_feature but has marginal distribution:

fig_list <- lapply(names(shap_values$mean_shap_score)[1:4], 
                   shap.plot.dependence, data_long = shap_long)
gridExtra::grid.arrange(grobs = fig_list, ncol = 4)



################SHAP force plot

#The SHAP force plot basically stacks these SHAP values for each observation, and show how the final output was obtained as a sum of each predictor’s attributions.

# choose to show top 4 features by setting `top_n = 4`, 
# set 6 clustering groups of observations.  
plot_data <- shap.prep.stack.data(shap_contrib = shap_values$shap_score, top_n = 4, n_groups = 4)

#plot_data <- shap.prep.stack.data(shap_contrib = shap_values$shap_score, top_n = 15)
# you may choose to zoom in at a location, and set y-axis limit using `y_parent_limit`  
shap.plot.force_plot(plot_data, zoom_in_location = 20, y_parent_limit = c(-0.8,0.8))


# plot the 6 clusters
shap.plot.force_plot_bygroup(plot_data)

#########################

#From this plot, we could know how SHAP values distributed on every single feature. 

#Moreover, `SHAPforxgboost` provides us `ggplot2` backend to make SHAP Dependence plots:
shap.plot.dependence(data_long = shap_long, 
                     x = 'PCK2', 
                     y = 'GOBP_REGULATION_OF_NAD_P_H_OXIDASE_ACTIVITY', 
                     color_feature = 'GOBP_REGULATION_OF_NAD_P_H_OXIDASE_ACTIVITY') + 
  labs('')

shap.plot.dependence(data_long = shap_long, 
                     x = 'FOXK2', 
                     y = 'GOBP_REGULATION_OF_NAD_P_H_OXIDASE_ACTIVITY', 
                     color_feature = 'GOBP_REGULATION_OF_NAD_P_H_OXIDASE_ACTIVITY') + 
  labs('')

shap.plot.dependence(data_long = shap_long, 
                     x = 'NQO1', 
                     y = 'GOBP_REGULATION_OF_NAD_P_H_OXIDASE_ACTIVITY', 
                     color_feature = 'GOBP_REGULATION_OF_NAD_P_H_OXIDASE_ACTIVITY') + 
  labs('')

shap.plot.dependence(data_long = shap_long, 
                     x = 'LDHA', 
                     y = 'GOBP_REGULATION_OF_NAD_P_H_OXIDASE_ACTIVITY', 
                     color_feature = 'GOBP_REGULATION_OF_NAD_P_H_OXIDASE_ACTIVITY') + 
  labs('')
######################################################################################
######################################################################################
shap.plot.dependence(data_long = shap_long, 
                     x = 'PCK2', 
                     y = 'GOMF_NADP_BINDING', 
                     color_feature = 'GOMF_NADP_BINDING') + 
  labs('')

shap.plot.dependence(data_long = shap_long, 
                     x = 'FOXK2', 
                     y = 'GOMF_NADP_BINDING', 
                     color_feature = 'GOMF_NADP_BINDING') + 
  labs('')

shap.plot.dependence(data_long = shap_long, 
                     x = 'NQO1', 
                     y = 'GOMF_NADP_BINDING', 
                     color_feature = 'GOMF_NADP_BINDING') + 
  labs('')

shap.plot.dependence(data_long = shap_long, 
                     x = 'LDHA', 
                     y = 'GOMF_NADP_BINDING', 
                     color_feature = 'GOMF_NADP_BINDING') + 
  labs('')

########################################
deg=read.csv('ELECTRON_TRANSPORT.txt',sep = '\t',header=T)
exp=read.csv('Exp_ALL.txt',sep = '\t',header=T)
#exp = read.table("Post_Combat.genesymbol.txt",header=T,sep="\t")
deg.exp = exp[match(deg$ID,exp$ID),]
write.table(deg.exp,"ELECTRON_TRANSPORT_EXP.txt.txt",col.names = T,row.names = F,quote=F,sep="\t")
#####################################################
##############################################

##########options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
#install.packages("ClassDiscovery")
library(ClassDiscovery) 
library(pheatmap) 
library(gplots) 
library(Seurat) 
library(gplots)
library("corrplot")
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

blue = "#307eb9"
purple = "#974e9e"
red = "#e41e25"
orange = "#f57f21"
yellow = "#f4ed35"
lblue = "#74d7ff"
red2 = "#ff1751"
green = "#4cb049"
grey="#b9afaf"
black="#1d1515"

mat <- read.csv("Mito_PathwayScore_BH.csv",row.names = 1,check.names = F,stringsAsFactors = F)
mat[1:3,1:3]

#data <- ScaleData(mat, verbose = FALSE)
#mat<-data
#data <- NormalizeData(mat,verbose = FALSE)	
#mat <- ScaleData(data, verbose = FALSE)	

annCol <- read.csv("GroupAnno.csv",row.names = 1,check.names = F,stringsAsFactors = F)
head(annCol)

annCol[is.na(annCol) | annCol == ""] <- "N/A"


annColors <- list()
annColors[["group"]] <- c("Control"="black","PAH"="Red")


annColors

library(Seurat)

mat3<-ScaleData(mat,features = NULL, vars.to.regress = NULL,
                latent.data = NULL,split.by = NULL,  model.use = "linear", 
                do.scale = TRUE,  do.center = TRUE,
                scale.max = 2,  block.size = 1000,  verbose = TRUE)

pheatmap(mat,
         scale = "row",
         color = COL2('BrBG'), # redblue(20) ,####redblue(20)   greenred(20)  COL2('BrBG')
         #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         annotation_col = annCol,
         cluster_rows=T,cluster_cols = F,
         annotation_colors = annColors,
         show_rownames = T, show_colnames = F)

################################B-H#########################################
library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)
library(limma)
rt<-read.table("Mito_PathwayScore_BH.txt",header=T,sep="\t",row.names=1)
dim(rt)
#differential  group = c(rep("a", 2095), rep("b", 2039))
class<-c(rep("a", 58), rep("b", 164))
#class<-c(rep("control",2095),rep("Is",2039))
design<-model.matrix(~factor(class))
colnames(design)<-c("control","Is")
fit<-lmFit(rt,design)
fit2<-eBayes(fit)

x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)

write.csv(x, "gsva_limma_Metab.csv", quote = F)

df <- read.csv("gsva_limma_Metab.csv")
head(df)

cutoff <- 1
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))

sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)


ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('palegreen3', 'snow3', 'dodgerblue4'), guide = FALSE) + 
  
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, 
             size = 0.3) + 
  
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID), color = group),
            size = 3, 
            hjust = "inward" ) +  
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1, label=ID, color = group),
            size = 3, hjust = "outward") +  
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score")+
  theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 0.6)) + 
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) 

