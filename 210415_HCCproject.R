
#dir.create("E:/Jeonju2021/Data/Research")
setwd("E:/Jeonju2021/Data/Research/Data")
library(GEOquery)
library(simpleaffy)
library(tidyverse)
library(RColorBrewer)
library(genefilter)
library(pheatmap)
library(ggrepel)
library(dplyr)
library(limma)
library(affy)
#RegParallel

BiocManager::install("monocle3")
BiocManager::install("limma")


n

#gse <- getGEO("GSE14520", GSEMatrix = TRUE)
#save(gse, file = "GSE14520.Rdata")
load("GSE14520.Rdata")
# we analyzing the first GEO and the data including two different array
gse <- gse[[1]]
length(gse)
show(gse)
pData(gse) ## print the sample information
fData(gse) ## print the gene annotation
exprs(gse) ## print the expression data
summary(exprs(gse))
exprs(gse) <- log2(exprs(gse))
boxplot(exprs(gse),outline=FALSE)
library(dplyr)
sampleInfo <- pData(gse)
str(sampleInfo)
sampleInfo <- dplyr::select(sampleInfo, characteristics_ch1,characteristics_ch1.1)

## Optionally, rename to more convenient column names
#sampleInfo <- rename(sampleInfo,group = "characteristics_ch1.1" , patient = "characteristics_ch1" )
#sampleInfo <- select(sampleInfo, group)
sampleInfo$characteristics_ch1
sampleInfo$group[sampleInfo$characteristics_ch1 == "tissue: Liver Non-Tumor Tissue"] <- "Normal"
sampleInfo$group[sampleInfo$characteristics_ch1 == "Tissue: Liver Non-Tumor Tissue"] <- "Normal"
sampleInfo$group[sampleInfo$characteristics_ch1 == "tissue: Liver Tumor Tissue"] <- "Tumor"
sampleInfo$group[sampleInfo$characteristics_ch1 == "Tissue: Liver Tumor Tissue"] <- "Tumor"

table(sampleInfo$group)

library(ggrepel)
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(exprs(gse)))

## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=group,label=paste(""))) + geom_point() + geom_text_repel()

cutoff <- median(exprs(gse), na.rm = T)
## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(gse) > cutoff
## Identify genes expressed in more than 2 samples
keep <- rowSums(is_expressed) > 2
## check how many genes are removed / retained.
table(keep)
## subset to just those expressed genes
gse <- gse[keep,]
#DEGs
design <- model.matrix(~0+sampleInfo$group)
colnames(design) <- c("Tumour", "Normal")

library(limma)

## subset to just those expressed genes
fit <- lmFit(exprs(gse), design)
head(fit$coefficients)
contrasts <- makeContrasts(Tumour - Normal, levels=design)
## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
anno <- fData(gse)
fit2$genes <- anno
res_all1 <- topTable(fit2, adjust="fdr",sort.by = "P", n = Inf)
#res_all1 <- topTable(fit2, adjust="BH",sort.by = "P", n = Inf)


head(res_all1)

#res_all1 <- topTable(fit2,sort.by = "P", n = Inf)

#head(res_all1)

#length((res_all1$Gene.Symbol))


str(res_all1)


library(gtools)
res_all1$FC <- logratio2foldchange(res_all1$logFC, base = 2)
res_all1 <- subset(res_all1, res_all1$adj.P.Val < 0.05)
mRNAS_DEG_HCC14520 <- subset(res_all1,abs(res_all1$FC) >=1.0)

mRNAS_DEG_HCC14520 =  mRNAS_DEG_HCC14520[order(mRNAS_DEG_HCC14520$Gene.Symbol,mRNAS_DEG_HCC14520$adj.P.Val),]
mRNAS_DEG_HCC14520=mRNAS_DEG_HCC14520[!duplicated(mRNAS_DEG_HCC14520$Gene.Symbol),]
length((mRNAS_DEG_HCC14520$Gene.Symbol))

up=subset(mRNAS_DEG_HCC14520,FC > 0)
down=subset(mRNAS_DEG_HCC14520,FC < 0)
length((mRNAS_DEG_HCC14520$Gene.Symbol))
length((up$Gene.Symbol))
length((down$Gene.Symbol))


length((res_all1$Gene.Symbol))


table(res_all1$"P.Value" <0.05)
table(res_all1$"adj.P.Val" <0.05)



library(dplyr)
library(readr)
features <- fData(gse)
#View(features)
### Look at the features data frame and decide the names of the columns you want to keep
features <- dplyr::select(features,ID,'Gene Symbol',ENTREZ_GENE_ID)
full_output <- cbind(features,exprs(gse))
#full_output<- full_output[full_output$ID %in% mRNAS_DEG_HCC14520$ID,]
#full_output <- subset(full_output, select = -c(25:100))
head(full_output)
str(mRNAS_DEG_HCC14520)

dir.create("D:/Jeonju2021/Data/Research/Plosone/Revised/210405")
setwd("E:/Jeonju2021/Data/Research/Plosone/Revised/210405")
save(full_output, file = "mRNA_HCC14520.Rdata")
save(mRNAS_DEG_HCC14520, file = "mRNAS_DEG_HCC14520.Rdata")
write.csv(full_output,"mRNA_HCC14520.csv")





RB <- c("CDK1","RFC4", "KIF14", "PRIM1","NEK2", "TOP2A", "RRM2", "CCNB1", "DBF4","KIF14","CG018","CENPA","H2AFX") 
RB <- data.frame(RB)
expr_df = full_output
expr_df <- expr_df[expr_df$ "Gene Symbol" %in% RB$RB,]
expr_df$"Gene Symbol"
write.csv(expr_df,"mRNA_HCC14520_top10.csv")







#WGCNA
library(dynamicTreeCut)
library(stats)
library(fastcluster)
#biocLite("GO.db")
library(WGCNA)
library(lattice)
library(latticeExtra)

#Read file
rt=full_output
rt = subset(rt, select = -c(1,3) )
names(rt)[1] <- "genename"
#rt$genename
s<-res_all1
datExpr0.group =rt
rt<- datExpr0.group[rownames(datExpr0.group) %in% rownames(s),]
#rt=cbind(features,exprs(gse))#new gse saved to c2
rt=as.matrix(rt)
head(rt)
head(s)
ncol(rt)
nrow(rt)
rt=rt[,1:446]
rownames(rt)=rt[,1]
colnames(rt)
exp=rt[,2:ncol(rt)]
dimname=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow = nrow(exp),dimnames = dimname)
rt <- as.data.frame(t(rt)) #
##
sampleTree = hclust(dist(rt), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
#
powers = c(c(1:10), seq(from = 9, to=20, by=2))
# 
sft = pickSoftThreshold(rt, powerVector = powers, verbose = 5)
#
sizeGrWindow(9, 5)
par(mfrow = c(1,2)); 
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
#
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#
softPower <- sft$powerEstimate 
adjacency = adjacency(rt, power = softPower);
k <- softConnectivity(datE=rt,power=softPower) 
sizeGrWindow(10, 5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")

##

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM       
hierTOM = hclust(as.dist(dissTOM),method="average");

# 
geneTree = hclust(as.dist(dissTOM), method = "average");

saveRDS(object = geneTree,file = "geneTree.RDS",compress = FALSE)
geneTree = readRDS(file = "geneTree.RDS")

memory.limit(1000000)
# Plot the resulting clustering tree (dendrogram)
#windows()
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.01);
# 
minModuleSize = 30;
# 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# 
nSelect = 400 
# For reproducibility, we set the random seed 
set.seed(10); 
select = sample(16074, size = nSelect)
selectTOM = dissTOM[select, select]
dynamicColors=labels2colors(dynamicMods)
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster. 
selectTree = hclust(as.dist(selectTOM), method = "average") 
selectColors = dynamicColors[select]; 
# Open a graphical window 
sizeGrWindow(9,9) 
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot 
plotDiss = selectTOM^softPower; 
diag(plotDiss) = NA; 
TOMplot(plotDiss, 
        selectTree, 
        selectColors, 
        main = "Network heatmap plot, selected genes") 

dynamicColors=labels2colors(dynamicMods)
MEList = moduleEigengenes(rt, colors = dynamicColors)
MEs = MEList$eigengenes
#
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
plotEigengeneNetworks(MEs, 
                      "Eigengene adjacency heatmap", 
                      marHeatmap = c(3,4,2,2), 
                      plotDendrograms = FALSE, 
                      xLabelsAngle = 90) 
#clustering
plot(METree, 
     main = "Clustering of module eigengenes",
     xlab = "", 
     sub = "")
#
MEDissThres = 0.2
abline(h=MEDissThres, col = "red")
merge_modules = mergeCloseModules(rt, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge_modules$colors;
mergedMEs = merge_modules$newMEs;
## to plot dynamic cut tree
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dataTraits <- read.csv("dataTraits3.csv",header=TRUE,row.names = 1)
dataTraits=as.matrix(dataTraits)
dataTrait<-dataTraits
dataTrait<-dataTrait[1:445,]
ncol(dataTrait)
nrow(dataTrait)
ncol(mergedMEs)
nrow(mergedMEs)
cor_ADR <- signif(WGCNA::cor(dataTrait,mergedMEs,use="p",method="pearson"),5)
p.values <- corPvalueStudent(cor_ADR,nSamples=445)
Freq_MS_max_cor <- which.max(abs(cor_ADR))
Freq_MS_max_p <- which.min(p.values)
GS1 <- as.numeric(WGCNA::cor(dataTrait,rt,use="p",method="pearson"))
GeneSignificance <- abs(GS1)
ModuleSignificance <- tapply(GeneSignificance,mergedColors,mean,na.rm=T)
Find_max_ModuleSign <- which.max(ModuleSignificance)
Find_max_ModuleSign 

#Target module
module = "darkmagenta";
probes = colnames(rt)
inModule = (mergedColors==module);
modProbes = probes[inModule]; 
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
#Network visulization via cytoscape
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = mergedColors[inModule]);
rtmod <- rt[inModule, inModule];
colorh1 = mergedColors[inModule]
ADJ1=abs(cor(rtmod,use="p"))^softPower 
Alldegrees1=intramodularConnectivity(ADJ1, colorh1) 
rtKME=signedKME(rt, mergedMEs, outputColumnName="MM.")
head(rtKME)
FilterGenes_spe = ((GeneSignificance > 0.2) & (abs(rtKME["MM.darkmagenta"])>0.8)) 
table(FilterGenes_spe)
trait_hubGenes_spe <- colnames(rt)[FilterGenes_spe] 
head(trait_hubGenes_spe)
trait_hubGenes_spe

str(full_output)

full_output<- full_output[full_output$ID %in% mRNAS_DEG_HCC14520$ID,]
full_output<- full_output[full_output$"Gene Symbol" %in% trait_hubGenes_spe,]

write.csv(full_output,"HCC_full_output.csv")


rt<- rt[rt$genename %in% CON.OR$X,]








trait_hubGenes_spe=t(trait_hubGenes_spe)
m=t(s[,1])
Gene=trait_hubGenes_spe%in%m
sum(trait_hubGenes_spe%in%m==TRUE)#(163)
head(trait_hubGenes_spe)
#### getting intra module connectivity 
CON.OR<-Alldegrees1
CON.OR <- CON.OR[order(CON.OR$kWithin, decreasing = T),]
CON.OR<-CON.OR[1:10,]
#write.csv(CON.OR,"CONOR20.csv")
CON.OR = read.csv("CONOR20.csv",header=TRUE)
CON.OR <- CON.OR[CON.OR$X %in% trait_hubGenes_spe,]
rt=read.csv("1_mRNA_HCC.csv",header = T,check.names = F)
rt = subset(rt, select = -c(1,3) )
names(rt)[1] <- "genename"
rt<- rt[rt$genename %in% CON.OR$X,]
rt$genename
#savefile
dir.create("D:/Jeonju2021/Data/Research/Plosone/Revised/WGCNA/2021FDR/210401")
setwd("D:/Jeonju2021/Data/Research/Plosone/Revised/WGCNA/2021FDR/210401")
write.csv(full_output,"1_mRNA_HCC.csv")
write.csv(DEgene,"2_DEG_HCC.csv")
write.csv(trait_hubGenes_spe,"3_GSMM_DEG_HCC.csv")
write.csv(Alldegrees1,"4_Kin_DEG_HCC.csv")
write.csv(CON.OR,"5_top20_DEG_HCC.csv")







### Network analysis
##Required libraries
library("readxl")
library(qgraph)
library(huge)
library("psychTools")
library(reshape2)
library(NetworkComparisonTest)
library(bootnet)
library(data.table)


rt=read.csv("1_mRNA_HCC.csv",header = T,check.names = F)
rt = subset(rt, select = -c(1,3) )
names(rt)[1] <- "genename"
head(rt)
rt2<- rt[rt$genename %in% CON.OR$X,]
write.csv(rt2,"5_top10_DEG_HCC.csv")

RB=read.csv("5_top10_DEG_HCC2.csv",header = T,check.names = F)


setwd("E:/Jeonju2021/Data/Research/Plosone/Revised/210405")
RB=read.csv("a1.csv",header = T,check.names = F)
RB=read.csv("a2.csv",header = T,check.names = F)
head(RB)
#RB<- RB[RB$Group ==1,]


RB <- RB[3:12]

#Numeric data
RB=sapply(RB,as.numeric)

#Nonparanormal Transformation
data1<- huge.npn(RB)

net1 <- estimateNetwork(data1,default = "EBICglasso",corMethod = "cor_auto",tuning = 0.5)

#Plot
L1 <- averageLayout(net1)


group.item <- list("Hub genes" = 1:3,
                   "Other" = 4:10)
color=c( "red","lightblue")
#color=c("lightblue", "lightgreen", "orange","yellow","lightyellow","green", "purple","red")
layout(t(1:1))
P2 <-plot(net1, groups = group.item,cut = 0.03,negDashed=F,layout = 'spring',color=color,legend=F)
P2
#centrality_auto(net1)

library(tidyverse)
library(stringi)
P1 <-centralityPlot(net1, standardized = TRUE)
Cent1 <- centralityTable(net1, standardized = TRUE)
table.data1 <- data.frame("Item" = Cent1[Cent1$measure == "Strength",3],
                          "Strength" = Cent1[Cent1$measure == "Strength",5])


library(ggplot2)
library(gridExtra)
P1 <-P1 +scale_colour_Publication()+ theme_Publication()
P2 <-P2 +scale_colour_Publication()+ theme_Publication()

grid.arrange(P1,P2,nrow=1)
gridExtra::grid.arrange(P1,P2, ncol = 2)







###SURVIVAL ANALYSIS
library("emmeans")
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library(miRBaseConverter)
library(enrichR)
library(VennDiagram)
library(SpidermiR)
library(DESeq2)
library(clusterProfiler)
library(wordcloud)
library(org.Hs.eg.db)
library(multcomp) 
library(Coxnet) 


setwd("D:/Jeonju2021/Data/Research/Cancer/IBD_HCC/210406")
tcga_data=load("gbmProteinExpression.rda")
tcga_data=data
#tcga_data=data
### check some information regarding to the clinical data 
clinical = as.data.frame(tcga_data@colData)
### we just check sum information 
str(clinical)
table(clinical$patient)
table(clinical$sample)
table(clinical$definition)
length(unique(clinical$patient))

#DEG

# we are only interested in the "Primary solid Tumor" cases for survival
clin_df = clinical[clinical$definition == "Primary solid Tumor",
                   c("patient",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "gender",
                     "tumor_stage")]



limma_pipeline = function(
  tcga_data,
  condition_variable,
  reference_group=NULL){
  
  design_factor = colData(tcga_data)[, condition_variable, drop=T]
  
  group = factor(design_factor)
  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}
  
  design = model.matrix(~ group)
  
  dge = DGEList(counts=assay(tcga_data),
                samples=colData(tcga_data),
                genes=as.data.frame(rowData(tcga_data)))
  
  # filtering
  keep = filterByExpr(dge,design)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)
  
  # Normalization (TMM followed by voom)
  dge = calcNormFactors(dge)
  v = voom(dge, design, plot=TRUE)
  
  # Fit model to data given design
  fit = lmFit(v, design)
  fit = eBayes(fit)
  
  # Show top genes
  topGenes = topTable(fit, coef=ncol(design), number=100, sort.by="p")
  
  return(
    list(
      voomObj=v, # normalized data
      fit=fit, # fitted model and statistics
      topGenes=topGenes # the 100 most differentially expressed genes
    )
  )
}


limma_res = limma_pipeline(
  tcga_data=tcga_data,
  condition_variable="definition",
  reference_group="Solid Tissue Normal"
)

# Save the data as a file, if you need it later, you can just load this file
# instead of having to run the whole pipeline again
saveRDS(object = limma_res,
        file = "limma_res.RDS",
        compress = FALSE)



#Survival Analysis
# extract clinical data
tcga_data=load("gbmProteinExpression.rda")
tcga_data=data
clinical = tcga_data@colData

dim(clinical)
# Transpose and make it into a matrix object
d_mat = as.matrix(t(limma_res$voomObj$E))

# As before, we want this to be a factor
d_resp = as.factor(limma_res$voomObj$targets$definition)



# create a new boolean variable that has TRUE for dead patients
# and FALSE for live patients
clin_df$deceased = clin_df$vital_status == "Dead"

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

# show first 10 samples
head(clin_df)

Surv(clin_df$overall_survival, clin_df$deceased)
Surv(clin_df$overall_survival, clin_df$deceased) ~ clin_df$gender
# fit a survival model
fit = survfit(Surv(overall_survival, deceased) ~ gender, data=clin_df)
print(fit)
str(clin_df)












RB <- c("CDK1","RFC4", "KIF14", "PRIM1","NEK2", "TOP2A", "RRM2", "CCNB1") 
RB <- data.frame(RB)
setwd("D:/Jeonju2021/Data/Research/Cancer/IBD_HCC/210406")
load("mRNAS_coexp1.Rdata")
expr_df = res_all_adj
expr_df <- expr_df[expr_df$external_gene_name %in% RB$RB,]
expr_df$external_gene_name

gene_id = expr_df[3, "ensgenes"]
gene_name1 = expr_df[3,"external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name1))
p5 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.715 \n p = 0.002",size =5)
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)

gene_id = expr_df[7, "ensgenes"]
gene_name = expr_df[7, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p6 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.769 \n p = 0.001",size =5)

gene_id = expr_df[1, "ensgenes"]
gene_name = expr_df[1, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p7 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 2.057 \n p < 0.0001",size =5)

gene_id = expr_df[4, "ensgenes"]
gene_name = expr_df[4, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p8 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.958 \n p = 0.0002",size =5)


gene_id = expr_df[6, "ensgenes"]
gene_name1 = expr_df[6,"external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name1))
p1 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.751 \n p = 0.0017",size =5)
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)



gene_id = expr_df[5, "ensgenes"]
gene_name = expr_df[5, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p2 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.812 \n p = 0.0009",size =5)

gene_id = expr_df[2, "ensgenes"]
gene_name = expr_df[2, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p3 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.627 \n p = 0.006",size =5)

gene_id = expr_df[8, "ensgenes"]
gene_name = expr_df[8, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p4 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.549 \n p = 0.013",size =5)




RB <- c("CDK1","RFC4", "KIF14", "PRIM1","NEK2", "TOP2A", "RRM2", "CCNB1") 
RB <- data.frame(RB)
setwd("D:/Jeonju2021/Data/Research/Cancer/IBD_HCC/210406")
load("mRNAS_coexp1.Rdata")
expr_df = res_all_adj
expr_df <- expr_df[expr_df$external_gene_name %in% RB$RB,]
expr_df$external_gene_name

gene_id = expr_df[3, "ensgenes"]
gene_name1 = expr_df[3,"external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
p5 <- ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name1))
ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "",size =5)
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)

gene_id = expr_df[7, "ensgenes"]
gene_name = expr_df[7, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
p6 <- ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "",size =5)

gene_id = expr_df[1, "ensgenes"]
gene_name = expr_df[1, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
p7 <-ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
 ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "",size =5)

gene_id = expr_df[4, "ensgenes"]
gene_name = expr_df[4, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
p8 <- ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "",size =5)


gene_id = expr_df[6, "ensgenes"]
gene_name1 = expr_df[6,"external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
p1 <- ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name1))
ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "",size =5)
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)



gene_id = expr_df[5, "ensgenes"]
gene_name = expr_df[5, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
p2 <- ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "",size =5)

gene_id = expr_df[2, "ensgenes"]
gene_name = expr_df[2, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
p3 <- ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "",size =5)

gene_id = expr_df[8, "ensgenes"]
gene_name = expr_df[8, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
p4 <- ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "",size =5)






RB <- c("CDK1","RFC4", "KIF14", "PRIM1","NEK2", "TOP2A", "RRM2", "CCNB1") 
RB <- data.frame(RB)
setwd("D:/Jeonju2021/Data/Research/Cancer/IBD_HCC/210406")
load("mRNAS_coexp1.Rdata")
expr_df = res_all_adj
expr_df <- expr_df[expr_df$external_gene_name %in% RB$RB,]
expr_df$external_gene_name

gene_id = expr_df[3, "ensgenes"]
gene_name1 = expr_df[3,"external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name1))
p5 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.715 \n p = 0.002",size =10)
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)

gene_id = expr_df[7, "ensgenes"]
gene_name = expr_df[7, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p6 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.769 \n p = 0.001",size =5)

gene_id = expr_df[1, "ensgenes"]
gene_name = expr_df[1, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p7 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 2.057 \n p < 0.0001",size =5)

gene_id = expr_df[4, "ensgenes"]
gene_name = expr_df[4, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p8 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.958 \n p = 0.0002",size =5)


gene_id = expr_df[6, "ensgenes"]
gene_name1 = expr_df[6,"external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name1))
p1 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.751 \n p = 0.0017",size =5)
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)



gene_id = expr_df[5, "ensgenes"]
gene_name = expr_df[5, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p2 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.812 \n p = 0.0009",size =5)

gene_id = expr_df[2, "ensgenes"]
gene_name = expr_df[2, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p3 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.627 \n p = 0.006",size =5)

gene_id = expr_df[8, "ensgenes"]
gene_name = expr_df[8, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p4 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.549 \n p = 0.013",size =5)









gridExtra::grid.arrange(p1, p2,p3, p4, ncol = 2)
gridExtra::grid.arrange(p5, p6,p7, p8, ncol = 2)

p1
p2
p3 
p4
p5
p6
p7
p8




gene_id = expr_df[3, "ensgenes"]
gene_name = expr_df[3, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")


fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p5 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.715 \n p = 0.002",size =5)











RB <- c("NEK2","KIF14","TOP2A","CCNB1","RFC4","CDK1","RRM2","PRIM1") 
RB <- data.frame(RB)
setwd("D:/Jeonju2021/Data/Research/Cancer/IBD_HCC/210406")
load("mRNAS_coexp1.Rdata")
expr_df = res_all_adj
expr_df <- expr_df[expr_df$external_gene_name %in% RB$RB,]
expr_df$external_gene_name





gene_id = expr_df[1, "ensgenes"]
gene_name = expr_df[1, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$NEK2 = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")

gene_id = expr_df[2, "ensgenes"]
gene_name = expr_df[2, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$KIF14 = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")

gene_id = expr_df[3, "ensgenes"]
gene_name = expr_df[3, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$TOP2A = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")

gene_id = expr_df[4, "ensgenes"]
gene_name = expr_df[4, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$CCNB1 = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")


gene_id = expr_df[5, "ensgenes"]
gene_name = expr_df[5, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$RFC4 = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")

gene_id = expr_df[6, "ensgenes"]
gene_name = expr_df[6, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$CDK1 = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")

gene_id = expr_df[7, "ensgenes"]
gene_name = expr_df[7, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$RRM2 = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")

gene_id = expr_df[8, "ensgenes"]
gene_name = expr_df[8, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$PRIM1 = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")

df=clin_df
df$deceased[df$deceased == "TRUE"] <- 1
df$deceased[df$deceased == "FALSE"] <- 2
str(df$deceased)

head(clin_df)
coxdata=df
coxdata$time=coxdata$overall_survival
coxdata$status=coxdata$deceased
lung=coxdata


head(lung)


covariates <- c("NEK2","KIF14","TOP2A","CCNB1","RFC4","CDK1","RRM2","PRIM1")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = lung)})
ggforest(model)
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         wald.test<-signif(x$wald["test"], digits=3)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=3);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)




covariates <- c("NEK2","KIF14","TOP2A","CCNB1","RFC4","CDK1","RRM2","PRIM1")

fit.coxph = coxph(Surv(overall_survival, deceased) ~ NEK2, data=clin_df)
fit.coxph = coxph(Surv(overall_survival, deceased) ~ KIF14, data=clin_df)
fit.coxph = coxph(Surv(overall_survival, deceased) ~ TOP2A, data=clin_df)
fit.coxph = coxph(Surv(overall_survival, deceased) ~ CCNB1, data=clin_df)
fit.coxph = coxph(Surv(overall_survival, deceased) ~ RFC4, data=clin_df)
fit.coxph = coxph(Surv(overall_survival, deceased) ~ CDK1, data=clin_df)
fit.coxph = coxph(Surv(overall_survival, deceased) ~ RRM2, data=clin_df)
fit.coxph = coxph(Surv(overall_survival, deceased) ~ PRIM1, data=clin_df)

cox.zph(fit.coxph)
cox.zph(fit.coxph, transform="km")


library(glmnet)
library(survival)

coxdata=lung

#coxdata[1:10,10:20]
#coxdata=coxdata[, 10:20]
head(coxdata)

coxdata<- na.omit(coxdata) 
res4<-glmnet(coxdata[, 10:18], Surv(coxdata$time, coxdata$status), family="cox", lambda='lambda.min')
coef(res4)







cvfit<-cv.glmnet(data.matrix(coxdata[,1:ncol(coxdata)]),Surv(coxdata$time, coxdata$status), family = "cox",alpha=1)
plot(cvfit)
c<-coef(cvfit ,s='lambda.min')

glmnet(data.matrix(coxdata[,1:ncol(coxdata)]),Surv(coxdata$time, coxdata$status), family = "cox")

res4<-glmnet(Z[, 1:30], Y, family="cox", lambda=0)
res4$beta[1:10,1]

predict(fit, type = "coef")

set.seed(1213)
N=100;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]
ty=rexp(N,exp(xb))
tcens=rbinom(n=N,prob=.3,size=1) # censoring indicator
y=cbind(time=ty,status=1-tcens)
fiti=Coxnet(x,y,penalty="Lasso") # Lasso











Coxnet=function(x, y, Omega=NULL, penalty=c("Lasso","Enet", "Net"), alpha=1, lambda=NULL, nlambda=50, rlambda=NULL, nfolds=1, foldid=NULL, inzero=TRUE, adaptive=c(FALSE,TRUE), aini=NULL, isd=FALSE, 
                ifast=TRUE, keep.beta=FALSE, thresh=1e-6, maxit=1e+5) 
  {#fcall=match.call()
  penalty=match.arg(penalty)
  if (penalty=="Lasso") {
    penalty="Enet"
    alpha=1
  }
  
  if (penalty=="Net" & is.null(Omega)) {
    penalty="Enet"
    cat("Enet was performed as no input of Omega")
  }
  
  fit=switch(penalty,
             "Enet"=coxEnet(x,y,alpha,lambda,nlambda,rlambda,nfolds,foldid,inzero,adaptive[1],aini,isd,ifast,keep.beta,thresh,maxit),
             "Net"=coxNet(x,y,Omega,alpha,lambda,nlambda,rlambda,nfolds,foldid,inzero,adaptive,aini,isd,ifast,keep.beta,thresh,maxit))
  
  #fit$call=fcall
  class(fit)="Coxnet"
  return(fit)
}







fit <- Coxnet(x = data.matrix(coxdata[,1:ncol(coxdata)]),
              y = Surv(coxdata$time, coxdata$status),penalty="Lasso")


install.packages("Coxnet")


set.seed(1213)
N=100;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]
ty=rexp(N,exp(xb))
tcens=rbinom(n=N,prob=.3,size=1)  # censoring indicator
y=cbind(time=ty,status=1-tcens)

fiti=Coxnet(x,y,penalty="Lasso",nlambda=10,nfolds=10) # Lasso
# attributes(fiti)





coef(fit, s = 0.05)

help(predict.glmnet)


predict(fit,newx=x[1:5,],s=c(0.01,0.005))
predict(fit,type="coef")




setwd("D:/Jeonju2021/Data/Research/Plosone/Revised/210416")
write.csv(lung,"hccsur2.csv")









library(h2o)
h2o.init()
lung.hex = as.h2o(survival::lung)

covariates = c("age", "sex",  "ph.karno", "ph.ecog", "wt.loss")
lung.coxph.model = h2o.coxph(x=covariates, event_column = "status", training_frame = lung.hex,
                             stop_column = "time")

# print model coefficients
lung.coxph.model@model$coefficients_table


res.cox <- coxph(Surv(time, status) ~ age + sex + ph.karno + ph.ecog + wt.loss, data =  lung)
summary(res.cox)




covariates = c("a", "b")
lung.coxph.model = h2o.coxph(x=covariates, event_column = "deceased", training_frame = df,
                             stop_column = "overall_survival")

# print model coefficients
lung.coxph.model@model$coefficients_table


res.cox <- coxph(Surv(time, status) ~ age + sex + ph.karno + ph.ecog + wt.loss, data =  lung)
summary(res.cox)





str(expr_df)
dir.create("D:/Jeonju2021/Data/Research/Plosone/Revised/210416")
setwd("D:/Jeonju2021/Data/Research/Plosone/Revised/210416")
write.csv(clin_df,"D:/Jeonju2021/Data/Research/ML/bank-full/clin_df2.csv")

dir.create("D:/Jeonju2021/Data/Research")




gene_id = expr_df[, "ensgenes"]
gene_name = expr_df[, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")


str(clin_df)



head(clin_df)

clin_df1 = clin_df["overall_survival", "deceased", "gene"]

fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)


covariates = c("KIF14", "RFC4", "CDK1" ,"PRIM1")


lung.coxph.model = h2o.coxph(x=clin_df$external_gene_name, event_column = clin_df$deceased, data=clin_df,
                             stop_column = overall_survival)








ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)


fit.coxph1 = coxph(Surv(overall_survival, deceased) ~ NEK2, data=clin_df)
fit.coxph2 = coxph(Surv(overall_survival, deceased) ~ KIF14, data=clin_df)
fit.coxph3 = coxph(Surv(overall_survival, deceased) ~ TOP2A, data=clin_df)
fit.coxph4 = coxph(Surv(overall_survival, deceased) ~ CCNB1, data=clin_df)
fit.coxph5 = coxph(Surv(overall_survival, deceased) ~ RFC4, data=clin_df)
fit.coxph6 = coxph(Surv(overall_survival, deceased) ~ CDK1, data=clin_df)
fit.coxph7 = coxph(Surv(overall_survival, deceased) ~ RRM2, data=clin_df)
fit.coxph8 = coxph(Surv(overall_survival, deceased) ~ PRIM1, data=clin_df)
library(multcomp) 
g1<-glht(fit.coxph1,linfct=mcp(NEK2="Tukey"),test=adjusted ("bonferroni"))
g2<-glht(fit.coxph2,linfct=mcp(KIF14="Tukey"),test=adjusted ("bonferroni"))
g3<-glht(fit.coxph3,linfct=mcp(TOP2A="Tukey"),test=adjusted ("bonferroni"))
g4<-glht(fit.coxph4,linfct=mcp(CCNB1="Tukey"),test=adjusted ("bonferroni"))
g5<-glht(fit.coxph5,linfct=mcp(RFC4="Tukey"),test=adjusted ("bonferroni"))
g6<-glht(fit.coxph6,linfct=mcp(CDK1="Tukey"),test=adjusted ("bonferroni"))
g7<-glht(fit.coxph7,linfct=mcp(RRM2="Tukey"),test=adjusted ("bonferroni"))
g8<-glht(fit.coxph8,linfct=mcp(PRIM1="Tukey"),test=adjusted ("bonferroni"))



g1<-summary(g1)$test$pvalues
g2<-summary(g2)$test$pvalues
g3<-summary(g3)$test$pvalues
g4<-summary(g4)$test$pvalues
g5<-summary(g5)$test$pvalues
g6<-summary(g6)$test$pvalues
g7<-summary(g7)$test$pvalues
g8<-summary(g8)$test$pvalues


c=cbind(g1,g2,g3,g4,g5,g6,g7,g8)

setwd("D:/Jeonju2021/Data/Research/Plosone/Revised/210416")
write.csv(c,"bonf.csv")

write.csv(res,"bonf1.csv")





RB <- c("NEK2", "TOP2A", "RRM2", "CCNB1") 
RB <- data.frame(RB)

load("mRNAS_coexp1.Rdata")
expr_df = res_all_adj
expr_df <- expr_df[expr_df$external_gene_name %in% RB$RB,]
expr_df = expr_df[order(expr_df$external_gene_name,expr_df$log2FoldChange),]

gene_name = expr_df[c(2,4,1,3),"external_gene_name"]
#write.csv(gene_name,"immune4.csv")

gene_id = expr_df[2, "ensgenes"]
gene_name1 = expr_df[2,"external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name1))
p1 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 0.674 \n p = 0.026",size =5)
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)



gene_id = expr_df[4, "ensgenes"]
gene_name = expr_df[4, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p2 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.475 \n p = 0.028",size =5)

gene_id = expr_df[1, "ensgenes"]
gene_name = expr_df[1, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p3 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.994 \n p = 0.0001",size =5)

gene_id = expr_df[3, "ensgenes"]
gene_name = expr_df[3, "external_gene_name"]
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]
median_value = median(clin_df$gene_value)
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)
ggsurv=ggsurvplot(fit, data=clin_df, pval = F,pval.method = F, risk.table=T, title=paste(gene_name))
fit.coxph = coxph(Surv(overall_survival, deceased) ~ gene, data=clin_df)
summary(fit.coxph)
glht_rx<-glht(fit.coxph,linfct=mcp(gene="Tukey"),test=adjusted("bonferroni"))
summary(glht_rx)
p4 <- ggsurv$plot +ggplot2::annotate("text",x = Inf, y = Inf,
                                     vjust = 1.2, hjust = 1.5,label = "HR = 1.459 \n p = 0.032",size =5)
gridExtra::grid.arrange(p1, p2,p3, p4, ncol = 2)
