# Code and Software statement

## BBT-BC-Serum sample-Prediction-Model

This is a Medical Prediction Model which can be used to predict the benignity or malignancy of the breast tumor patients through their serum.

### System requirements

The following package / library versions were used in this study:

python (version 3.9.18)
joblib (version  1.2.0)
numpy (version 1.26.0)
scikit-learn (version 1.2.1)

### Installation guide



### Demo



### Instructions for use




## WGCNA

Soft clustering of time series gene expression data

### System requirements

The following package / library versions were used in this study:

R (version 4.2.3)
WGCNA (version 1.72-1)

### Installation guide

install.packages("WGCNA")

### Demo
library(dplyr)
WGCNA_func<-function(dataExpr,type,corType,isislog, mingene){
  library(WGCNA)
  library(reshape2)
  library(stringr)
  options(stringsAsFactors = FALSE)
  robustY = ifelse(corType=="pearson",T,F)
  dataExpr <- t(dataExpr)
  if (isislog){
    dataExpr <- as.data.frame(log2(dataExpr+1))
  }

  
  gsg = goodSamplesGenes(dataExpr, verbose = 3)
  if (!gsg$allOK){
    dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
  }
  
  nGenes = ncol(dataExpr)
  nSamples = nrow(dataExpr)
  
  sampleTree = hclust(dist(dataExpr), method = "average")
  
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(dataExpr,powerVector= powers,
                          networkType=type, verbose=5)
  power = sft$powerEstimate
  pdf(file = paste('./','WGCNA_Softthreshold1_',type,corType,isislog,mingene,'.pdf'),onefile = FALSE)
  p1<-plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=1,col="red")
  abline(h=0.85,col="red")
print(p1)
dev.off()
pdf(file = paste('./','WGCNA_Softthreshold2_',type,corType,isislog,mingene,'.pdf'),onefile = FALSE)
p2<-plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=1, col="red")
print(p2)
dev.off()
  if (is.na(power)){
    power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                   ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                 ifelse(type == "unsigned", 6, 12))       
                   )
    )
  }
  cor <- WGCNA::cor
  net = blockwiseModules(dataExpr
                         ,power = power
                         , maxBlockSize = nGenes
                         , TOMType = type, minModuleSize = mingene
                         # ,pearsonFallback = Fallback
                         # ,networkType = network
                         , reassignThreshold = 0, mergeCutHeight = 0.25
                         , numericLabels = TRUE, pamRespectsDendro = FALSE
                         , saveTOMs=TRUE, corType = corType 
                         ,loadTOMs=TRUE
                         # ,deepSplit = Split
                         ,saveTOMFileBase = paste0('exprMat', ".tom")
                         , verbose = 3)
  cor<-stats::cor
  moduleLabels = net$colors
  moduleColors = labels2colors(moduleLabels)
  pdf(file = paste('./','WGCNA_clusterdendrogram_',type,corType,isislog,mingene,'.pdf'),onefile = FALSE)
  p3<-plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  print(p3)
  dev.off()
  MEs = net$MEs
  
  MEs_col = MEs
  colnames(MEs_col) = paste0("ME", labels2colors(
    as.numeric(str_replace_all(colnames(MEs),"ME",""))))
  MEs_col = orderMEs(MEs_col)
  pdf(file = paste('./','WGCNA_EigengeneNetworks_',type,corType,isislog,mingene,'.pdf'),onefile = FALSE)
  p4<-plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                        marDendro = c(3,3,2,4),
                        marHeatmap = c(3,4,2,2), plotDendrograms = T,
                        xLabelsAngle = 90)
print(p4)
dev.off()  

  TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
  load(net$TOMFiles[1], verbose=T)
  TOM <- as.matrix(TOM)
  dissTOM = 1-TOM
  plotTOM = dissTOM^7
  diag(plotTOM) = NA
  pdf(file = paste('./','WGCNA_plotTOM_',type,corType,isislog,mingene,'.pdf'),onefile = FALSE)
  p5<- TOMplot(plotTOM, net$dendrograms, moduleColors,main = "Network heatmap plot, all genes")
  print(p5)
  dev.off()  
  probes = colnames(dataExpr)
  dimnames(TOM) <- list(probes, probes)
  pdf(file = paste('./','WGCNA_cytoscape_',type,corType,isislog,mingene,'.pdf'),onefile = FALSE)
  cyt = exportNetworkToCytoscape(TOM,
                                 edgeFile = paste(dataExpr, ".edges.txt", sep=""),
                                 nodeFile = paste(dataExpr, ".nodes.txt", sep=""),
                                 weighted = TRUE, threshold = 0,
                                 nodeNames = probes, nodeAttr = moduleColors)
  print(cyt)
  dev.off()
  if (corType=="pearson") {
    modTraitCor = cor(MEs_col, clinical, use = "p")
    modTraitP = corPvalueStudent(modTraitCor, nSamples)
  } else {
    modTraitCorP = bicorAndPvalue(MEs_col, clinical, robustY=robustY)
    modTraitCor = modTraitCorP$bicor
    modTraitP   = modTraitCorP$p
  }
  write.csv(modTraitCor, paste('./',type,corType,isislog,mingene,'correlation.csv'))
  write.csv(modTraitP, paste('./',type,corType,isislog,mingene,'pvalue.csv'))
  gene_colors <- data.frame(genes = colnames(dataExpr),
                            color = moduleColors)
  print(gene_colors)
  write.csv(gene_colors, paste('./',type,corType,isislog,mingene,'gene_color.csv'))

if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}
write.csv(geneModuleMembership, paste('./',type,corType,isislog,mingene,'correlation1.csv'))
write.csv(MMPvalue, paste('./',type,corType,isislog,mingene,'pvalue1.csv'))

if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(dataExpr, clinical, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, clinical, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}
write.csv(geneTraitCor, paste('./',type,corType,isislog,mingene,'correlation2.csv'))
write.csv(geneTraitP, paste('./',type,corType,isislog,mingene,'pvalue2.csv'))

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
pdf(file = paste('./','WGCNA_',type,corType,isislog,mingene,'.pdf'),onefile = FALSE)
ax <- labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(clinical), 
                     yLabels = colnames(MEs_col), 
                     cex.lab = 0.5, 
                     ySymbols = colnames(MEs_col), colorLabels = FALSE, 
                     colors = blueWhiteRed(50), 
                     textMatrix = textMatrix, setStdMargins = FALSE, 
                     cex.text = 0.2, zlim = c(-0.3,0.3),
                     main = paste("Module-trait relationships"))
print(ax)
dev.off()
}

## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.

setwd("/Users/wgcna")

dataExpr <- read.csv('data.csv', row.names=1, header=T)
clinical <- read.csv('clinical.csv', row.names=1, header=T)

for (corType in c('pearson')){
  for (type in c('signed')){
      for (isislog in c(TRUE)){
          for (mingene in  c(50)){
            try(WGCNA_func(dataExpr,type = type,corType = corType,isislog = isislog,mingene=mingene),silent = F)
        }
      }
    }
}

### Instructions for use




## Mfuzz

Soft clustering of time series gene expression data

### System requirements

The following package / library versions were used in this study:

R (version 4.2.3)
Mfuzz (version 2.60.0)

### Installation guide

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Mfuzz")

### Demo

library("Mfuzz")
data(yeast)
yeast.r <- filter.NA(yeast, thres=0.25)
yeast.f <- fill.NA(yeast.r,mode="mean")
tmp <- filter.std(yeast.f,min.std=0)
yeast.s <- standardise(yeast.f)
cl <- mfuzz(yeast.s,c=16,m=1.25)
mfuzz.plot(yeast.s,cl=cl,mfrow=c(4,4),time.labels=seq(0,160,10))
m1 <- mestimate(yeast.s)
O <- overlap(cl)
Ptmp <- overlap.plot(cl,over=O,thres=0.05)
cl3 <- mfuzz(yeast.s,c=10,m=1.25)
mfuzz.plot(yeast.s,cl=cl3,mfrow=c(3,4))
O3 <- overlap(cl3)
overlap.plot(cl3,over=O3,P=Ptmp,thres=0.05)

### Instructions for use


