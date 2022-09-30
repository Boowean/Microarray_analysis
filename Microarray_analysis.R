library(GEOquery)
library(R.utils)
library(oligo)
library(oligoClasses)
library(tidyverse)
library(ggplot2)
library(affy)
library(hgu133plus2.db)
library(org.Hs.eg.db)
library(limma)
library(RColorBrewer)
gse50737 <- getGEO('GSE50737',GSEMatrix=TRUE,AnnotGPL=TRUE,destdir="./")
gse50737 <- gse50737[[1]]
#filePaths <- getGEOSuppFiles('GSE50737')
#untar("GSE50737/GSE50737_RAW.tar", list = TRUE)
#gse50737_celdata <- read.celfiles(list.celfiles('GSE50737',full.names=TRUE,listGzipped=TRUE))
varLabels(gse50737)

gse50737$supplementary_file
pd <- pData(gse50737)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1) 
gse50737_celdata <- read.celfiles(paste0('GSE50737/',pd$cel_file),phenoData=phenoData(gse50737))
#pData(gse50737)[["treatment:ch1"]]
pData(gse50737)[,c("geo_accession","treatment:ch1")]

hist(gse50737)
#oligo_normalised <- normalize(raw_nobg,method='quantile',which='pm')

gse50737_eset <- oligo::rma(gse50737_celdata)
design <- model.matrix( ~ gse50737_eset[['treatment:ch1']])
fit <- lmFit(gse50737_eset,design)
fitted.ebayes <- eBayes(fit)
topTable(fitted.ebayes)

interesting_genes <- topTable(fitted.ebayes,number=Inf,p.value = 0.05,lfc=2)
eset_of_interest <- gse50737_eset[rownames(interesting_genes),]
#heatmap(exprs(eset_of_interest))
heatmap(exprs(eset_of_interest),
        labCol=gse50737_eset[['culture medium:ch1']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))
