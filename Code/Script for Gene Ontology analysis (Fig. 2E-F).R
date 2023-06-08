###################################################################### 
# R script to perform Gene Ontology analysis for up and Down reglalated genes- Cohen et al 2023
###################################################################### 
library(clusterProfiler)
library(Seurat)
library(tidyverse)
library(readxl)
library(openxlsx)
library(assertthat)
library(magrittr)
library(extrafont)
library(reshape2)
library(org.Mm.eg.db)
library(AnnotationDbi)


#load data
load("C:/dbdb single cell DATA_cr5.0/dbdb_Seurat_Final_Object")
pops <- as.character(dbdb@active.ident)


###################################################################### 
#perform GO analysis for Up-regulated genes
###################################################################### 
#enrichGO R function 
#https://rdrr.io/bioc/clusterProfiler/src/R/enrichGO.R
# over reprentation analysis using hypergeometric test 
# - https://yulab-smu.github.io/clusterProfiler-book/chapter2.html#over-representation-analysis
#- https://yulab-smu.github.io/clusterProfiler-book/chapter5.html

#steps to perform DE analysis for each cell population
#1. load excel fils with DE results and apply pct threshold
#2.filter DE genes basd on Pvalue<0.01
#3.obtain a lit of UR(or DR) genes and their Entraz ids
#4.use enrich GO function to perform GO analysis


setwd("C:/dbdb single cell DATA_cr5.0/GO analysis")
Cluster.id<-names(table(pops))


for(i in 1:length(names(table(pops)))){
  DE.results<-read_xlsx("C:/dbdb single cell DATA_cr5.0/DE_MAST_cpmDetRate_psignif_CR5/DE_MASTcmpDetRate_psignif_v6.xlsx",sheet=i)
  # all genes considered in this dataset usually around 15,000-22,000
  # we use this all.genes as the backgroud gene list for GO over-representation analysis
  all.genes<-row.names(dbdb@assays$RNA)
  length(all.genes)
  DE.results<-DE.results[DE.results$pct_nonzero_dbhplus>=0.1|DE.results$pct_nonzero_dbdb>=0.1,]
  DE.results<-DE.results[DE.results$pval<=0.01,]
  DE.genes<-DE.results$gene[DE.results$UR_in_dbdb==1]
  print(names(table(pops))[i])
  print(DE.results$cluster[1])
  
  # Convert gene symbols to ENTREZ gene ids
  #gene id conversion using org.Mm.eg.db
  
  DE.genes.df<-AnnotationDbi::select(org.Mm.eg.db, keys=DE.genes, columns='ENTREZID', keytype='SYMBOL')
  all.genes.df<-AnnotationDbi::select(org.Mm.eg.db, keys=all.genes, columns='ENTREZID', keytype='SYMBOL')
  
  enrichGO.results<-enrichGO(
    DE.genes.df$ENTREZID,OrgDb = "org.Mm.eg.db",
    universe = all.genes.df$ENTREZID,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 0.5,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.5,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
  
  #use simplify to remove redundant GO terms
  #https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
  # GO terms that have semantic similarity higher than `cutoff` are treated as redundant terms
  # the results should be strip down before using simplify otherwise it takes too long to finish.
  #- https://github.com/YuLab-SMU/clusterProfiler/issues/28
  enrichGO.results.2 <- clusterProfiler::simplify (enrichGO.results,cutoff = 0.7, by = "p.adjust", select_fun = min)
  go.result<-as.data.frame(enrichGO.results.2@result)
  
  #filter GO terms by min genes mapped
  go.result<-go.result[go.result$Count>=2,]  
  
  #save files
  write.csv(go.result,paste("C:/dbdb single cell DATA_cr5.0/GO analysis",Cluster.id[i],".csv"))
  
} 

#all GO terms into a excel sheet
library(openxlsx)
GO_terms<- createWorkbook("DE.genes")

for(i in 1:length( Cluster.id)){
  
  
  go.results<-read_csv(paste("C:/dbdb single cell DATA_cr5.0/GO analysis",Cluster.id[i],".csv"))
  
  print(Cluster.id[i])
  addWorksheet(GO_terms,  sheetName=paste(Cluster.id[i]))
  
  writeData(GO_terms,sheet = paste(Cluster.id[i]),go.results,rowNames = FALSE)            
  saveWorkbook(GO_terms,"C:/dbdb single cell DATA_cr5.0/GO analysis/UR_GO_analysis.xlsx",overwrite = TRUE)
}


````
###################################################################### 
#perform GO analysis for Down-regulated genes
###################################################################### 

setwd("C:/dbdb single cell DATA_cr5.0/GO analysis")
Cluster.id<-names(table(pops))

for(i in 1:length(names(table(pops)))){
  DE.results<-read_xlsx("C:/dbdb single cell DATA_cr5.0/DE_MAST_cpmDetRate_psignif_CR5/DE_MASTcmpDetRate_psignif_v6.xlsx",sheet=i)
  # all genes considered in this dataset usually around 15,000-22,000
  # we use this all.genes as the backgroud gene list for GO over-representation analysis
  all.genes<-row.names(dbdb@assays$RNA)
  length(all.genes)
  DE.results<-DE.results[DE.results$pct_nonzero_dbhplus>=0.1|DE.results$pct_nonzero_dbdb>=0.1,]
  DE.results<-DE.results[DE.results$pval<=0.01,]
  DE.genes<-DE.results$gene[DE.results$DR_in_dbdb==1]
  print(names(table(pops))[i])
  print(DE.results$cluster[1])
  
  # Convert gene symbols to ENTREZ gene ids
  #gene id conversion using org.Mm.eg.db
  DE.genes.df<-AnnotationDbi::select(org.Mm.eg.db, keys=DE.genes, columns='ENTREZID', keytype='SYMBOL')
  all.genes.df<-AnnotationDbi::select(org.Mm.eg.db, keys=all.genes, columns='ENTREZID', keytype='SYMBOL')
  
  enrichGO.results<-enrichGO(
    DE.genes.df$ENTREZID,OrgDb = "org.Mm.eg.db",
    universe = all.genes.df$ENTREZID,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 0.5,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.5,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
  #use simplify to remove redundant GO terms
  #https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
  # GO terms that have semantic similarity higher than `cutoff` are treated as redundant terms
  # the results should be strip down before using simplify otherwise it takes too long to finish.
  #- https://github.com/YuLab-SMU/clusterProfiler/issues/28
  enrichGO.results.2 <- clusterProfiler::simplify (enrichGO.results,cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  go.result<-as.data.frame(enrichGO.results.2@result)
  
  #filter GO terms by min genes mapped
  go.result<-go.result[go.result$Count>=2,]
  
  write.csv(go.result,paste("C:/dbdb single cell DATA_cr5.0/GO analysis",Cluster.id[i],".csv"))
  
} 

#all GO terms into a excel sheet
library(openxlsx)
GO_terms<- createWorkbook("DE.genes")

for(i in 1:length( Cluster.id)){
  
  
  go.results<-read_csv(paste("C:/dbdb single cell DATA_cr5.0/GO analysis",Cluster.id[i],".csv"))
  
  print(Cluster.id[i])
  addWorksheet(GO_terms,  sheetName=paste(Cluster.id[i]))
  
  writeData(GO_terms,sheet = paste(Cluster.id[i]),go.results,rowNames = FALSE)            
  saveWorkbook(GO_terms,"C:/dbdb single cell DATA_cr5.0/GO analysis/DR_GO_analysis.xlsx",overwrite = TRUE)
}

