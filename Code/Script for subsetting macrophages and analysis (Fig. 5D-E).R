#############################################################################################
# Name : R script to perform macrophages analysis
#############################################################################################

library(tidyverse)
library(readxl)
library(dplyr)
library(openxlsx)
library(assertthat)
library(magrittr)
library(extrafont)
library(reshape2)
library(readr)
library(Seurat)
library(ggplot2)
library(devtools)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(viridis)

#load data
load("C:/dbdb single cell DATA_cr5.0/dbdb_Seurat_Final_Object.RData")
Idents(dbdb) <- "RNA_snn_res.1.2"
#rename clusters with sub population names
new.cluster.ids <- c("Fibroblast1", "Fibroblast2", "Endocardial1", "EC1", "Mac1", "Fibroblast3", "B-cells", "Fibroblast4", "LECs", "Fibroblast5", "Pericytes", "Endocardial2", "Mac2", "Mac3", "EC2", "T-cells", "EC3", "Fibroblast7", "DC-like", "NK-cells", "SMCs", "gd T-cells", "Granulocytes", "Fibroblast8", "Fibroblast9", "Schwann cells", "Epicardial")
names(new.cluster.ids) <- levels(dbdb)
dbdb <- RenameIdents(dbdb, new.cluster.ids)

#subset healthy and dbdb cells
dbdb_alone <- subset(x = dbdb, subset= (trt == "dbdb"))
dbh_alone <- subset(x = dbdb, subset = (trt == "dbhplus"))


# Macrophage subsetting and marker gene analysis
dbh_Mac_subset <- subset(x = dbh_alone,
                     idents = c("Mac1", "Mac2", "Mac3"))

dbdb_MacSubset <- subset(x = dbdb_alone,
                       idents = c("Mac1", "Mac2", "Mac3"))

merged_MacSubset <- subset(x = dbdb,
                           idents = c("Mac1", "Mac2", "Mac3"))



### cell number table ###
cellnumber_table <- as.data.frame.matrix(table(dbdb@active.ident, dbdb@meta.data$trt))
dbdb_cellnumber<- cellnumber_table %>% dplyr::select(dbdb)
dbh_cellnumber <- cellnumber_table %>% dplyr::select(dbhplus)
colnames(dbh_cellnumber) <- "cellnumber"
colnames(dbdb_cellnumber) <- "cellnumber"
dbdb_cellnumber$condition <- "dbdb"
dbh_cellnumber$condition <- "dbh"
dbdb_cellnumber$cluster <- row.names(dbdb_cellnumber)
dbh_cellnumber$cluster <- row.names(dbh_cellnumber)
row.names(dbdb_cellnumber) <- NULL
row.names(dbh_cellnumber) <- NULL


#Figure 5D
Genes1 = c("Csf1r", "Mrc1", "H2-Ab1", "C1qa", "C1qb", "Adgre1")
VlnPlot(merged_MacSubset, features=Genes1, pt.size = 0, stack=T, flip=T, 
        cols = c("tomato", "palevioletred2", "darkorchid3", "springgreen4", "green") +
          theme(legend.position = 'none')) +
  scale_fill_viridis(option = "A", discrete = T)



Genes2 = c("Csf1r", "Mrc1", "H2-Ab1", "Adgre1", "Apoe", "Tlr4" , "Abca1", "Sra1", "Acat1", "Msr1")
VlnPlot(dbh_Mac_subset, features=Genes2, pt.size = 0, stack=T, flip=T, 
        cols = c("tomato", "#74C476", "dodgerblue"))+
  scale_fill_viridis(option = "H", discrete = T) 


DimPlot(merged_MacSubset, reduction = "tSNE",label=TRUE, pt.size = 1.3, label.size=5, split.by = "trt",
        cols = c("tomato", "#74C476", "dodgerblue"))

DimPlot(merged_MacSubset, reduction = "tSNE",label=TRUE, pt.size = 1.7, label.size=5, 
        cols = c("#FB6A4A", "#E6550D", "orangered"))


#Subset Macrophages, find marker genes and perform GO analysis 

#run findallmarkers function
all.genes <- FindAllMarkers(merged_MacSubset)

Mac1_2_subset <-subset(x = dbh_alone,
                       idents = c("Mac1", "Mac2"))

Mac1_2 <- FindAllMarkers(Mac1_2_subset)

# Convert gene symbols to ENTREZ gene ids
#gene id conversion using org.Mm.eg.db

Mac1_2.df<-AnnotationDbi::select(org.Mm.eg.db, keys=Mac1_2$gene, columns='ENTREZID', keytype='SYMBOL')
all.genes.df<-AnnotationDbi::select(org.Mm.eg.db, keys=all.genes$gene, columns='ENTREZID', keytype='SYMBOL')

enrichGO.results<-enrichGO(
  Mac1_2.df$ENTREZID,OrgDb = "org.Mm.eg.db",
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
go.result<-as.data.frame(enrichGO.results@result)

#filter GO terms by min genes mapped
go.result<-go.result[go.result$Count>=2,]  

#save files
write.csv(go.result,paste("C:/dbdb single cell DATA_cr5.0/GO analysis/Macrophage GO terms", ".csv"))

#Figure 5E
barplot(enrichGO.results.2, color = "pvalue", order=TRUE, showCategory = 20)

write.csv(go.result, file =  "C:/dbdb single cell DATA_cr5.0/Macrophages in dbdb/GO analysis Macrophages in dbdb.csv")


