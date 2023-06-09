#############################################################################################
# Name : R script to perform fibroblasts analysis
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


#tsne visualization
DimPlot(dbdb, reduction = "tSNE",   
        cols.highlight = "#DE2D26",
        split.by = "trt", 
        pt.size = 1.3,
        cols = c("#C7E9C0", "#74C476", "orchid1", "dodgerblue1", "#FB6A4A", "#A1D99B", "mediumorchid2", "#41AB5D", "#636363", "#41AB5D", "firebrick", "orchid3", "#E6550D", "orangered", "dodgerblue3", "purple1", "lightslateblue", "#006D2C", "#FECC5C", "pink2", "maroon1", "plum3", "aquamarine", "#238B45", "#00441B", "darkorange", "hotpink"),
        ncol = 1,
        label = FALSE,
        repel = TRUE) + NoLegend() + theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())



### FIBROBLAST SUBSETTING AND MARKER GENE NALYSIS

Fibro_subset <- subset(x = dbdb,
                       idents = c("Fibroblast1", "Fibroblast2", "Fibroblast3", "Fibroblast4","Fibroblast5", "Fibroblast7","Fibroblast8", "Fibroblast9"))

dbdb_Fibro_subset <- subset(x = dbdb_alone, 
                      idents = c("Fibroblast1", "Fibroblast2", "Fibroblast3", "Fibroblast4","Fibroblast5", "Fibroblast7","Fibroblast8", "Fibroblast9"))

dbh_Fibro_subset <- subset(x = dbh_alone, 
                      idents = c("Fibroblast1", "Fibroblast2", "Fibroblast3", "Fibroblast4","Fibroblast5", "Fibroblast7","Fibroblast8", "Fibroblast9"))

#run find all markers function
dbdb_fibroblasts <- FindAllMarkers(Fibro_subset)


#Figure 4D
### merged violin plots ### 
Genes = c("Col1a1", "Col8a1", "Col14a1", "Meox1", "Postn", "Wif1", "Cilp", "Ccn2")
VlnPlot(Fibro_subset, features=Genes, pt.size = 0, stack=T, flip=T, 
        cols = c("tomato", "yellow2", "palevioletred2", "darkorchid3", "brown2", "springgreen4", "dodgerblue", "limegreen", "purple3")) +
  theme(legend.position = 'none')
