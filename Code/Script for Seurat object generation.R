###################################################################### 
# Name: R script to perform preliminary analysis to generate 'dbdb' seurat object
###################################################################### 

library(dplyr)
library(Seurat)
library(reticulate)
library(biomaRt)
library(ggplot2)
source("utilities.R")

#load cellranger outputs into Seurat
dbhplus<- Read10X(data.dir = "C:/dbdb single cell DATA_cr5.0/cell ranger files/dbhplus_CR5_premrna/filtered_feature_bc_matrix")
dbhplus <- CreateSeuratObject(counts = dbhplus, project = "DbDb", min.cells = 3, min.features = 200)
dbhplus$trt<-"dbhplus"


dbdb<- Read10X(data.dir = "C:/dbdb single cell DATA_cr5.0/cell ranger files/dbdb_CR5_premrna/filtered_feature_bc_matrix")
dbdb<- CreateSeuratObject(counts = dbdb, project = "DbDb", min.cells = 3, min.features = 200)
dbdb$trt<-"dbdb"

# Data pre-processing  

#calculate mt percentage
dbhplus[["percent.mt"]] <- PercentageFeatureSet(dbhplus, pattern = "^mt-")
VlnPlot(dbhplus, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

dbdb[["percent.mt"]] <- PercentageFeatureSet(dbdb, pattern = "^mt-")
VlnPlot(dbdb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#merge suerat objects
samples <- c('dbhplus', 'dbdb')
dbdb <- merge(x=dbhplus, y=dbdb,add.cell.ids = samples, project="DbDb2019")

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(dbdb, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dbdb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#QC filtering
dbdb <- subset(dbdb, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10 & nCount_RNA > 700)

# Normalizing the data
dbdb<- NormalizeData(dbdb, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features (feature selection)
dbdb<- FindVariableFeatures(dbdb, selection.method = "vst", nfeatures = 2000)

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(dbdb), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(dbdb)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
CombinePlots(plots = list(plot2))

#cell-cycle regression
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl",host="apr2019.archive.ensembl.org")
mouse <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl",host="apr2019.archive.ensembl.org") 
genesmus <- getLDS(attributes="hgnc_symbol", filters="hgnc_symbol",
                   values=cc.genes$s.genes, mart=human,
                   attributesL="mgi_symbol", martL=mouse,
                   uniqueRows=TRUE)

# Add a column onto dat that gives mouse gene symbol,
# but do not include any gene without 1:1 orthology
genesmus2 <- dplyr::rename(genesmus, "human"="HGNC.symbol",
                           "mouse"="MGI.symbol") %>% group_by(human) %>%
  mutate(nmouse=n()) %>% group_by(mouse) %>%
  mutate(nhuman=n()) %>%
  dplyr::filter(nmouse == 1, nhuman == 1) %>%
  dplyr::select(-nmouse, -nhuman)

etrans <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "chromosome_name"),
                filter="mgi_symbol", values=genesmus2$mouse, mart=mouse)
genesmus3 <- inner_join(genesmus2, etrans, by=c("mouse"="mgi_symbol")) %>%
  dplyr::filter(!grepl("PATCH", chromosome_name)) %>%
  dplyr::filter(!grepl("CHR_MM", chromosome_name)) %>%
  dplyr::select(-chromosome_name)
s.genes<-genesmus3$mouse

mouse <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl",host="apr2019.archive.ensembl.org")  
genesmus <- getLDS(attributes="hgnc_symbol", filters="hgnc_symbol",
                   values=cc.genes$g2m.genes, mart=human,
                   attributesL="mgi_symbol", martL=mouse,
                   uniqueRows=TRUE)
# Add a column onto dat that gives mouse gene symbol,
# but do not include any gene without 1:1 orthology
genesmus2 <- dplyr::rename(genesmus, "human"="HGNC.symbol",
                           "mouse"="MGI.symbol") %>% group_by(human) %>%
  mutate(nmouse=n()) %>% group_by(mouse) %>%
  mutate(nhuman=n()) %>%
  dplyr::filter(nmouse == 1, nhuman == 1) %>%
  dplyr::select(-nmouse, -nhuman)

etrans <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "chromosome_name"),
                filter="mgi_symbol", values=genesmus2$mouse, mart=mouse)
genesmus3 <- inner_join(genesmus2, etrans, by=c("mouse"="mgi_symbol")) %>%
  dplyr::filter(!grepl("PATCH", chromosome_name)) %>%
  dplyr::filter(!grepl("CHR_MM", chromosome_name)) %>%
  dplyr::select(-chromosome_name)
g2m.genes<-genesmus3$mouse

#calculate cell cycle scores
dbdb<- CellCycleScoring(
  object = dbdb,
  g2m.features  = g2m.genes,
  s.features = s.genes,set.ident=TRUE)

#Scaling the data
dbdb<-ScaleData(dbdb)
dbdb <- RunPCA(dbdb, features = VariableFeatures(object = dbdb))
DimPlot(dbdb, reduction = "pca")

dbdb<- RunPCA(dbdb, features = c(s.genes, g2m.genes))
DimPlot(dbdb)

#remove unwanted variation by regressing out heterogeneity associated with cell cycle stage and mitochondial contamination.
dbdb<-ScaleData(dbdb,vars.to.regress=c("percent.mt","S.Score","G2M.Score"))
# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
dbdb<- RunPCA(dbdb, features = c(s.genes, g2m.genes))
DimPlot(dbdb)

#Perform linear dimensional reduction
dbdb <- RunPCA(dbdb, features = VariableFeatures(object = dbdb))
DimPlot(dbdb, reduction = "pca")

#Determine the ‘dimensionality’ of the dataset
ElbowPlot(dbdb, ndims= 50)

#select number of Pcs and resolusion
dbdb <- FindNeighbors(dbdb, dims = 1:40,verbose = FALSE)
dbdb <- FindClusters(dbdb, resolution = 1.2,verbose = FALSE)

#Run non-linear dimensional reduction tSNE for 2D visualization
dbdb <- RunTSNE(dbdb, reduction="pca",dims=1:40, tsne.method="Rtsne", reduction.name="tSNE", reduction.key="tSNE_",
  seed.use=1) 
DimPlot(dbdb, reduction = "tSNE",label=TRUE, label.size=5) + NoLegend()
