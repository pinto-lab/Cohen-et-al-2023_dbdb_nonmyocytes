
library(dplyr)
library(Seurat)
library(reticulate)
library(biomaRt)
library(ggplot2)
source("utilities.R")

#load object : object name=dbdb
load("C:/dbdb single cell DATA_cr5.0/dbdb_Seurat_Final_Object.RData")

# find markers for every cluster compared to all remaining cells, report only the positive ones
dbdb.markers <- FindAllMarkers(dbdb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dbdb.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
save(dbdb, dbdb.markers, file="DBDB2019")
write.csv(dbdb.markers,file="C:/dbdb single cell DATA_cr5.0/dbdb.markers.csv",row.names = FALSE,quote = FALSE)

#check for known marker gene expression
FeaturePlot(dbdb, "Ddah1", col=c('lightgray', 'red'), split.by='trt', sort.cell=TRUE,reduction="tSNE")
FeaturePlot(dbdb, "Cd74", col=c('lightgray', 'red'), split.by='trt', sort.cell=TRUE,reduction="tSNE")
FeaturePlot(dbdb, "Csf1r", col=c('lightgray', 'red'), split.by='trt', sort.cell=TRUE,reduction="tSNE")

#get top 10 markers for each cluster
top10 <-dbdb.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)


#add cluster annotaions - braod cell populations
new.cluster.ids <- c("Fibroblasts", "Fibroblasts", "Endocardial", "ECs", "Macrophages", "Fibroblasts", "B-cells", "Fibroblasts", "LECs", "Fibroblasts", "Pericytes", "Endocardial", "Macrophages", "Macrophages", "ECs", "T-cells", "ECs", "Fibroblasts", "DC-like", "NK cells", "SMCs", "T-cells", "Granulocytes", "Fibroblasts", "Fibroblasts", "Schwann cells", "Epicardial")
names(new.cluster.ids) <- levels(dbdb)
dbdb <- RenameIdents(dbdb, new.cluster.ids)


#function to extract top marker genes for each cell population
topNmarker <- function(object, markers, n){
  top100 <- markers %>% group_by(cluster) %>% top_n(100, p_val)
  topN<-NULL
  pops <- as.character(object@active.ident)
  for(i in 1:length(unique(pops))){
    topN.temp<-markers[markers$cluster==unique(top100$cluster)[i],][1:n,]
    topN<-rbind(topN, topN.temp)
  }
  topN<-topN[!is.na(topN$gene),]
  return(topN)
}

#extracting top 5 marker genes (ordered by p val)
top5 <- topNmarker(dbdb, dbdb.markers, 5)

# Top marker genes dotplot - ordered by pval
DotPlot(dbdb, features = unique(top5$gene),cols = c("#403891ff", "#f68f46ff"))+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=10))+
  xlab("Gene")+ylab("Cluster")+ theme(legend.position="top")
