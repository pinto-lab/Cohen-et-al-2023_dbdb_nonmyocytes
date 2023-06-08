#############################################################################################
# Name : R script to perform differential gene expression analysis using MASTcpmDetRate
#############################################################################################

library(dplyr)
library(Seurat)
library(ranger)
library(assertthat)
library(magrittr)
library(methods)   
library(tidyverse)
library(edgeR)
library(MAST)
library(readxl)
library(openxlsx)

# load seurat object with major cluster identities
dbdb <- readRDS("C:/dbdb single cell DATA_cr5.0/dbdb_Seurat_Object.rds")
Idents(dbdb) <- "RNA_snn_res.1.2"

#rename clusters - add sub populations
new.cluster.ids <- c("Fibroblast1", "Fibroblast2", "Endocardial1", "EC1", "Mac1", "Fibroblast3", "B-cells", "Fibroblast4", "LECs", "Fibroblast5", "Pericytes", "Endocardial2", "Mac2", "Mac3", "EC2", "T-cells", "EC3", "Fibroblast7", "DC-like", "NK-cells", "SMCs", "gd T-cells", "Granulocytes", "Fibroblast8", "Fibroblast9", "Schwann cells", "Epicardial")
names(new.cluster.ids) <- levels(dbdb)
dbdb <- RenameIdents(dbdb, new.cluster.ids)

#table summary -number of cells per each population
pops <- as.character(dbdb@active.ident)
table(pops)

#defining treatment group in meta data
dbdb@meta.data$trt <- ifelse(dbdb@meta.data$trt == "dbhplus", "dbhplus", "dbdb")

#store cluster identities a meta data slot
dbdb@meta.data$mojor.clusters<-pops

#run_MASTcpmDetRate  :function to perform DE analysis 

run_MASTcpmDetRate <- function(L) {
  timing <- system.time({
    stopifnot(all(names(L$grp) == colnames(L$count)))
    cdr <- scale(colMeans(L$count > 0))  # scales exactly with nGene
    dge <- edgeR::DGEList(counts = L$count)
    dge <- edgeR::calcNormFactors(dge)
    cpms <- edgeR::cpm(dge)
    # pass the data into MAST. cData is cell data
    sca <- FromMatrix(exprsArray = log2(cpms + 1), 
                      cData = data.frame(wellKey = names(L$grp), 
                                         grp = L$grp, cdr = cdr))
    zlmdata <- zlm(~cdr + grp, sca)  # ~1.5 min for a cluster with ~400 cells
    mast <- lrTest(zlmdata, "grp")
  })
  list(timing = timing,
       #res = mast,
       result.details = "hurdle model",
       df = data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"],
                       row.names = names(mast[, "hurdle", "Pr(>Chisq)"])))
}

#ectract cluster lables for each cell from the seurat object
clusters <- rownames_to_column(dbdb@meta.data[, "major.clusters", drop=FALSE], "cell") %>% 
  deframe()
uclust <- unique(unname(clusters)) %>% sort()

#extract count data matrix from the seurat object
counts <- as.data.frame(dbdb@assays$RNA@counts) 
counts<-counts[, unlist(dbdb@assays$RNA@data@Dimnames[2])]

#extract treatment group information from the seurat object
cond <- dbdb@meta.data[, "trt", drop=FALSE] %>% rownames_to_column("cell") %>%
  deframe()

#perform DE analysis using run_MASTcpmDetRate function
DE_results <- list()
for (i in 1:length(uclust)) {
  cells <- names(clusters)[clusters == uclust[i]]
  grp <- cond[cells]
  X <- counts[, cells]
  DE_results[[i]] <- run_MASTcpmDetRate(list(grp=grp, count=X))
  DE_results[[i]]$cluster <- uclust[i]
  print(uclust[i])
}


p <- do.call('cbind', lapply(DE_results, "[[", "df", drop=F))
colnames(p) <- sapply(DE_results, "[[", "cluster")
p_tidy <- rownames_to_column(p, "gene") %>%
  gather("cluster", "pval", -gene)
p_signif <- p_tidy %>% arrange(pval) 
group_by(p_tidy,cluster)

#saved as a csv file
write.csv(p_signif,"C:/dbdb single cell DATA_cr5.0/DE_MAST_cpmDetRate_psignif_CR5/DE_MASTcmpDetRate_psignif_v9.csv")
#read csv file
p_signif<-read.csv("C:/dbdb single cell DATA_cr5.0/DE_MAST_cpmDetRate_psignif_CR5/DE_MASTcmpDetRate_psignif_v9.csv")[,-1]

#calculate average gene expression gene
genes <- unique(p_signif$gene)
avg <- AverageExpression(dbdb, features=genes)
nonzero <- apply(counts, 1, function(vals) 
  tapply(vals, dbdb@meta.data$major.clusters, function(x) mean(x > 0))) %>% t()

avg.df<-as.matrix(avg$RNA)
colnames(avg.df)<-c(as.character(colnames(avg.df)))
avg.df.melt<-reshape2::melt(avg.df)
colnames(avg.df.melt)<-c("gene","cluster","avg_expression_overall")

nonzero<-as.matrix(nonzero)
colnames(nonzero)<-c(as.character(colnames(nonzero)))
nonzero.melt<-reshape2::melt(nonzero)
colnames(nonzero.melt)<-c("gene","cluster","pct_nonzero_overall")


nonzero.melt[cbind(p_signif$gene, p_signif$cluster),]

p_signif<-inner_join(p_signif,avg.df.melt,by=c("gene","cluster"))
p_signif<-inner_join(p_signif,nonzero.melt,by=c("gene","cluster"))


cond <- rownames_to_column(dbdb@meta.data[, "trt", drop=F], "cell") %>% deframe()
for (name in c("dbhplus", "dbdb")) {
  dat <- SubsetData(dbdb, cells=names(cond)[cond == name])
  avg <- AverageExpression(dat, features=genes)
  
  counts <- as.data.frame(dat@assays$RNA@counts) 
  counts<-counts[, unlist(dat@assays$RNA@data@Dimnames[2])]
  
  nonzero <- apply(counts, 1, function(vals) 
    tapply(vals, dat@meta.data$major.clusters, function(x) mean(x > 0))) %>% t()
  
  
  avg.df<-as.matrix(avg$RNA)
  colnames(avg.df)<-c(as.character(colnames(avg.df)))
  avg.df.melt<-reshape2::melt(avg.df)
  colnames(avg.df.melt)<-c("gene","cluster",paste("avg_expression_",name,sep=""))
  
  nonzero<-as.matrix(nonzero)
  colnames(nonzero)<-c(as.character(colnames(nonzero)))
  nonzero.melt<-reshape2::melt(nonzero)
  colnames(nonzero.melt)<-c("gene","cluster",paste("pct_nonzero_",name,sep=""))
  
  p_signif<-inner_join(p_signif,avg.df.melt,by=c("gene","cluster"))
  p_signif<-inner_join(p_signif,nonzero.melt,by=c("gene","cluster"))
}


#write De results into a dataframe
DE.table1<-NULL
for(i in 1:length(uclust)){
  DE.table<-p_signif[p_signif$cluster==paste(uclust[i]),]
  #applying pct cut-off; percentage of cells expressing a given gene
  DE.table<-DE.table[DE.table$pct_nonzero_dbhplus>=0.1|DE.table$pct_nonzero_dbdb>=0.1,]
  #adjust for multiple hypothesis testing
  DE.table$BH.pvalue<-p.adjust( DE.table$pval,method = "BH")
  #apply p value cut-off to identify DE genes
  DE.table<-DE.table[DE.table$pval<=0.05,]
  DE.table1<-rbind(DE.table,DE.table1)
}
p_signif1<-DE.table1

#write DE results into a csv file
write.csv(p_signif,"C:/dbdb single cell DATA_cr5.0/DE_MAST_cpmDetRate_psignif_CR5/DE_MASTcmpDetRate_psignif_v9.csv",row.names = FALSE)

towrite <- c( split(p_signif1, p_signif1$cluster))

#comple csv files into a multiple-sheet excel file
library(openxlsx)
DE.genes<- createWorkbook("DE.genes")
for(i in 1:length(uclust)){
  DE.table<-towrite[[i]]
  DE.table$log2FC<-log2(DE.table$avg_expression_dbdb/DE.table$avg_expression_dbhplus)
  DE.table$UR_in_dbdb<-ifelse(DE.table$log2FC>0,1,0)
  DE.table$DR_in_dbdb<-ifelse(DE.table$log2FC<0,1,0)
  DE.table<-DE.table[DE.table$pct_nonzero_dbhplus>=0.1|DE.table$pct_nonzero_dbdb>=0.1,]
  addWorksheet(DE.genes,  sheetName=paste(DE.table$cluster[i]))
  writeData(DE.genes,sheet = paste(DE.table$cluster[i]),DE.table,rowNames =FALSE)            
  saveWorkbook(DE.genes,"C:/dbdb single cell DATA_cr5.0/DE_MAST_cpmDetRate_psignif_CR5/DE_MASTcmpDetRate_psignif_v10.xlsx",overwrite = TRUE)
  }
