############################################################################################  
# Name: R script to perform bulk-RNA-seq analyses
############################################################################################ 

library(rJava)
library(UpSetR)
library(tidyverse)
library(venneuler)
library(grid)
library(reshape)
library(ggplot2)
library(readxl)
library(scales)
library(eulerr)
library(readr)
library(VennDiagram)
library("ggvenn")
library("ggVennDiagram")
library(enrichplot)
library(DOSE)

##import datasets into global environment###
Kesherwani <- read.csv("C:/dbdb single cell DATA_cr5.0/Bulk RNA seq/Kesherwani_GSE66575_RNAseq_genes_expression_original.csv")
View(Kesherwani)

dbdb_scRNAseq <- read.csv("C:/dbdb single cell DATA_cr5.0/DE_MAST_cpmDetRate_psignif_CR5/Original sub-clustered DE table ALSO FOR BULK/For_venn_DE_gene_list_wo_Rps.csv")
View(dbdb_scRNAseq)

dbdb_bulk <- read_csv("C:/dbdb single cell DATA_cr5.0/Bulk RNA seq/dbdb_bulk_added_pval.csv")
View(dbdb_bulk)

#for GO analysis#
raw_Akita <- Kesherwani$gene
write.csv(raw_Akita, file = "C:/dbdb single cell DATA_cr5.0/Bulk RNA seq/Bulk RNA DE gene lists/Total_Akita_genes.csv")
raw_dbdb_bulk <- dbdb_bulk$GeneName
write.csv(raw_dbdb_bulk, file = "C:/dbdb single cell DATA_cr5.0/Bulk RNA seq/Bulk RNA DE gene lists/raw_dbdb_bulk.csv")
raw_dbdb_scRNAseq <- dbdb_scRNAseq$gene
write.csv(raw_dbdb_scRNAseq, file = "C:/dbdb single cell DATA_cr5.0/Bulk RNA seq/Bulk RNA DE gene lists/Total_Akita_genes.csv")


##filter gene list based on pval (P<0.001)##
scRNAseq_dbdb_filtered <- filter(dbdb_scRNAseq, pval<=0.001)
Kesherwani_filter <- filter(Kesherwani, p_value<=0.01)
dbdb_bulk_filtered <- filter(dbdb_bulk, pvalue<=0.01)

##remove redundant entries/genes##
scRNAseq_dbdb_gene_list <- as.character(unique(scRNAseq_dbdb_filtered$gene))
Akita_gene_list <- as.character(unique(Kesherwani_filter$gene))
dbdb_bulk_gene_list <- as.character(unique(dbdb_bulk_filtered$GeneName))

overlap=calculate.overlap(
  x=list(
    "A"=c(scRNAseq_dbdb_gene_list),
    "B"=c(Akita_gene_list),
    "C"=c(dbdb_bulk_gene_list)
  )
)

x <- list(
  A = sample(scRNAseq_dbdb_gene_list), 
  B = sample(Akita_gene_list), 
  C = sample(dbdb_bulk_gene_list)
)

ggvenn(
  x, 
  fill_color = c("dodgerblue", "olivedrab2", "springgreen"),
  stroke_size = 0.5, set_name_size = 4
)

ggVennDiagram(x, label_alpha = 0.5,
              category.names = c("dbdb scRNAseq", "Akita bulk", "dbdb bulk"))

overlapping_list1 <- intersect(Akita_gene_list, scRNAseq_dbdb_gene_list)
overlapping_list2 <- intersect(Akita_gene_list, dbdb_bulk_gene_list)
overlapping_list3 <- intersect(dbdb_bulk_gene_list, scRNAseq_dbdb_gene_list)

write.csv(overlapping_list1, file = "C:/dbdb single cell DATA_cr5.0/Bulk RNA seq/Overlapping_gene_files/Akita&scRNAseqdbdb_intersect_v2.csv")
write.csv(overlapping_list2, file = "C:/dbdb single cell DATA_cr5.0/Bulk RNA seq/Overlapping_gene_files/Akita_vs_dbdb_BULK_v2.csv")
write.csv(overlapping_list3, file = "C:/dbdb single cell DATA_cr5.0/Bulk RNA seq/Overlapping_gene_files/dbdb_bulk_vs_scRNAseq_intersect_v2.csv")

##overlapping between all datasets##
common_diabetes_genes <- print(overlap$a5)
write.csv(overlap$a5, file = "C:/dbdb single cell DATA_cr5.0/Bulk RNA seq/Overlapping_gene_files/Major_intersect_v2.csv")

# Set the chart data
expressionInput <- c('Akita' = Akita_gene_list, 'dbdb_scRNAseq' = scRNAseq_dbdb_gene_list,
                      'dbdb_bulk' = dbdb_bulk_gene_list) 
# Create the Venn diagram
plot(ggvenn(x, fill_color = "dodgerblue", "coral", "springgreen", stroke_size = 0.5, set_name_size = 3))


##values MUST be proportional for this package to work correctly, see below for proportional values of above raw vals###
##A = dbdb scRNAseq, B = Akita, C = dbdb bulk##
myExpVenn <- venneuler(c(A = 0.56,B = 0.12, C = 0.24, 
                         "A&B" = 0.03, "A&C" = 0.04,
                         "B&C" = 0.010,
                         "A&B&C" = 0.01))

## above values must add up to '1' for the venn to work##
par(cex=1.2)
plot(myExpVenn, main = "Common genes in the diabetic heart")
grid.text(
  "Cohen et al 2021",
  x = 0.9,
  y = 0.2,
  gp = gpar(
    fontsize = 10,
    fontface = 3
  )
) 


write.csv(Akita_gene_list, file="C:/dbdb single cell DATA_cr5.0/Bulk RNA seq/Overlapping_gene_files/Akita_gene_list.csv", row.names = FALSE)
write.csv(dbdb_bulk_gene_list, file="C:/dbdb single cell DATA_cr5.0/Bulk RNA seq/Overlapping_gene_files/dbdb_bulk_gene_list.csv", row.names = FALSE)
write.csv(scRNAseq_dbdb_gene_list, file="C:/dbdb single cell DATA_cr5.0/Bulk RNA seq/Overlapping_gene_files/scRNAseq_dbdb_gene_list.csv", row.names = FALSE)

########GO ANALYSIS ON OVERLAPPING GENES####
knitr::opts_chunk$set(echo = TRUE)
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
library(dplyr)

#assign overlaps to name#
##overlapping genes##
common_diabetes_genes <- overlap$a5
Akita_vs_dbdbscRNAseq <- overlap$a2
Akita_vs_dbdb_bulk <- overlap$a6
dbdbscRNAseq_vs_dbdb_bulk <- overlap$a4
##unique genes##
dbdb_scRNAseq_unique_genes <- overlap$a1
Akita_unique_genes <- overlap$a3
dbdb_bulk_unique_genes <- overlap$a7

# Convert gene symbols to ENTREZ gene ids
#gene id conversion using org.Mm.eg.db

raw_Akita2 <- as.data.frame(as.character(raw_Akita))
raw_dbdb_bulk2 <- as.data.frame(raw_dbdb_bulk)
raw_dbdb_scRNAseq2 <- as.data.frame(raw_dbdb_scRNAseq)

names(raw_Akita2) <- "Genes"
names(raw_dbdb_bulk2) <- "Genes"
names(raw_dbdb_scRNAseq2) <- "Genes"

all.genes<-rbind(raw_Akita2, raw_dbdb_bulk2, raw_dbdb_bulk2)
all.genes$Genes <- as.character(all.genes$Genes)

all.genes <- (unique(all.genes))

Akita_vs_dbdb_bulk.df<-AnnotationDbi::select(org.Mm.eg.db, keys=Akita_vs_dbdb_bulk, columns='ENTREZID', keytype='SYMBOL')
all.genes.df<-AnnotationDbi::select(org.Mm.eg.db, keys=all.genes$Genes, columns='ENTREZID', keytype='SYMBOL')

enrichGO.results<-enrichGO(
  Akita_vs_dbdb_bulk.df$ENTREZID,OrgDb = "org.Mm.eg.db",
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

#use simplify to remove redundant GO term
#https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
# GO terms that have semantic similarity higher than `cutoff` are treated as redundant terms
# the results should be strip down before using simplify otherwise it takes too long to finish.
#- https://github.com/YuLab-SMU/clusterProfiler/issues/28

enrichGO.results.2 <- clusterProfiler::simplify (enrichGO.results,cutoff = 0.7, by = "p.adjust", select_fun = min)

go.result<-as.data.frame(enrichGO.results.2@result)


# enrichGO.results.2@result$Count <- filter(enrichGO.results.2@result$Count >= 2)
enrichGO.results.2_filter2 <- gsfilter(enrichGO.results.2, by = 'Count', min = 2)


#filter GO terms by min genes mapped
go.result<-go.result[go.result$Count>=2,]  

#save files
write.csv(go.result,paste("C:/dbdb single cell DATA_cr5.0/Bulk RNA seq/Wo_Rps_genes_GO_analysis_of_all_genes/GO analysis_Akita_vs_dbdb_bulk.csv"))

### CHANGING COLOURS OF BARPLOT FOR COMMON GO TERMS IN THE DIABETIC HEART ###

barplot.enrichResult <- function(height, x="Count", color='p.adjust', showCategory=20, font.size=12, title="", ...) {
  ## use *height* to satisy barplot generic definition
  ## actually here is an enrichResult object.
  object <- height
  
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
  }
  else if (x == "count" || x == "Count") {
    x <- "Count"
  }
  
  df <- fortify(object, showCategory=showCategory, by=x, ...)
  
  if(colorBy %in% colnames(df)) {
    p <- ggplot(df, aes_string(x = "Description", y = x, fill = colorBy)) +
      theme_dose(10) +
      scale_fill_continuous(low="indianred2", high="gray", name = "p-value", guide=guide_colorbar(reverse=TRUE))
  } else {
    p <- ggplot(df, aes_string(x = "Description", y = x, fill = "Description")) +
      theme_dose(10) +
      theme(legend.position="none")
  }
  p + geom_bar(stat = "identity",width = 0.75)+ ylim(0, 10) + coord_flip() +
    ggtitle(title) + xlab(NULL) + ylab(NULL)+theme_bw()+ theme(axis.title.x = element_blank(), 
                                                               axis.title.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12))
  
}

barplot(enrichGO.results.2_filter2, color = "pvalue", order=TRUE, showCategory = 20)
 

