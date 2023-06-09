##################################################################################################################
# Name   : Script to generate dot plot for top 10 down-regulated genes
# Discription: This script will plot a dotplot for a merged list of top 10 DR genes identified in each cell population
#> with MASTcpmDetRate with 0.1 pct cutoff and pval<=0.001. 
##################################################################################################################

library(Seurat)
library(readxl)
library(openxlsx)
library(assertthat)
library(magrittr)
library(extrafont)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(dplyr)


#load data
load("C:/dbdb single cell DATA_cr5.0/dbdb_Seurat_Final_Object.RData")

# dbdb and dbh tSNE initial
DimPlot(dbdb, reduction = "tSNE",  pt.size = 1.5, label = TRUE) +NoLegend()

# Add cell name to meta data in "_cellname_"
dbdb@meta.data[["_cellname_"]] <- Idents(object = dbdb)

# Add cell name to meta data in "_cellnameBroad_"
dbdb@meta.data[["_cellnameBroad_"]] <- Idents(object = dbdb)

# Switch cell identity 
Idents(dbdb) <- "_cellname_"

# Extracting and Sorting the cluster names from the object
pops <- as.character(dbdb@active.ident)
my_levels <- sort(names(table(pops)))

# Leveling the cluster names in object based on sorted cluster names
dbdb@active.ident <- factor(x = dbdb@active.ident, levels = my_levels)

# Subset seurat object based on condition and save into a list (For fetching UR or DR data)
sub_object_list <- list()
for(i in 1:length(unique(dbdb@meta.data$trt))){
  name <- unique(dbdb@meta.data$trt)[i]
  sub_object_list[[name]] <- SubsetData(object = dbdb, subset.name = "trt",accept.value = name)
  }


#### Top10 UR/DR dot plot ####
# Variable for input DE data
path_de_datasheet = "C:/dbdb single cell DATA_cr5.0/DE_MAST_cpmDetRate_psignif_CR5/DE_MASTcmpDetRate_psignif_v6.xlsx"

# Variables for filter
n = 10 # top n genes
pct_nonzero = 0.1 # Non-zero percentage
p_value = 0.001 # p-value threshold
FC_threshold = 2 # FC threshold

# Variables for statistics 
# NOTE-1: FC = average expression condition1 / average expression condition_2 
#         (e.g. FC for UR_in_dbdb =  avg_expression_dbdb / avg_expression_dbhplus)
# NOTE-2: condition_1 and condition_2 need to be changed based on different regulation 
condition_1 = "dbhplus"  
condition_2 = "dbdb"
regulation = "DR_in_dbdb"  # UR or DR in condition_1
pct_nonzero_cond_1 = paste0("pct_nonzero_", condition_1) # Feild name of Non-zero percentage of condition 1
pct_nonzero_cond_2 = paste0("pct_nonzero_", condition_2) # Feild name of Non-zero percentage of condition 2
avg_expression_cond_1 = paste0("avg_expression_", condition_1) # Feild name of average expression of condition 1
avg_expression_cond_2 = paste0("avg_expression_", condition_2) # Feild name of average expression of condition 2

# object variable for fetching data based on different regulation 
object = sub_object_list[[condition_1]]

# set ggplot2 parameters for the dotplot 
x.lab.rot = T
plot.legend = T 
dot.scale = 8
do.return = T
scale.by = "radius"
col.min = 0
col.max = 1
scale.min = NA
scale.max = NA
dot.min = 0
dot.scale = 8
FC_dot_color ="blue"

# create a list of top 10 UR genes 
top10<-NULL
for(i in 1:length(names(table(pops)))){ 
  # Read in xlsx one sheet a time
  DE.results <- read_xlsx(path_de_datasheet, sheet=names(table(pops))[i])
  
  # Check if the cell name match with input DE table
  if(DE.results$cluster[1]==names(table(pops))[i]){
    txt = paste0(">>MATCH: ", DE.results$cluster[1], "==", names(table(pops))[i] )
    print(txt)
  }else if(DE.results$cluster[1]!=names(table(pops))[i]){
    txt = paste0("   >>", DE.results$cluster[1], "!=", names(table(pops))[i])
    print("ERROR: Cell name doesn't match!")
    print(txt)
  }else{
    txt = paste0("   >>", names(table(pops))[i])
    print("ERROR: Cell name not found!")
    print(txt)
  }
  
  # Selecting significant genes based on given gene percentage, p-value, regulation, and FC 
  DE.results <- DE.results[DE.results[[pct_nonzero_cond_1]]>= pct_nonzero|DE.results[[pct_nonzero_cond_2]]>= pct_nonzero,]
  DE.results <- DE.results[DE.results$pval <= p_value,]
  UR_or_DR.genes <- DE.results[DE.results[[regulation]] == 1,][1:n,]
  UR_or_DR.genes$cell_type <- names(table(pops))[i]
  UR_or_DR.genes$FC <- (UR_or_DR.genes[[avg_expression_cond_1]]/UR_or_DR.genes[[avg_expression_cond_2]])
  UR_or_DR.genes$FC[UR_or_DR.genes$FC > 6] <- 6
  UR_or_DR.genes <- UR_or_DR.genes[abs(UR_or_DR.genes$FC) >= FC_threshold,] #added FC therchold 21/01/2021
  top10<-rbind(top10,UR_or_DR.genes)
}

# merged list of top 10 UR genes in dbdb treated cells
UR_or_DR.gene.list<-unique(top10$gene)


genes.plot = sort(as.character(UR_or_DR.gene.list),decreasing = TRUE)

#scale.func to scale dot size
scale.func <- switch(EXPR = scale.by, size = scale_size, 
                     radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
# PercentAbove function to Calculate the percentage of a vector above some threshold
PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}

# create data frame for the dotplot
data.to.plot <- data.frame(FetchData(object = object, vars = genes.plot))
colnames(x = data.to.plot) <- genes.plot
data.to.plot$cell <- rownames( data.to.plot)
data.to.plot$id <- object@active.ident
data.to.plot <- data.to.plot %>% gather(key = genes.plot, 
                                        value = expression, -c(cell, id))
data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% 
  summarize(avg.exp = mean(expm1(x = expression)), pct.exp = PercentAbove(x = expression, 
                                                                          threshold = 0))
data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>% 
  mutate(avg.exp.scale = scale(x = avg.exp)) %>% mutate(avg.exp.scale = MinMax(data = avg.exp.scale, 
                                                                               max = col.max, min = col.min))
data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot, 
                                  levels = rev(x = genes.plot)) 

data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
data.to.plot$pct.exp <- data.to.plot$pct.exp * 100

data.to.plot<-data.to.plot[order(data.to.plot$id,decreasing = TRUE),]

colnames(data.to.plot)<-c("cell_type","gene","avg.exp","pct.exp","avg.exp.scale")

# obtain average gene expression of genes in UR_or_DR.gene.list in each cell population
top<-NULL
for(i in 1:length(names(table(pops)))){ 
  DE.results <- read_xlsx(path_de_datasheet, sheet=names(table(pops))[i])
  DE.results <- DE.results[DE.results[[pct_nonzero_cond_1]]>= pct_nonzero|DE.results[[pct_nonzero_cond_2]]>= pct_nonzero,]
  UR_or_DR.genes <- DE.results[DE.results$gene%in%UR_or_DR.gene.list,]
  UR_or_DR.genes$cell_type <- names(table(pops))[i]
  UR_or_DR.genes$FC <- (UR_or_DR.genes[[avg_expression_cond_1]]/UR_or_DR.genes[[avg_expression_cond_2]])
  UR_or_DR.genes$FC[UR_or_DR.genes$FC>6] <- 6
  top<-rbind(top,UR_or_DR.genes)
}

# merge 2 datasets 
data.to.plot$cell_type<-as.factor(data.to.plot$cell_type)
top$cell_type<-as.factor(top$cell_type)
data.to.plot.new <- data.to.plot %>% inner_join(top, by=c("gene","cell_type"))

# mark UR genes with p-value <0.001
#####data.to.plot.new$sig<-ifelse(data.to.plot.new$pval<=0.001,1,0)
data.to.plot.new$sig<-ifelse(data.to.plot.new$pval<=p_value & data.to.plot.new[[regulation]]==1, 1, 0) # corrected 21/01/2021

# plot data with ggplot2
p <-ggplot(data = data.to.plot.new, mapping = aes(x = gene, 
                                                 y = cell_type)) + 
  geom_point(mapping = aes(size = FC ,alpha=avg.exp.scale,stroke = 1),colour=FC_dot_color) +theme_bw()+
  scale.func(range = c(0, dot.scale), limits = c(scale.min, 
                                                 scale.max)) + theme(axis.title.x = element_blank(), 
                                                                     axis.title.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p <- p + geom_point(data=data.to.plot.new[data.to.plot.new$sig=="1",], aes(x=gene,y=cell_type),colour="red",size=0.4)

p <- p + theme(axis.text.x = element_text(angle = 90, 
                                          vjust = 1))+ theme(axis.text=element_text(size=10),
                                                               axis.title=element_text(size=10))
p


##################################################################################################################
