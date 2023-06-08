#### Below is the script for generating figure shown in Figure 3A, performed using output generated using CellPhoneDb.  ###


library(pheatmap)


heatmaps_plot = function(count_matrix_ofn, meta_file, pvalues_file, count_filename, log_filename, count_network_filename, interaction_count_filename, count_network_separator, interaction_count_separator, show_rownames = T, show_colnames = T,
                         scale="none", cluster_cols = T, border_color='white', cluster_rows = T, fontsize_row=11,
                         fontsize_col = 11, main = '',treeheight_row=0, family='Arial', treeheight_col = 0,
                         col1 = "", col2 = '', col3 = '', meta_sep=',', pvalues_sep='\t', pvalue=0.001){

  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)
  
  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  intr_pairs = all_intr$interacting_pair
  # all_intr = all_intr[,-c(1:11)]
  
  
  split_sep = '\\|'
  join_sep = '|'
  
  pairs1_all = unique(meta[,2])
  # print(pairs1_all)
  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
      pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
  new_pair1 =NULL
  for(i in 1:length(pairs1)){
    if(length(unlist(gregexpr('c2_', pairs1[i]))) !=2 && unlist(gregexpr('c2_', pairs1[i]))>0){
      new_pair1 = rbind(new_pair1, pairs1[i])
    }}
  new_pair1 <- as.character(new_pair1)
  pairs1 = new_pair1
  
  tissue_heart =NULL
  heart_tissue =NULL
  for(i in 1:length(pairs1)){
    if(unlist(gregexpr('c2_', pairs1[i])) > 1){
      tissue_heart = rbind(tissue_heart, pairs1[i])
    }else{
      heart_tissue = rbind(heart_tissue, pairs1[i])
    }
  }
  tissue_heart <- as.character(tissue_heart)
  heart_tissue <- as.character(heart_tissue)
  
  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  
  pairs1 = tissue_heart
  tissue_heart_table = NULL
  heart_tissue_table = NULL
  for(i in 1:length(pairs1))
  {tryCatch({
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]

    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    print(pairs1[i])
    tissue_heart_intr = all_intr[which(all_intr[,pairs1[i]]<=pvalue),]
    tissue_heart_intr = tissue_heart_intr[which(tissue_heart_intr$receptor_a=="False" &
                                                  tissue_heart_intr$receptor_b=="True"),]
    tissue_heart_intr = tissue_heart_intr[,c(1:11, which(colnames(tissue_heart_intr)==pairs1[i]))]
    tissue_heart_intr$cell_pair = pairs1[i]
    names(tissue_heart_intr)[names(tissue_heart_intr) == pairs1[i]] <- 'pval'
    tissue_heart_table = rbind(tissue_heart_table, tissue_heart_intr)
    
    n1 = tissue_heart_intr$interacting_pair
    
    pairs_rev = paste(p2, p1, sep=join_sep)
    print(pairs_rev)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    
    heart_tissue_intr = all_intr[which(all_intr[,pairs_rev]<=pvalue),]
    heart_tissue_intr = heart_tissue_intr[which(heart_tissue_intr$receptor_a=="True" &
                                                  heart_tissue_intr$receptor_b=="False"),]
    heart_tissue_intr = heart_tissue_intr[,c(1:11, which(colnames(heart_tissue_intr)==pairs_rev))]
    heart_tissue_intr$cell_pair = pairs_rev
    names(heart_tissue_intr)[names(heart_tissue_intr) == pairs_rev] <- 'pval'
    heart_tissue_table = rbind(heart_tissue_table, heart_tissue_intr)
    
    n2 = heart_tissue_intr$interacting_pair
    
    if(p1!=p2)
      count1 = length(unique(n1))+length(unique(n2))
      # count1 = length(unique(c(n1, n2)))
    else
      count1 = length(unique(n1))
    
    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  },error=function(e){})
  }
  
  all_count = all_count[-1,]
  row.names(all_count) = NULL
  all_count = as.data.frame(all_count)
  all_count$pair = paste0(all_count$SOURCE,"|",all_count$TARGET)
  missing_pair = setdiff(tissue_heart, all_count$pair)
  missing_pair = as.data.frame(missing_pair)
  library(dplyr)
  library(tidyr)
  missing_pair = missing_pair %>% separate(missing_pair, c("SOURCE", "TARGET"), sep ='\\|')
  missing_pair$count = 0
  missing_pair$pair = paste0(missing_pair$SOURCE,"|",missing_pair$TARGET)
  all_count = rbind(all_count, missing_pair)
  write.csv(all_count, count_matrix_ofn,  row.names = F)
  write.csv(heart_tissue_table, "//fs1/home/IHsu/Projects/R/2022/proj_adiposeCellulome/data/cellphonedb/wat_dbdb_heart/stringent_heart_tissue_table.csv",  row.names = F)
  write.csv(tissue_heart_table, "//fs1/home/IHsu/Projects/R/2022/proj_adiposeCellulome/data/cellphonedb/wat_dbdb_heart/stringent_tissue_heart_table.csv",  row.names = F)
  
 
  
  
 # write.table(all_count, count_network_filename, sep=count_network_separator, quote=F, row.names = F)
  
  #######   count interactions
  
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))
    
  }
  
  
  if (any(count1)>0)
  {
    count_matrix = matrix(count1, nrow=length(unique(meta[,2])), ncol=length(unique(meta[,2])))
    rownames(count_matrix)= unique(meta[,2])
    colnames(count_matrix)= unique(meta[,2])
    
    all_sum = rowSums(count_matrix)
    all_sum = cbind(names(all_sum), all_sum)
    # write.table(all_sum, file=interaction_count_filename, quote=F, sep=count_network_separator, row.names=F)
    # write.table(count_matrix, count_matrix_ofn, na = "", row.names = TRUE, col.names = TRUE, sep=",")
    count_matrix <- read.csv(count_matrix_ofn, sep=",",check.names=F)
    count_matrix <- data.matrix(count_matrix)
    # print(count_matrix)
    col.heatmap <- colorRampPalette(c(col1,col3 ))( 1000 )
    count_matrix[lower.tri(count_matrix)]<-NA
    
    pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col,filename = count_filename)
    
    # pheatmap(log(count_matrix+1), show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             # border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             # main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col)
    
    
    } else {
    stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
    }
  # return(count_matrix)
}


#################

# meta_file = "./data/cellphoneDB_example_data/dbdb_metadata_table.csv"
# meta_file = "./data/cellphoneDB_example_data/dbhplus_metadata_table.csv"
meta_file = "//fs1/home/IHsu/Projects/R/2022/proj_adiposeCellulome/data/cellphonedb/wat_dbdb_heart/220606_watdbdb_heart_metadata_table.csv"

# pvalues_file = "./data/cellphoneDB_example_data/dbdb_pvalues.txt"
# pvalues_file = "./data/cellphoneDB_example_data/dbh_pvalues.txt"
pvalues_file = "//fs1/home/IHsu/Projects/R/2022/proj_adiposeCellulome/data/cellphonedb/wat_dbdb_heart/pvalues.txt"


count_fileneame = "//fs1/home/IHsu/Projects/R/2022/proj_adiposeCellulome/data/cellphonedb/wat_dbdb_heart/stringent_heatmap.pdf"

count_matrix_ofn = "//fs1/home/IHsu/Projects/R/2022/proj_adiposeCellulome/data/cellphonedb/wat_dbdb_heart/stringent_count_matrix.csv"

# col1="#c13c81ff"
col1="white"
# col2="#9324a3ff"
col2=""
# col3="#0d1687ff"
col3="orangered"
count_matrix <- heatmaps_plot(count_matrix_ofn, meta_file, pvalues_file, count_fileneame, cluster_cols = F, 
                             col1=col1, col3 = col3, cluster_rows =F)
heatmaps_plot(count_matrix_ofn, meta_file, pvalues_file, count_fileneame, cluster_cols = F, 
              col1=col1, col3 = col3, cluster_rows =F)

# count_matrix <- read.csv(diff_count_matrix_ofn, sep=",",check.names=F)
