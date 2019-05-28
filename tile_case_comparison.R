library(GenomicRanges)
library(cnvex)
library(tidyverse)
library(rlist)
library(gplots)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

set.seed(1024)

#set the file directories
setwd("/mctp/users/mcieslik/proj/study/cptac3/data/cnvex-hl5/")

#get a list of folders
files = list.files()

#results file
tile_case = list()
colnames_list = ""
#reading each file, process it, save the n.cov information, delete the file
for (file_num in 1:length(files)-1){
  file = readRDS(paste("./",files[file_num],"/",files[file_num],".rds",sep = ""))
  tiles = as.data.frame(file$tile)
  tile_case[file_num] = tiles %>% select(n.cov)
  colnames_list[file_num] = files[file_num]
  print(paste("File processed: ",file_num,"/",length(files),sep = ""))
}

#convert a list to data frame
tile_case = as.data.frame(matrix(unlist(tile_case),ncol = length(tile_case), byrow = FALSE))
colnames(tile_case) = colnames_list

#write the file and data frames
write.table(tile_case, "/home/noshadh/Codes/Germline_variation_detection/Tile_case_hl5.tsv", sep = "\t", col.names = TRUE)

#look at the density of coverages for each segment
as.data.frame(tile_case %>% select(`C3N-01522-hl5`)) %>%
  ggplot(aes(x = `C3N-01522-hl5`))+
  geom_density()+
  xlim(0,1)+ylim(0,5)

#make a Hierarchical Clustering heatmap
tile_case = t(tile_case)
heatmap.2(as.matrix(t(tile_case[229000:230000,])),density.info="none", trace="none")




#cluster the data using kmeans for 2 clusters
  significant_tiles = data.frame(0,0.0)
  counter = 1
  for (i in 156804:dim(tile_case)[2]){
    ptm = Sys.time()
    try({k = kmeans(tile_case[,i],centers = 2, nstart = 100)})
    tile_case_clust = as.data.frame(tile_case[,i]) %>% mutate(cluster = k$cluster)
    clust1 = tile_case_clust %>% filter(cluster == 1)
    clust2 = tile_case_clust %>% filter(cluster == 2)
    stat_test = t.test(clust1,clust2,var.equal = FALSE, paired = FALSE)
    if (stat_test$p.value < 0.000001){
      significant_tiles[counter,] = c(i,stat_test$p.value)
      counter = counter+1
    }
    if (i == 1){
      time_length = Sys.time() - ptm  
    }
    
    print(paste("Tile processed: ",i,"/",dim(tile_case)[2],"------- Time remaining:",time_length * (dim(tile_case)[2] - i)/60,sep = ""))
  }
  colnames(significant_tiles) = c('tile', 'pvalue')
  
  sig_tiles_total = rbind(significant_tiles,significant_tiles2)
  
  #look at the pvalues histogram
  sig_tiles_total %>% ggplot(aes(x = log10(pvalue)))+
    geom_density() + 
    xlim(-11,0)

  #look at the position of tiles with lowest pvalue
  selected_tiles = sig_tiles_total %>% filter(sig_tiles_total$pvalue < 1e-6) %>% arrange(pvalue)
  #top selected distribution on genome
  selected_tiles %>% ggplot(aes(x = tile))+
    geom_density() + 
    xlim(0,dim(tile_case)[2])+ theme_minimal()
    
  #find the position on chromosome
  tile_ranges = as.data.frame(ranges(file$tile))
  tile_ranges %>% filter(start < 30420200 & end > 30420200)

    






#note: here we will find best number of cluster and use anova for statistics value
#cluster the data using kmeans
for (i in 1:dim(tile_case)[2]){
  #finding best number of clusters using silhouette
  best_clust_num = fviz_nbclust(as.data.frame(tile_case[,i]), kmeans, method = "silhouette")
  best_clust_num = as.numeric(best_clust_num$data %>% filter(y == max(best_clust_num$data$y)) %>% select(clusters))
  k = kmeans(tile_case[,i],centers = best_clust_num, nstart = 100)
  
  print(paste("Tile processed: ",i,"/",dim(tile_case)[2],sep = ""))
}
euci_dist = get_dist(tile_case[,i])

library(pheatmap)
annotation_sam_pheatmap=data.frame(cluster=factor(k$cluster))
pheatmap(as.matrix((tile_case[,i])),annotation=t(annotation_sam_pheatmap),cluster_cols=F,fontsize_row=5,fontsize_col=3)


#optimal number of clustering

fviz_nbclust(as.data.frame(tile_case[,i]), kmeans, method = "silhouette")

gap_stat <- clusGap(as.data.frame(tile_case[,i]), FUN = kmeans, nstart = 100,
                    K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

#temp
data_tile = as.data.frame(tile_case[,i]) %>% mutate(clustering = k$cluster)




  #idea: 
#cluster for each tile
#test if significantly different
#if yes, output which samples in each cluster, p-value, mean of each cluster
#make sure the tiles are from targeted sites