library(GenomicRanges)
library(cnvex)
library(tidyverse)
library(rlist)
library(gplots)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(mclust)

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
  time_length = 0
  significant_tiles = data.frame(0,0.0)
  significant_tiles = sig_tiles_total
  counter = 1
  for (i in 1:dim(tile_case)[2]){
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
    
    time_length = 0.3*(time_length) + 0.7*((Sys.time() - ptm) * (dim(tile_case)[2] - i)/60)
    
    print(paste("Tile processed: ",i,"/",dim(tile_case)[2],"------- Time remaining:",time_length,sep = ""))
  }
  #add column names
  colnames(significant_tiles) = c('tile', 'pvalue')
  #only keep the rows that has been filled up
  significant_tiles = significant_tiles %>% slice(1:counter)
  #sig_tiles_total = significant_tiles
  
  #tile_case_clust %>% ggplot(aes(x = `tile_case[, i]`, y = cluster))+geom_point()  #look at the position of points in clusters
  #look at the pvalues histogram
  significant_tiles %>% ggplot(aes(x = log10(pvalue)))+
    geom_density() + 
    xlim(-13,0)

  #look at the position of tiles with lowest pvalue
  #selected_tiles = significant_tiles %>% filter(significant_tiles$pvalue < 1e-8) %>% arrange(pvalue)
  #look at the position of tiles in the top 1% (in terms of lowest p-value)
  selected_tiles = significant_tiles %>% arrange(pvalue) %>% slice(1:round(dim(significant_tiles)[1]*.01))
  #top selected distribution on genome
  selected_tiles %>% ggplot(aes(x = tile))+
    geom_density() + 
    xlim(0,dim(tile_case)[2])+ theme_minimal()
    
  #find the position on chromosome
  tile_ranges = as.data.frame(ranges(file$tile))
  tile_ranges %>% filter(start < 30420200 & end > 30420200)

  
  
#start working on the bimodal idea
###########################
  #look at the distribution of tiles (the already finded to be significant)
  tile1 = as.data.frame(tile_case[,172])
  colnames(tile1) = "coverage"
  #for trimodal: tile 111968
  #for bimodal: tile 276
  #for unimodal: tile 1663
  tile1.gmm = Mclust(tile1)
  summary(tile1.gmm)
  tile1 = tile1 %>% mutate(cluster = tile1.gmm$classification)
  tile1$cluster = as.factor(tile1$cluster)
  tile1 %>% ggplot(aes(x = `coverage`, y = cluster))+geom_point()  #look at the position of points in clusters
  tile1 %>% ggplot(aes(x = coverage))+geom_density(alpha = 0.4)+theme_minimal()
  tile1 %>% ggplot(aes(x = coverage, fill = cluster))+geom_density(alpha = 0.4)+theme_minimal()
  #here we have two important information: likelihood for fitting unimodal, clusters. we also have the distributions too
###########################

  #interesting how results of this part, are tiles close to each other, that kindda match with hypothesis, take a look at their positions
#cluster the data using multimodal fitting
  time_length = 0
  #significant_tiles = data.frame(0,0.0)
  significant_tiles = sig_tiles_total %>% mutate(0) %>% mutate(1)
  counter = 1
  for (i in 1:dim(tile_case)[2]){
    ptm = Sys.time()
    if (sum(tile_case[,i]) > 0.1){ #if all the values are zero, or we only have 1 non-zero value, Mclust can't function
      try({k = Mclust(tile_case[,i],verbose = FALSE)})
      s
     
      tile_case_clust = as.data.frame(tile_case[,i]) %>% mutate(cluster = k$classification)
      #clust1 = tile_case_clust %>% filter(cluster == 1)
      #clust2 = tile_case_clust %>% filter(cluster == 2)
      #stat_test = t.test(clust1,clust2,var.equal = FALSE, paired = FALSE)
      
      #only keep ones with more than one component (multimodals)
      if (length(unique(k$classification)) > 1){
        significant_tiles[counter,] = c(i,k$loglik, k$bic, length(unique(k$classification)))
        counter = counter+1
      }
    }
    
    
    if (i %% 10 == 0){
      print(paste("Tile processed: ",i,"/",dim(tile_case)[2],"------- Time remaining:",time_length/10,sep = ""))
      time_length = 0
    }else{
      time_length = time_length + (Sys.time() - ptm) * (dim(tile_case)[2] - i)/60
    }
  }
  #add column names
  colnames(significant_tiles) = c('tile', 'loglik', "bic", 'modal_num')
  #only keep the rows that has been filled up
  significant_tiles = significant_tiles %>% slice(1:counter-1)
  #sig_tiles_total = significant_tiles
  
  #tile_case_clust %>% ggplot(aes(x = `tile_case[, i]`, y = cluster))+geom_point()  #look at the position of points in clusters
  #look at the pvalues histogram
  significant_tiles %>% ggplot(aes(x = loglik))+
    geom_density() + 
    xlim(-20,500)
  
  #look at the position of tiles with lowest pvalue
  #selected_tiles = significant_tiles %>% filter(significant_tiles$pvalue < 1e-8) %>% arrange(pvalue)
  #look at the position of tiles in the top 1% (in terms of lowest p-value)
  selected_tiles = significant_tiles %>% arrange(loglik) %>% slice(1:round(dim(significant_tiles)[1]*.01))
  #top selected distribution on genome
  selected_tiles %>% ggplot(aes(x = tile))+
    geom_density(alpha = 0.6) + 
    xlim(0,300000)+ theme_minimal()
  
significant_tiles_backup = significant_tiles

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