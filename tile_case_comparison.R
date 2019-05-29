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
  tile1 = as.data.frame(tile_case[,252709])
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
  selected_tiles %>% ggplot(aes(x = tile, fill = as.factor(modal_num)))+
    geom_histogram(binwidth = 100,alpha = 0.6,) + 
    xlim(0,dim(significant_tiles)[1])+ theme(line = )
  
significant_tiles_backup = significant_tiles

#look at the position of results in genome
  #make a data frame of all the tile positions:
  tile_pos = as.data.frame(file$tile) %>% select(seqnames,start,end,strand,target,arm,seg)
  tile_pos = tile_pos %>% mutate(tile = row_number())
  
  tile_modal = left_join(tile_pos, significant_tiles, by = "tile")
  
  tile_modal %>% filter(seqnames == 'chr8' & modal_num > 5)
  #get each chromosome tile boundary:
  tile_modal %>% group_by(seqnames) %>% summarise(min(tile),max(tile))
  
  #plot the positions
  (tile_modal %>% filter(modal_num > 1)) %>% ggplot()+
    geom_histogram(binwidth = 200,alpha = 0.6,aes(x = tile, fill = ..count.. > 20)) + 
    xlim(0,287509)+ theme_minimal()+scale_fill_manual(values=c("#FFFFFF", "red"))+
    
    
    geom_vline(xintercept = 24896,linetype = "dashed",size = 0.1)+ geom_text(x = 24896, y = 0, label = "chr1", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 49116,linetype = "dashed",size = 0.1)+ geom_text(x = 49116, y = 0, label = "chr2", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 68946,linetype = "dashed",size = 0.1)+ geom_text(x = 68946, y = 0, label = "chr3", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 87968,linetype = "dashed",size = 0.1)+ geom_text(x = 87968, y = 0, label = "chr4", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 106122,linetype = "dashed",size = 0.1)+ geom_text(x = 106122, y = 0, label = "chr5", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 123203,linetype = "dashed",size = 0.1)+ geom_text(x = 123203, y = 0, label = "chr6", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 139138,linetype = "dashed",size = 0.1)+ geom_text(x = 139138, y = 0, label = "chr7", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 153652,linetype = "dashed",size = 0.1)+ geom_text(x = 153652, y = 0, label = "chr8", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 167492,linetype = "dashed",size = 0.1)+ geom_text(x = 167492, y = 0, label = "chr9", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 180872,linetype = "dashed",size = 0.1)+ geom_text(x = 180872, y = 0, label = "chr10", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 194381,linetype = "dashed",size = 0.1)+ geom_text(x = 194381, y = 0, label = "chr11", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 207709,linetype = "dashed",size = 0.1)+ geom_text(x = 207709, y = 0, label = "chr12", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 219146,linetype = "dashed",size = 0.1)+ geom_text(x = 219146, y = 0, label = "chr13", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 229851,linetype = "dashed",size = 0.1)+ geom_text(x = 229851, y = 0, label = "chr14", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 240051,linetype = "dashed",size = 0.1)+ geom_text(x = 240051, y = 0, label = "chr15", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 249085,linetype = "dashed",size = 0.1)+ geom_text(x = 249085, y = 0, label = "chr16", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 257411,linetype = "dashed",size = 0.1)+ geom_text(x = 257411, y = 0, label = "chr17", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 265449,linetype = "dashed",size = 0.1)+ geom_text(x = 265449, y = 0, label = "chr18", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 271311,linetype = "dashed",size = 0.1)+ geom_text(x = 271311, y = 0, label = "chr19", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 277756,linetype = "dashed",size = 0.1)+ geom_text(x = 277756, y = 0, label = "chr20", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 282427,linetype = "dashed",size = 0.1)+ geom_text(x = 282427, y = 0, label = "chr21", size = 4, angle = 90, vjust = -0.4, hjust= 0)+
    geom_vline(xintercept = 287509,linetype = "dashed",size = 0.1)+ geom_text(x = 287509, y = 0, label = "chr22", size = 4, angle = 90, vjust = -0.4, hjust= 0)
#    +geom_vline(xintercept = 303114,linetype = "dashed",size = 0.1)+ geom_text(x = 303114, y = 0, label = "chrX", size = 4, angle = 90, vjust = -0.4, hjust= 0)
#    +geom_vline(xintercept = 308837,linetype = "dashed",size = 0.1)+ geom_text(x = 308837, y = 0, label = "chrY", size = 4, angle = 90, vjust = -0.4, hjust= 0)


#adding fragile sites to the previous plot:
  fragile = read_tsv("/home/noshadh/Codes/Germline_variation_detection/Fragile_sites", col_names = FALSE)
  fragile = fragile %>% select(X1,X2,X3,X4)
  colnames(fragile) = c('chr','start','end','site_name')
  #map to tiles
  fragile = fragile %>% mutate(tile_start = 0) %>% mutate(tile_end = 0)
  for (i in 1:dim(fragile)[1]){
    #here, for every member of fragile, we will find the tile. O(N^2), BAD BAD BAD !!!
    chr = fragile$chr[i]
    chr_tiles = tile_modal %>% filter(seqnames == chr)
    
    start_tile_on_chr = round(fragile[i,]$start/10000)+1
    end_tile_on_chr = round(fragile[i,]$end/10000)
    
    if (fragile[i,]$end > max(chr_tiles$end)){
      end_tile_on_chr = dim(chr_tiles)[1]
    }
    fragile[i,5] = chr_tiles %>% slice(start_tile_on_chr) %>% select(tile)
    fragile[i,6] = chr_tiles %>% slice(end_tile_on_chr) %>% select(tile)
    
    print(paste("Fragile site number",i,"has been processed..."))
  }
  
  #till here has been tested and works
  
  
  #convert the positions, to 
  s_positions = fragile %>% select(pos = tile_start) %>% mutate(type = "start")
  e_positions = fragile %>% select(pos = tile_end) %>% mutate(type = "end")
  positions = rbind(s_positions, e_positions)
  
  p = ggplot()
  for (i in 5:7){
    p = p + geom_segment(aes(x=s_positions$pos[i],xend=e_positions$pos[i],y=1,yend=1))
  }
  p = p + geom_segment(aes(x=s_positions$pos[i],xend=e_positions$pos[i],y=18,yend=18)) + geom_segment(aes(x=s_positions$pos[i+6],xend=e_positions$pos[i+6],y=18,yend=18))
  p = p + geom_segment(aes(x=s_positions$pos[i+23],xend=e_positions$pos[i+23],y=18,yend=18))
  p
  ggplot()+geom_segment(aes(x=s_positions$pos,xend=e_positions$pos,y=18,yend=18))
    geom_vline(xintercept = s_positions$pos, aes(fill = 'red'))+
    geom_vline(xintercept = e_positions$pos)












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





  #idea: 
#cluster for each tile
#test if significantly different
#if yes, output which samples in each cluster, p-value, mean of each cluster
#make sure the tiles are from targeted sites