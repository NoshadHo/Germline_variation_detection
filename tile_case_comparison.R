library(GenomicRanges)
library(cnvex) #this is a local library, might not work on other machines
library(tidyverse)
library(rlist)
library(gplots)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(mclust)
library(gridExtra)
library(foreach)
library(doParallel)
set.seed(1024)
numCores = detectCores()
registerDoParallel(numCores-1)

#all the data are saved in K-means_Multimodal_fitting_data.RData, K-means data (first step) are available in Data_first.RData
##READ THE FILES-------------------------------------------------------------------------------------------------------------
#set the file directories
setwd("/mctp/users/mcieslik/proj/study/cptac3/data/cnvex-hl5/")
#get a list of folders
files = list.files()
#results file
tile_case = list()
colnames_list = ""
#reading each file, process it, save the n.cov information, delete the file
#here we use list, to improve the performance and get results faster
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

##INITIAL ANALYSIS-------------------------------------------------------------------------------------------------------------
#make a Hierarchical Clustering heatmap
#the heatmap is so big, we need to look at the samll chunk of tiles
tile_case = t(tile_case)
dev.off()
heatmap.2(as.matrix(t(tile_case[22900:23000,])),density.info="none", trace="none")

##K-MEANS METHOD---------------------------------------------------------------------------------------------------------------
#cluster the data using kmeans for 2 clusters
  time_length = 0
  
  #use any dataframe with the dimentions of 300000      2 (or whatever number of tiles we have)
  significant_tiles = as.data.frame(matrix(data = c(0,0),nrow = 300000,ncol = 2)) #here I made a predefined data.frame to improve the performance, but it is not defined yet
  counter = 1 #counter for the output file
  for (i in 1:dim(tile_case)[2]){
    ptm = Sys.time()
    try({k = kmeans(tile_case[,i],centers = 2, nstart = 100)})
    tile_case_clust = as.data.frame(tile_case[,i]) %>% mutate(cluster = k$cluster) #add corresponding clusters for each sampleto the data.frame
    clust1 = tile_case_clust %>% filter(cluster == 1)
    clust2 = tile_case_clust %>% filter(cluster == 2)
    stat_test = t.test(clust1,clust2,var.equal = FALSE, paired = FALSE) #find the significance
    if (stat_test$p.value < 0.000001){ #select only significant ones
      significant_tiles[counter,] = c(i,stat_test$p.value)
      counter = counter+1
    }
    time_length = 0.3*(time_length) + 0.7*((Sys.time() - ptm) * (dim(tile_case)[2] - i)/60) #estimate the remainig time, FUN :))
    print(paste("Tile processed: ",i,"/",dim(tile_case)[2],"------- Time remaining:",time_length,sep = ""))
  }
  #add column names to data.frame
  colnames(significant_tiles) = c('tile', 'pvalue')
  #only keep the rows that has been filled up
  significant_tiles = significant_tiles %>% slice(1:counter)
  #Optional analysis that can be used for each tile (quality control)
  tile_case_clust %>% ggplot(aes(x = `tile_case[, i]`, y = cluster))+geom_point()  #look at the position of points in clusters
  
  #look at the pvalues histogram
  significant_tiles %>% ggplot(aes(x = log10(pvalue)))+
    geom_density() + 
    xlim(-13,0)
  #look at the position of tiles in the top 1% (in terms of lowest p-value)
  selected_tiles = sig_tiles_total %>% arrange(pvalue) %>% slice(1:round(dim(significant_tiles)[1]*0.01))
  
  #PLOTING
  #top selected distribution on genome
  selected_tiles %>% ggplot()+
    geom_histogram(binwidth = 100,aes(x = tile, fill = "#01e0cf"))+
    xlim(0,287509)+ theme_classic()+
    geom_vline(xintercept = 24896,linetype = "dashed",size = 0.3)+
    geom_text(x = 24896, y = -1, label = "chr1", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 49116,linetype = "dashed",size = 0.3)+
    geom_text(x = 49116, y = -1, label = "chr2", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 68946,linetype = "dashed",size = 0.3)+
    geom_text(x = 68946, y = -1, label = "chr3", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 87968,linetype = "dashed",size = 0.3)+
    geom_text(x = 87968, y = -1, label = "chr4", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 106122,linetype = "dashed",size = 0.3)+
    geom_text(x = 106122, y = -1, label = "chr5", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 123203,linetype = "dashed",size = 0.3)+
    geom_text(x = 123203, y = -1, label = "chr6", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 139138,linetype = "dashed",size = 0.3)+
    geom_text(x = 139138, y = -1, label = "chr7", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 153652,linetype = "dashed",size = 0.3)+
    geom_text(x = 153652, y = -1, label = "chr8", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 167492,linetype = "dashed",size = 0.3)+
    geom_text(x = 167492, y = -1, label = "chr9", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 180872,linetype = "dashed",size = 0.3)+
    geom_text(x = 180872, y = -1, label = "chr10", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 194381,linetype = "dashed",size = 0.3)+
    geom_text(x = 194381, y = -1, label = "chr11", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 207709,linetype = "dashed",size = 0.3)+
    geom_text(x = 207709, y = -1, label = "chr12", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 219146,linetype = "dashed",size = 0.3)+
    geom_text(x = 219146, y = -1, label = "chr13", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 229851,linetype = "dashed",size = 0.3)+
    geom_text(x = 229851, y = -1, label = "chr14", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 240051,linetype = "dashed",size = 0.3)+
    geom_text(x = 240051, y = -1, label = "chr15", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 249085,linetype = "dashed",size = 0.3)+
    geom_text(x = 249085, y = -1, label = "chr16", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 257411,linetype = "dashed",size = 0.3)+
    geom_text(x = 257411, y = -1, label = "chr17", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 265449,linetype = "dashed",size = 0.3)+
    geom_text(x = 265449, y = -1, label = "chr18", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 271311,linetype = "dashed",size = 0.3)+
    geom_text(x = 271311, y = -1, label = "chr19", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 277756,linetype = "dashed",size = 0.3)+
    geom_text(x = 277756, y = -1, label = "chr20", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 282427,linetype = "dashed",size = 0.3)+
    geom_text(x = 282427, y = -1, label = "chr21", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 287509,linetype = "dashed",size = 0.3)+
    geom_text(x = 287509, y = -1, label = "chr22", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))
    

  
##MULTIMODAL METHOD---------------------------------------------------------------------------------------------------------------

#start working on the bimodal idea, test if it works for one tiles
#whenever we want to get informations about one specific tile, we should use this chunk
###########################
  #look at the distribution of tiles (the already finded to be significant)
  tile1 = as.data.frame(tile_case[,74014])
  colnames(tile1) = "coverage"
  #for trimodal: tile 111968 in lh5
  #for bimodal: tile 276 in lh5
  #for unimodal: tile 1663 in lh5
  tile1.gmm = Mclust(tile1) #multimodal fitting
  summary(tile1.gmm)
  tile1 = tile1 %>% mutate(cluster = tile1.gmm$classification) #add clusters to data.frame
  tile1$cluster = as.factor(tile1$cluster)
  tile1 %>% ggplot(aes(x = `coverage`, y = cluster))+geom_point()  #look at the position of points in clusters
    
  tile1 %>% ggplot(aes(x = coverage))+geom_density(alpha = 0.4)+theme_minimal() #look at the distribution of coverage in a specific tile as a whole

  tile1 %>% ggplot(aes(x = coverage, fill = cluster))+geom_density(alpha = 0.4)+theme_minimal() #look at the distribution of coverage in a specific tile for each modal
###########################
#if you want to try some other dataset instead of tile_case, do this:
  #-it should be in the form of row:variables(features) col:tiles or samples
  #and then run this:
  right_form_of_data = as.data.frame(t(purified_tile_cov_gc_normalized))
  tile_case = right_form_of_data
  rm(right_form_of_data,purified_tile_cov_gc_normalized,rotated_tile_cov_gc_normalized, eigen_vectors)
  #FITTING FOR ALL THE SAMPLES
  #interesting how results of this part, are tiles close to each other, that kindda match with hypothesis, take a look at their positions :)
  #cluster the data using multimodal fitting
  PTIME = system.time({
  significant_tiles = foreach(i = 1:dim(tile_case)[2]) %dopar% {
    if (sum(tile_case[,i]) > 0.1){ #if all the values are zero, or we only have 1 non-zero value, Mclust can't function
      try({k = Mclust(tile_case[,i],verbose = FALSE)})
      #only keep ones with more than one component (multimodals)
      if (length(unique(k$classification)) > 1){
        (c(i,k$loglik, k$bic, length(unique(k$classification)), min(table(k$classification))))
      }
    }
  }
  })
  significant_tiles = as.data.frame(do.call(rbind,significant_tiles))
  #add column names
  colnames(significant_tiles) = c('tile', 'loglik', "bic", 'modal_num', 'min_clust_size')
  dim(significant_tiles)
  #have a back up after a time-consuming loop
  significant_tiles_2to4pc = significant_tiles
  #look at the pvalues histogram
  significant_tiles %>% ggplot(aes(x = loglik))+
    geom_density() + 
    xlim(-20,500)+ggtitle("loglik after GC correction")
  #look at the position of tiles in the top 1% (in terms of lowest p-value)
  selected_tiles = significant_tiles %>% arrange(loglik) %>% slice(1:round(dim(significant_tiles)[1]*.01))
  #}

  
  #write the selected tiles:
  write.table(significant_tiles, "/home/noshadh/Codes/Germline_variation_detection/Selected_volatile_tiles.tsv", sep = "\t", col.names = TRUE)
    
#GET READY FOR MULTIMODAL PLOT AND THEN PLOTING ---------------------------------------------------------------------------
#make a general data.frame-join information we have, list of selected tiles, with their information on tile_case data.frame
  
  #look at the position of results in genome
  #make a data frame of all the tile positions: #here we choose one of the files (original files containing Rds info, cnvex output)
  tile_pos = as.data.frame(file$tile) %>% select(seqnames,start,end,strand,target,arm,seg) #extracting information of each tile
  tile_pos = tile_pos %>% mutate(tile = row_number()) #add the indexing of each tile (indexing in the file is relative to chromosomes)
  
  tile_modal = left_join(tile_pos, significant_tiles, by = "tile")
  
  #get each chromosome tile boundary:
  tile_modal %>% group_by(seqnames) %>% summarise(min(tile),max(tile)) #it is optional, we will use it for setting chr boundaries in plot
  
  #make the tile modal binary (i.e marking the even/odd chromosomes for plot)
  tile_modal = tile_modal %>% mutate(evenChr = case_when(seqnames == 'chr2' ~ TRUE,
                                                         seqnames == 'chr4' ~ TRUE,
                                                         seqnames == 'chr6' ~ TRUE,
                                                         seqnames == 'chr8' ~ TRUE,
                                                         seqnames == 'chr10' ~ TRUE,
                                                         seqnames == 'chr12' ~ TRUE,
                                                         seqnames == 'chr14' ~ TRUE,
                                                         seqnames == 'chr16' ~ TRUE,
                                                         seqnames == 'chr18' ~ TRUE,
                                                         seqnames == 'chr20' ~ TRUE,
                                                         seqnames == 'chr22' ~ TRUE,
                                                         TRUE ~ FALSE))
  
  
#FRAGILE SITES
#adding fragile sites to the previous plot:
  fragile = read_tsv("/home/noshadh/Codes/Germline_variation_detection/Fragile_sites", col_names = FALSE)
  fragile = fragile %>% select(X1,X2,X3,X4) #only need these values
  colnames(fragile) = c('chr','start','end','site_name')
  #map to tiles
  fragile = fragile %>% mutate(tile_start = 0) %>% mutate(tile_end = 0)
  for (i in 1:dim(fragile)[1]){
    chr = fragile$chr[i] #select the chromosome that fragile site is on
    chr_tiles = tile_modal %>% filter(seqnames == chr) #select all the tiles on that fragile site
    
    start_tile_on_chr = round(fragile[i,]$start/10000)+1
    end_tile_on_chr = round(fragile[i,]$end/10000)
    
    if (fragile[i,]$end > max(chr_tiles$end)){ #there might be some fragile sites that would fall outside the chromosome targets boundaries that we have
      end_tile_on_chr = dim(chr_tiles)[1]
    }
    fragile[i,5] = chr_tiles[start_tile_on_chr,] %>% select(tile)
    fragile[i,6] = chr_tiles[end_tile_on_chr,] %>% select(tile)
    
    print(paste("Fragile site number",i,"has been processed..."))
  }
  
  #convert the positions, to 
  s_positions = fragile %>% select(pos = tile_start) %>% mutate(type = "start")
  e_positions = fragile %>% select(pos = tile_end) %>% mutate(type = "end")
  positions = rbind(s_positions, e_positions) #didn't used it at the end
  
  p = p + geom_segment(aes(x=s_positions$pos,xend=e_positions$pos,y=1,yend=1)) #this p plot, has the boundaries for fragile sites, to be added to the final plot
  
  #plot the positions
   ggplot()+
    geom_histogram(na.rm = TRUE,data = (tile_modal %>% filter(modal_num < 8)),binwidth = 100,alpha = 0.6,aes(x = tile, color = as.factor(evenChr))) + 
    xlim(0,287509)+ theme_classic()+
     geom_hline(yintercept = 8)+scale_fill_manual(values=c("#e1e6cf", "dark green"))+
    geom_vline(xintercept = 24896,linetype = "dashed",size = 0.3)+
    geom_text(x = 24896, y = -4.5, label = "chr1", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 49116,linetype = "dashed",size = 0.3)+
     geom_text(x = 49116, y = -4.5, label = "chr2", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 68946,linetype = "dashed",size = 0.3)+
     geom_text(x = 68946, y = -4.5, label = "chr3", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 87968,linetype = "dashed",size = 0.3)+
     geom_text(x = 87968, y = -4.5, label = "chr4", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 106122,linetype = "dashed",size = 0.3)+
     geom_text(x = 106122, y = -4.5, label = "chr5", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 123203,linetype = "dashed",size = 0.3)+
     geom_text(x = 123203, y = -4.5, label = "chr6", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 139138,linetype = "dashed",size = 0.3)+
     geom_text(x = 139138, y = -4.5, label = "chr7", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 153652,linetype = "dashed",size = 0.3)+
     geom_text(x = 153652, y = -4.5, label = "chr8", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 167492,linetype = "dashed",size = 0.3)+
     geom_text(x = 167492, y = -4.5, label = "chr9", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 180872,linetype = "dashed",size = 0.3)+
     geom_text(x = 180872, y = -4.5, label = "chr10", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 194381,linetype = "dashed",size = 0.3)+
     geom_text(x = 194381, y = -4.5, label = "chr11", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 207709,linetype = "dashed",size = 0.3)+
     geom_text(x = 207709, y = -4.5, label = "chr12", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 219146,linetype = "dashed",size = 0.3)+
     geom_text(x = 219146, y = -4.5, label = "chr13", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 229851,linetype = "dashed",size = 0.3)+
     geom_text(x = 229851, y = -4.5, label = "chr14", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 240051,linetype = "dashed",size = 0.3)+
     geom_text(x = 240051, y = -4.5, label = "chr15", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 249085,linetype = "dashed",size = 0.3)+
     geom_text(x = 249085, y = -4.5, label = "chr16", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 257411,linetype = "dashed",size = 0.3)+
     geom_text(x = 257411, y = -4.5, label = "chr17", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 265449,linetype = "dashed",size = 0.3)+
     geom_text(x = 265449, y = -4.5, label = "chr18", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 271311,linetype = "dashed",size = 0.3)+
     geom_text(x = 271311, y = -4.5, label = "chr19", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 277756,linetype = "dashed",size = 0.3)+
     geom_text(x = 277756, y = -4.5, label = "chr20", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 282427,linetype = "dashed",size = 0.3)+
     geom_text(x = 282427, y = -4.5, label = "chr21", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
    geom_vline(xintercept = 287509,linetype = "dashed",size = 0.3)+
     geom_text(x = 287509, y = -4.5, label = "chr22", size = 4, angle = 90, vjust = -0.4, hjust= 0,aes(na.rm = TRUE))+
     geom_segment(aes(x=s_positions$pos,xend=e_positions$pos,y=0,yend=0))+ggtitle("3PC equal to zero") #geom segment add fragile sites
    #if add two lines below to figure, chrx and chry will be added
        #    +geom_vline(xintercept = 303114,linetype = "dashed",size = 0.1)+ geom_text(x = 303114, y = 0, label = "chrX", size = 4, angle = 90, vjust = -0.4, hjust= 0)
#    +geom_vline(xintercept = 308837,linetype = "dashed",size = 0.1)+ geom_text(x = 308837, y = 0, label = "chrY", size = 4, angle = 90, vjust = -0.4, hjust= 0)
   