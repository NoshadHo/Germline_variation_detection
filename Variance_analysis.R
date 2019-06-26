################################################################################################################################################
## Author: Noshad Hosseini                                                                                                                    ##
## Date: 06-12-2019                                                                                                                           ##
## Description: We are going to first plot the pca of our tile-coverage data, then we will find it's principle components, using linear       ##
##              Algebra to remove those from data (removing the variation, which is caused by the noise, it is a denoising thing)             ##
################################################################################################################################################
#PREPRATION-------------------------------------------------------------------------------------------------------------------------------------
#Libraries:
library(GenomicRanges)
library(cnvex) #this is a local library, might not work on other machines
library(tidyverse)
library(rlist)
library(gplots)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(mclust)
library(gridExtra)
library(corrr)
library(plotly)
library(rlang)
library(limma) #bioconductor lib
library(gridExtra)
library(ggfortify)
library(mclust)
library(foreach)
library(doParallel)
library(Rcpp)
library(fpc)
library(dbscan)
set.seed(1024)
numCores = detectCores()
registerDoParallel(numCores-1)
set.seed(1024)

#Data load: #we will use data loaded in tiles_Correlation_Analysis
base::load('/home/noshadh/Codes/Germline_variation_detection/K-means_Multimodal_fitting_data_gcNormlized.RData')
sex = read_rds(path = "/home/noshadh/Codes/Germline_variation_detection/sex.rds")
data = as.data.frame(t(tile_cov_gc_normalized))

#CACULATE VARIANCE FOR EACH TILE, FOR EACH SEX--------------------------------------------------------------------------------------------------
variance_sex = function(data, sex){ #data here is tile_coverage matrix, with 110 rows
  time1 = system.time({ #we break data set to data sets with the size equal to step_size, then run them parallel.
    #processing big data self will be too time consuming.
  STEP_SIZE = 5000
  variance_list = foreach(i = 0:(round(dim(data)[2]/STEP_SIZE)-1)) %dopar% {
    subdata = data[,(i*STEP_SIZE+1):(min((i+1)*STEP_SIZE,dim(data)[2]))]
    subdata %>% summarise_all(sd)
  }
  })
  variance_df = as.data.frame(t(do.call(cbind,variance_list)))
  print(time1)
  return(variance_df)
}

variance_df = variance_sex(as.data.frame(t(purified_tile_cov_gc_normalized)),"m")
variance_df = variance_df %>% dplyr::mutate(tile = blacklist_removed_tile_list$tile)
variance_df_bl_6 = variance_df
ggplot()+geom_point(aes(x = variance_df$tile,y = variance_df$V1),size = 0.4)+theme_minimal()+geom_vline(xintercept = 287509,linetype = "dashed",size = 0.3)+
  geom_vline(xintercept = 303114,linetype = "dashed",size = 0.3)+ylim(0,2)+
  labs(title = paste("PC removed",sc_num),x = "Genome tiles", y = "Variance")

#LOOK AT THE DIFFERENCE OF VARIANCE EACH STEP OF PC REMOVAL--------------------------------------------------------------------------------------
  variance_diff_0_1 = variance_df_0 - variance_df_1
  variance_diff_1_2 = variance_df_1 - variance_df_2
  variance_diff_2_3 = variance_df_2 - variance_df_3
  variance_diff_3_4 = variance_df_3 - variance_df_4
  variance_diff_4_5 = variance_df_4 - variance_df_5
  variance_diff_5_6 = variance_df_5 - variance_df_6

  
  TILES_NUMBER = 50
  ggplot()+geom_point(aes(x = 1:TILES_NUMBER,y = (variance_diff_bl_0_1 %>% arrange(desc(V1)))[1:TILES_NUMBER,]),size = 0.7)+
    theme_minimal()+
    labs(title = "0-1 difference",x = "Genome tiles", y = "Variance")
  
  TILES_NUMBER = 120
  end = dim(variance_diff_0_1)[1]
  ggplot()+geom_point(aes(x = end:(end-TILES_NUMBER),y = (variance_diff_0_1 %>% arrange(desc(V1)))[end:(end-TILES_NUMBER),]),size = 0.7)+
    theme_minimal()+
    labs(title = "0-1 difference",x = "Genome tiles", y = "Variance")
  #results of this part can be find in "Principal Component noise reduction results" google sheet
  
  # for blacklist
  variance_diff_bl_0_1$V1 = variance_df_bl_0$V1 - variance_df_bl_1$V1
  variance_diff_bl_0_1$tile = variance_df_bl_0$tile
  variance_diff_bl_1_2$V1 = variance_df_bl_1$V1 - variance_df_bl_2$V1
  variance_diff_bl_1_2$tile = variance_df_bl_1$tile
  variance_diff_bl_2_3$V1 = variance_df_bl_2$V1 - variance_df_bl_3$V1
  variance_diff_bl_2_3$tile = variance_df_bl_2$tile
  variance_diff_bl_3_4$V1 = variance_df_bl_3$V1 - variance_df_bl_4$V1
  variance_diff_bl_3_4$tile = variance_df_bl_3$tile
  variance_diff_bl_4_5$V1 = variance_df_bl_4$V1 - variance_df_bl_5$V1
  variance_diff_bl_4_5$tile = variance_df_bl_4$tile
  variance_diff_bl_5_6$V1 = variance_df_bl_5$V1 - variance_df_bl_6$V1
  variance_diff_bl_5_6$tile = variance_df_bl_5$tile
  TILES_NUMBER = 50
  ggplot()+geom_point(aes(x = 1:TILES_NUMBER,y = (variance_diff_bl_0_1 %>% arrange(desc(V1)))[1:TILES_NUMBER,1]),size = 0.7)+
    theme_minimal()+
    labs(title = "0-1 difference",x = "Genome tiles", y = "Variance")
  
  TILES_NUMBER = 120
  end = dim(variance_diff_0_1)[1]
  ggplot()+geom_point(aes(x = end:(end-TILES_NUMBER),y = (variance_diff_4_5 %>% arrange(desc(V1)))[end:(end-TILES_NUMBER),1]),size = 0.7)+
    theme_minimal()+
    labs(title = "0-1 difference",x = "Genome tiles", y = "Variance")
  #
  
#FIND THE TILES, AND THEN ANNOTATE THEM, SEE IF WE CAN MAKE ANY SENSE OF THEM--------------------------------------------------------------------
  #set the values selected in the last part:
  TOP_0_1 = 8
  END_0_1 = 5
  TOP_1_2 = 9
  TOP_2_3 = 12
  TOP_3_4 = 4
  TOP_4_5 = 3
  TOP_5_6 = 3
  
  #add tile id as a column to datasets (FOR BLACKLIST REMOVED, DON'T DO THIS STEP)
  variance_diff_0_1 = variance_diff_0_1 %>% mutate(row_number())
  variance_diff_1_2 = variance_diff_1_2 %>% mutate(row_number())
  variance_diff_2_3 = variance_diff_2_3 %>% mutate(row_number())
  variance_diff_3_4 = variance_diff_3_4 %>% mutate(row_number())
  variance_diff_4_5 = variance_diff_4_5 %>% mutate(row_number())
  variance_diff_5_6 = variance_diff_5_6 %>% mutate(row_number())
  #select the number of tiles:
  variant_tiles_0_1 = variance_diff_0_1 %>% arrange(desc(V1)) %>% slice(1:TOP_0_1)
  variant_tiles_1_2 = variance_diff_1_2 %>% arrange(desc(V1)) %>% slice(1:TOP_1_2)
  variant_tiles_2_3 = variance_diff_2_3 %>% arrange(desc(V1)) %>% slice(1:TOP_2_3)
  variant_tiles_3_4 = variance_diff_3_4 %>% arrange(desc(V1)) %>% slice(1:TOP_3_4)
  variant_tiles_4_5 = variance_diff_4_5 %>% arrange(desc(V1)) %>% slice(1:TOP_4_5)
  variant_tiles_5_6 = variance_diff_5_6 %>% arrange(desc(V1)) %>% slice(1:TOP_5_6)
  variant_tiles_0_1_end = variance_diff_0_1 %>% arrange(V1) %>% slice(1:END_0_1)
  colnames(variant_tiles_0_1) = c("var_diff", "tile")
  colnames(variant_tiles_1_2) = c("var_diff", "tile")
  colnames(variant_tiles_2_3) = c("var_diff", "tile")
  colnames(variant_tiles_3_4) = c("var_diff", "tile")
  colnames(variant_tiles_4_5) = c("var_diff", "tile")
  colnames(variant_tiles_5_6) = c("var_diff", "tile")
  colnames(variant_tiles_0_1_end) = c("var_diff", "tile")

  #FOR BALCKLIST REMOVED DATAFRAME
  #select the number of tiles:
  variant_tiles_bl_0_1 = variance_diff_bl_0_1 %>% arrange(desc(V1)) %>% slice(1:TOP_0_1)
  variant_tiles_bl_1_2 = variance_diff_bl_1_2 %>% arrange(desc(V1)) %>% slice(1:TOP_1_2)
  variant_tiles_bl_2_3 = variance_diff_bl_2_3 %>% arrange(desc(V1)) %>% slice(1:TOP_2_3)
  variant_tiles_bl_3_4 = variance_diff_bl_3_4 %>% arrange(desc(V1)) %>% slice(1:TOP_3_4)
  variant_tiles_bl_4_5 = variance_diff_bl_4_5 %>% arrange(desc(V1)) %>% slice(1:TOP_4_5)
  variant_tiles_bl_5_6 = variance_diff_bl_5_6 %>% arrange(desc(V1)) %>% slice(1:TOP_5_6)
  variant_tiles_bl_0_1_end = variance_diff_bl_0_1 %>% arrange(V1) %>% slice(1:END_0_1)
  
    #this is how I looked at the reagions selected:
  #variant_tiles_4_5 %>% arrange(tile)
  
  t1 = intersect(variant_tiles_0_1$tile, variant_tiles_1_2$tile)
  t2 = intersect(t1, variant_tiles_2_3$tile)
  t3 = intersect(t2, variant_tiles_3_4$tile)
  t4 = intersect(t3, variant_tiles_4_5$tile)
  t5 = intersect(t4, variant_tiles_5_6$tile)
  
  
