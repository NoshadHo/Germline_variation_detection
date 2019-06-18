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
  step_size = 10000
  variance_list = foreach(i = 0:(round(dim(data)[2]/step_size)-1)) %dopar% {
    subdata = data[,(i*step_size+1):(min((i+1)*step_size,dim(data)[2]))]
    subdata %>% summarise_all(sd)
  }
  })
  variance_df = as.data.frame(t(do.call(cbind,variance_list)))
  return(variance_df)
}

variance_df = variance_sex(as.data.frame(t(sex_removed_tile_cov_gc)),"m")
variance_df_0 = variance_df
ggplot()+geom_point(aes(x = 1:dim(variance_df_1)[1],y = variance_df$V1),size = 0.4)+theme_minimal()+geom_vline(xintercept = 287509,linetype = "dashed",size = 0.3)+
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
  ggplot()+geom_point(aes(x = 1:TILES_NUMBER,y = (variance_diff_0_1 %>% arrange(desc(V1)))[1:TILES_NUMBER,]),size = 0.7)+
    theme_minimal()+
    labs(title = "0-1 difference",x = "Genome tiles", y = "Variance")
  
  TILES_NUMBER = 120
  end = dim(variance_diff_0_1)[1]
  ggplot()+geom_point(aes(x = end:(end-TILES_NUMBER),y = (variance_diff_0_1 %>% arrange(desc(V1)))[end:(end-TILES_NUMBER),]),size = 0.7)+
    theme_minimal()+
    labs(title = "0-1 difference",x = "Genome tiles", y = "Variance")
    