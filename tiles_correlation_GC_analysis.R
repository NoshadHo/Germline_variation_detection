################################################################################################################################################
## Author: Noshad Hosseini                                                                                                                    ##
## Date: 06-03-2019                                                                                                                           ##
## Description: FInd the correlation between tiles in the genome, to detect artifacts. We basically expect to see correlation between tiles   ##
##              no more than 1MB further away then the tile.                                                                                  ##
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
    set.seed(1024)

  #Data load: #we will use data loaded in tiles_Correlation_Analysis
    base::load('/home/noshadh/Codes/Germline_variation_detection/K-means_Multimodal_fitting_data.RData')
    tile_case = as.data.frame(t(tile_case))

  #make a tile-gc matrix (gc content is the same for all the patients)
    tile_gc = as.data.frame(file$tile) %>% select(gc)             #each row_number shows a tile and the value is gc value

#ANALYSIS FOR ONE PATIENT------------------------------------------------------------------------------------------------------------------------
  #make a data frame with tile num-coverage-gc_contetnt
    patient_num = 100
    tile_cov_gc = tile_case %>% select(cov = !!patient_num)
    tile_cov_gc = tile_cov_gc %>% mutate(gc = tile_gc$gc)
      
  #look at the gc-cov plot:
    tile_cov_gc %>% ggplot(aes(x = gc,y = cov))+
      geom_point()


