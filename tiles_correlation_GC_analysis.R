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
    library(limma) #bioconductor lib
    library(gridExtra)
    set.seed(1024)

  #Data load: #we will use data loaded in tiles_Correlation_Analysis
    base::load('/home/noshadh/Codes/Germline_variation_detection/K-means_Multimodal_fitting_data.RData')
    tile_case = as.data.frame(t(tile_case))
    tile_coverage = tile_case #it makes more sense to select this name for our variable
    rm(tile_case)

  #make a tile-gc matrix (gc content is the same for all the patients)
    tile_gc = as.data.frame(file$tile) %>% select(gc)             #each row_number shows a tile and the value is gc value

#ANALYSIS FOR ONE PATIENT------------------------------------------------------------------------------------------------------------------------
  #make a data frame with tile num-coverage-gc_contetnt
    patient_num = 23
    tile_cov_gc = tile_coverage %>% select(cov = !!patient_num)
    tile_cov_gc = tile_cov_gc %>% mutate(gc = tile_gc$gc)
    
  #look at the gc-cov plot:
    tile_cov_gc %>% ggplot(aes(x = gc,y = cov))+
      geom_point(size = 0.4)
  #look at the gc-cov for cov < 1 (less than 1000 points are above 1)
    p_raw = (tile_cov_gc %>% filter(cov < 1)) %>% ggplot(aes(x = gc,y = cov))+
      geom_point(size = 0.4)+ggtitle(paste("raw plot for patient", patient_num))

# GC NORMALIZATION-------------------------------------------------------------------------------------------------------------------------------
  #use leoss fit for normalization
    tile_cov_gc_normalized_227 = tile_coverage
    for (patient in 1:dim(tile_coverage)[2]){
      gc.residuals = limma::loessFit(y = tile_coverage[[patient]], x=tile_gc$gc)$residuals
      cov.offset = lm(gc.residuals~tile_coverage[[patient]])$coefficients[1]
      gc.residuals = gc.residuals - cov.offset
      tile_cov_gc_normalized_227[patient] = gc.residuals
      print(paste("patient is proccessed:",patient))
    }
    save.image("~/Codes/Germline_variation_detection/K-means_Multimodal_fitting_data_gcNormlized.RData")
    #Main results here is tile_cov_gc_normalized. the rest is for QC.
  #Quality control, ploting the same plot as last part, with new data  
    tile_cov_gc = tile_cov_gc_normalized_227 %>% select(cov = !!patient_num)
    tile_cov_gc = tile_cov_gc %>% mutate(gc = tile_gc$gc)
    #look at the gc-cov for cov < 1 (less than 1000 points are above 1)
    p_norm = (tile_cov_gc %>% filter(cov < 1)) %>% ggplot(aes(x = gc,y = cov))+
      geom_point(size = 0.4)+ggtitle(paste("normalized plot for patient", patient_num))
    
    grid.arrange(p_raw,p_norm)
    #SEEMS GOOD :) (Our objective here was to normalize coverage so regions with low GC and very high GC would have same coverage as middle of plot)
    