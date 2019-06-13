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
    set.seed(1024)
  
  #Data load: #we will use data loaded in tiles_Correlation_Analysis
    base::load('/home/noshadh/Codes/Germline_variation_detection/K-means_Multimodal_fitting_data_gcNormlized.RData')

#TILES_PCA PLOT---------------------------------------------------------------------------------------------------------------------------------
  #PLOTING TILES PCA USING GC-NORMALIZED COVERAGE DATA FOR NORMAL PATIENTS
    autoplot(prcomp((tile_cov_gc_normalized)), loadings = FALSE)+theme_minimal()
  #analysi the PC information
    pca_data = prcomp((tile_cov_gc_normalized))
    pc_variance = as.data.frame((pca_data$sdev)^2/sum((pca_data$sdev)^2))
    pc_cummulative_variation = cumsum(pc_variations)
    pc_variance = pc_variance %>% mutate(cumsum = pc_cummulative_variation$pc_variance)
    colnames(pc_variance) = c("variance", "cum_variance")
    rm(pc_cummulative_variation)
  #plot these values:
    pc_variance %>% ggplot()+
      geom_point(aes(y = variance, x = 1:dim(pc_variance)[1]))+theme_minimal()
  #By looking at these values, I decided to remove three of PC's
    num_PC = 3
    
#REMOVING THE COMPONENTS FROM NORMAL DATA-------------------------------------------------------------------------------------------------------
  #first method: find eigen vslues, rotate data to that space, make the value zero, rotate back
    #eigen_vectors = eigen(cov(tile_cov_gc_normalized))$vectors
    eigen_vectors = prcomp(tile_cov_gc_normalized, center = TRUE)$rotation
    rotated_tile_cov_gc_normalized = as.matrix(tile_cov_gc_normalized) %*% eigen_vectors
    rotated_tile_cov_gc_normalized[, 1:num_PC] <- 0
    purified_tile_cov_gc_normalized <- rotated_tile_cov_gc_normalized %*% t(eigen_vectors) # transpose of orthogonal matrix = inverse
    
  #QC: turn the space into pc subspace and plot
    QC_matrix = as.matrix(purified_tile_cov_gc_normalized) %*% eigen_vectors
    ggplot()+geom_point(aes(x = QC_matrix[,1], y = QC_matrix[,2]))
  #second method: map unto the subplane spaned by those PC's, map each vector onto that subplane, remove the values from the vector