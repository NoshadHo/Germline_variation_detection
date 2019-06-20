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
    autoplot(prcomp(t(sex_removed_tile_cov_gc_blacklist_newMask)), loadings = FALSE)+theme_minimal()
  #analysi the PC information
    pca_data = prcomp(t(tile_cov_gc_normalized))
    pc_variance = as.data.frame((pca_data$sdev)^2/sum((pca_data$sdev)^2))
    pc_cummulative_variation = cumsum(pc_variance)
    pc_variance = pc_variance %>% mutate(cumsum = pc_cummulative_variation$pc_variance)
    colnames(pc_variance) = c("variance", "cum_variance")
    rm(pc_cummulative_variation)
  #plot these values:
    pc_variance %>% ggplot()+
      geom_point(aes(y = variance, x = 1:dim(pc_variance)[1]))+theme_minimal()
  #By looking at these values, I decided to remove three of PC's
    num_PC = 6
    
#REMOVING THE COMPONENTS FROM NORMAL DATA-------------------------------------------------------------------------------------------------------
  #first method: find eigen vslues, rotate data to that space, make the value zero, rotate back
    #eigen_vectors = eigen(cov(tile_cov_gc_normalized))$vectors
    eigen_vectors = prcomp((sex_removed_tile_cov_gc_blacklist_newMask), center = TRUE)$rotation
    rotated_tile_cov_gc_normalized = as.matrix((sex_removed_tile_cov_gc_blacklist_newMask)) %*% eigen_vectors
    rotated_tile_cov_gc_normalized[, 1:num_PC] = 0
    purified_tile_cov_gc_normalized = rotated_tile_cov_gc_normalized %*% t(eigen_vectors) # transpose of orthogonal matrix = inverse
    
  #(OPTIONAL) Analysing the effect of principle reduction
    sum(purified_tile_cov_gc_normalized)
    #analysing column changes from pc1:3 removed compare to pc1:4 removed (expect to be pseudo uniformily)
    col_sum3 = as.data.frame(purified_tile_cov_gc_normalized) %>% summarise_all(funs(sum))  #sum of columns #change num_pc to 4 and try again for col_sum3
    col_changes34 = abs(col_sum3) - abs(col_sum4)                                           #sum of changes from removing pc4 from pc3removed case
    ggplot() + geom_point(aes(x = 1:dim(col_changes34)[2],y = t(col_changes34)))+theme_minimal()
    #analysing row changes from pc1:3 removed compare to pc1:4 removed (expect to be non-uniformly)
    row_sum3 = rowSums(x = purified_tile_cov_gc_normalized, dims = 1)
    row_changes34 = abs(row_sum2) - abs(row_sum3) 
    ggplot() + geom_point(aes(x = 1:length(row_changes34),y = (row_changes34)))+theme_minimal()+ggtitle("2 to 3")
  #QC: turn the space into pc subspace and plot
    QC_matrix = as.matrix(t(sex_removed_tile_cov_gc_blacklist)) %*% eigen_vectors #rotate it to the pca components
    ggplot()+geom_point(aes(x = QC_matrix[,1], y = QC_matrix[,5]))
  #second method: map unto the subplane spaned by those PC's, map each vector onto that subplane, remove the values from the vector
    
    
####SVD----------------------------------------------------------------------------------------------------------------------------------------
    #remove sex chromosomes:
    sex_removed_tile_cov_gc = tile_cov_gc_normalized[1:287509,]
    sex_removed_tile_cov_gc = sex_removed_tile_cov_gc %>% mutate(blacklist = blacklist$blacklist[1:287509])
    sex_removed_tile_cov_gc = sex_removed_tile_cov_gc %>% mutate(tile = row_number()) #to keep track of what tiles will remain
    sex_removed_tile_cov_gc_blacklist = sex_removed_tile_cov_gc %>% filter(blacklist  == 0)
    sex_removed_tile_cov_gc_blacklist = sex_removed_tile_cov_gc_blacklist %>% select(-blacklist)
    blacklist_removed_tile_list = sex_removed_tile_cov_gc_blacklist %>% select(tile)
    blacklist_removed_tile_list = blacklist_removed_tile_list %>% mutate(new_row = row_number())
    sex_removed_tile_cov_gc_blacklist = sex_removed_tile_cov_gc_blacklist %>% select(-tile)
    
    #ADDING NEW BLACKLIST MASK (THESE ARE THE REGIONS SELECTED BY VARIANCE ANALYSIS, THE REGIONS PC'S EFFECT THE MOST)
      sex_removed_tile_cov_gc_blacklist_newMask = (sex_removed_tile_cov_gc_blacklist %>% mutate(tile = blacklist_removed_tile_list$tile)) #first add the ORIGINAL ROW's to the matrix
      #now remove the tiles(rows) we want to be removed
      NEW_MASK = c(14320:14327,278577:278603,28188,96148,267422)  #it has the tile numbers
      sex_removed_tile_cov_gc_blacklist_newMask = sex_removed_tile_cov_gc_blacklist_newMask %>% filter(!(tile %in% NEW_MASK))
      blacklist_removed_tile_list_newMask = sex_removed_tile_cov_gc_blacklist_newMask %>% select(tile)
      blacklist_removed_tile_list_newMask = blacklist_removed_tile_list_newMask %>% mutate(new_row = row_number())
      sex_removed_tile_cov_gc_blacklist_newMask = sex_removed_tile_cov_gc_blacklist_newMask %>% select(-tile)
    #our matrix should have 110 rows (we want each point to be a patient and not a tile)
    svd = svd(sex_removed_tile_cov_gc_blacklist)
    #scree plot
    as.data.frame(svd$d) %>% ggplot()+
      geom_point(aes(y = svd$d, x = 1:length(svd$d)))+theme_minimal()+geom_line(aes(y = svd$d, x = 1:length(svd$d)))
    #we choose the number of sc to be deleted
    sc_num = 8
    
    svd$d[1:sc_num] = 0
    svd$d = diag(svd$d)
    purified_tile_cov_gc_normalized = svd$u %*% tcrossprod(svd$d,svd$v)
    
#look at the distribution of data before and after normalization
    plotDist = function(data, sample,start = 1,end = 308837){ #it should have 110 row and ...
      temp = as.data.frame((data[sample,])) %>% mutate(gr = 0) #we use this to just show the selected region on plot
      temp = temp %>% mutate(bl = if_else(file$tile$blacklist > 0,'blacklist','normal')) #bl is to show blacklisted regions
      temp$bl[257422:257422] = 'selected' #this shows the tile we are looking at
      row = as.data.frame((data[sample,start:end]))
      colnames(row) = "val"
      
      #row =row %>% filter(val < 1 & val > -1)
      ggplot()+geom_point(aes(x = 1:dim(row)[1],y = row$val,color = as.factor(temp[start:end,]$bl)),size = 0.8)+theme_linedraw()
      #+ylim(-10,45)
        #geom_point(aes(x = 1:dim(row)[1],y = row$val),size = 0.3)+theme_minimal()+ylim(0,10)
        #+geom_vline(xintercept = 287509,linetype = "dashed",size = 0.3)+
        #geom_vline(xintercept = 303114,linetype = "dashed",size = 0.3)+ylim(0,100)
    }
    
    plotDistPure = function(data, sample,start = 1,end = 308837){ #it should have 110 row and ...
      #temp = as.data.frame((data[sample,])) %>% mutate(gr = 0) #we use this to just show the selected region on plot
      #temp = temp %>% mutate(bl = if_else(file$tile$blacklist > 0,'blacklist','normal')) #bl is to show blacklisted regions
      #temp$bl[24952:24952] = 'selected' #this shows the tile we are looking at
      row = as.data.frame((data[sample,start:end]))
      colnames(row) = "val"
      
      #row =row %>% filter(val < 1 & val > -1)
      ggplot()+geom_point(aes(x = 1:dim(row)[1],y = row$val),size = 0.8)+theme_linedraw()+ylim(-10,45)
      #geom_point(aes(x = 1:dim(row)[1],y = row$val),size = 0.3)+theme_minimal()+ylim(0,10)
      #+geom_vline(xintercept = 287509,linetype = "dashed",size = 0.3)+
      #geom_vline(xintercept = 303114,linetype = "dashed",size = 0.3)+ylim(0,100)
    }
    
    plotDist(t(tile_cov_gc_normalized),65,1,12341)
    plotDistPure(t(purified_tile_cov_gc_normalized),sample,1,31804)
      
    
  #LOOK AT ONE CHROMOSOME LIKE THE THING WE HAVE ABOVE
    #chr1 long arm is from 12509 to 24896
    #chr21: short arm is from tile 277757 to 278857
    #how to find arms of chromosome tiles: 
      #temp = (as.data.frame(file$tile)) %>% mutate(tile = row_number())
      #(temp %>% group_by(arm) %>% slice(1))[34:37,]
    
    #analysis commands for a sample tile (14324 here)
    #variant_tiles_bl_8_9
    #file$tile[14324]
    #temp = (as.data.frame(file$tile)) %>% mutate(tile = row_number())
    #(temp %>% group_by(arm) %>% slice(1))[35:38,]
    #blacklist_removed_tile_list %>% filter(tile >= 12341) %>% slice(1)
    