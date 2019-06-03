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
    set.seed(1024)

  #Data load:
    #base::load('/home/noshadh/Codes/Germline_variation_detection/K-means_Multimodal_fitting_data.RData')
    #tile_case = t(tile_case)
  #we have saved it in tile_case_comparison.R file, we just load it from there
    tile_case = read_tsv('/home/noshadh/Codes/Germline_variation_detection/Tile_case_hl5.tsv',col_names = TRUE)
  #Add tile number as the first column
    tile_case = as.data.frame(tile_case) %>% add_column(tile = 0, .before = TRUE) #add index
    tile_case = tile_case %>% mutate(tile = row_number())
  
  #load the selected tiles from tile_Case_comparison.R file
    vol_tiles = read_tsv('/home/noshadh/Codes/Germline_variation_detection/Selected_volatile_tiles.tsv',col_names = TRUE)
    #vol_tiles = significant_tiles
    tile_order =vol_tiles$tile
  #join it with tile case to capture ncov:
    vol_tiles = vol_tiles %>% select(tile) #only select the ids
    vol_tiles = vol_tiles %>% left_join(tile_case,by = "tile")
    rownames(vol_tiles) = vol_tiles$tile
    vol_tiles = vol_tiles %>% select(-tile)
#CAPTURING THE CORRELATION COEFFICIENT----------------------------------------------------------------------------------------------------------
  #---FIrst we will try our approach on the selected tiles from tile_Case_comparison.R file---
  #predefinning the data.frame
    tile_tile_corr_selected = (matrix(nrow = dim(vol_tiles)[1],ncol = dim(vol_tiles)[1]))
  #calculating the correlation (this function calculate between columns, so we transpose the data.frame)
    tile_tile_corr_selected = cor(t(vol_tiles))
    tile_tile_corr_selected = as.data.frame(tile_tile_corr_selected) #this step is not neccessary
    colnames(tile_tile_corr_selected) = tile_order
    rownames(tile_tile_corr_selected) = tile_order
    #make the correlations under 0 to it's absoloute value
    tile_tile_corr_selected_pos = abs(tile_tile_corr_selected)
#PLOTS FOR  ONE TILE ------------------------------------------------------------------------------------------------------------------------------------------
  #a simple plot for a arbitary tile, to look at it's correlation with other tiles
    data_file = tile_tile_corr_selected_pos #to work with positiive/negetive file more easily
    tile1_corr = data_file[1000,]
    tile_nums = 1:dim(data_file)[1]
    tile1_corr = rbind(tile1_corr,tile_nums)
    #group points for each 10mb
    tile1_corr = rbind(tile1_corr,tile_nums)
    
    group = 1
    for (i in 1:dim(data_file)[1]){
      if (i > group*1000){
        group = group + 1
      }
      tile1_corr[3,i] = (group)
      
      print(i)
    }
    tile1_corr[3,] = as.factor(tile1_corr[3,])
    
    df = as.data.frame(t(tile1_corr[,20000:30000]))
    colnames(df) = c('value','index','group')
    p = df %>% ggplot(aes(x = group, y = value, fill = group)) +
        geom_boxplot(na.rm = TRUE) + theme(
          axis.text.y = element_blank(),
          axis.ticks = element_blank())
    ggplotly(p)

#PLOTS FOR  ALL TILES ON CHROMOSOMES ------------------------------------------------------------------------------------------------------------------------------------------
  #CHROMOSOME 19 IS FROM TILE 265449 TO 271311
    data_file = tile_case[265449:271311,] #to work with positiive/negetive file more easily
    data_file_corr = cor(t(data_file),method = "spearman") #find the correlation
    data_file_corr = abs(data_file_corr)
    colnames(data_file_corr) = 265449:271311
    rownames(data_file_corr) = 265449:271311
    heatmap.2(data_file_corr[200:900,200:900],density.info="none", trace="none",Rowv = NULL, Colv = NULL)
    corrplot(data_file_corr[1:200,1:200],type = 'lower')
    pheatmap::pheatmap(as.matrix(data_file_corr),cluster_rows = FALSE,cluster_cols = FALSE,
                       main = "CHR19 tile correlation")
  #specify each tile groups (we devide tiles into groups of 100 on genome, 1Mb)
    data_file_corr = data_file %>% mutate(group = 0)
    group = 1
    for (i in 1:dim(data_file)[1]){
      if (tile_order[i] > group*2500){
        group = group + 1
      }
      #data_file$group[i] = group #changed it with the line below to speed up
      data_file[i,42736] = group 
      print(i)
    }
  #take mean over all the patients for each tile (we will have 1 value for each tile)
  #a simple plot for a arbitary tile, to look at it's correlation with other tiles
    
    tile1_corr = data_file[1000,]
    tile_nums = 1:dim(data_file)[1]
    tile1_corr = rbind(tile1_corr,tile_nums)
  #group points for each 10mb
    tile1_corr = rbind(tile1_corr,tile_nums)
    
    group = 1
    for (i in 1:dim(data_file)[1]){
      if (i > group*1000){
        group = group + 1
      }
      tile1_corr[3,i] = (group)
      
      print(i)
    }
    tile1_corr[3,] = as.factor(tile1_corr[3,])
    
    df = as.data.frame(t(tile1_corr[,10000:20000]))
    colnames(df) = c('value','index','group')
    p = df %>% ggplot(aes(x = group, y = value, fill = group)) +
      geom_boxplot(na.rm = TRUE) + theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
    ggplotly(p)
    
    #10MB BOX-PLOTS ANALYSIS------------------------------------------------------------------------------------------------------------------------
#IN THIS PART, WE MEASURE EACH SELECTED TILE CORRELATION WITH:
#GROUP 1: ALL THE TILES IN THE DISTANCE OF <10MB
#GROUP 2: ALL THE TILES IN THE DISTANCE OF >10MB AND <20MB
#...
#DO THIS FOR ALL SELECTED TILES
#KEEP THE RESULTS IN A DATA.FRAME
#annotating where tiles belong (which box group)
#every 10mb (1000 tile) go into one group
    
