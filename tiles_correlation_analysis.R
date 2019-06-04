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
#BOXPLOTS FOR  ONE TILE ------------------------------------------------------------------------------------------------------------------------------------------
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

#HEATMAP PLOTS FOR  ALL TILES ON EACH CHROMOSOME------------------------------------------------------------------------------------------------------------------------------------------
  #CHROMOSOME 17 IS FROM TILE 249085 TO 257411
    data_file = tile_case[249085:257411,] #to work with positiive/negetive file more easily
    data_file_corr = cor(t(data_file),method = "spearman") #find the correlation
    data_file_corr = abs(data_file_corr)
    colnames(data_file_corr) = 249085:257411
    rownames(data_file_corr) = 249085:257411
    
    pheatmap::pheatmap(as.matrix(data_file_corr),cluster_rows = FALSE,cluster_cols = FALSE,
                       main = "CHR17 tile correlation")
    

#BOXPLOT FOR ALL TILES ON EACH CHROMOSOME-------------------------------------------------------------------------------------------------------------------------------------------------
  #find the correlation matrix on the chromosome
    tile_num = 257411-249085+1
    max_group_nums = floor(tile_num/100) + 1
    data_file = tile_case[249085:257411,] #to work with positiive/negetive file more easily
    data_file_corr = cor(t(data_file),method = "spearman") #find the correlation
    data_file_corr = as.data.frame(abs(data_file_corr))
    colnames(data_file_corr) = 1:tile_num
    rownames(data_file_corr) = 1:tile_num
    
  #Group points for each 1MB (100tiles)
    grouping_values = as.data.frame(matrix(nrow = max_group_nums, ncol = 0)) #249 groups in columns for each tile in row
    grouping_values = grouping_values %>% dplyr::mutate(group = 1:max_group_nums)
    for (tile in 1:tile_num){
      row_values = data_file_corr %>% slice(tile)
      row_values = as.data.frame(t(row_values))
      row_values = row_values %>% mutate(group = floor(abs((row_number()-tile)/100))+1)
      #now take a mean over values in each group (we represent each group in a tile as a point)
      colnames(row_values) = c("corr", 'group')
      row_values = row_values %>% group_by(group) %>% summarise(mean(corr))
      #grouping_values[row,] = t(row_values[2])

      grouping_values = grouping_values %>% left_join(row_values,by = "group")
      print(paste("Tile processed:",tile))
    }
    grouping_values_backup = grouping_values
    grouping_values = grouping_values %>% select(-group)
    colnames(grouping_values)[1:tile_num] = 1:tile_num
    
  #plot the boxplot
    #making the data ready
    group_numbers_to_plot = max_group_nums
    df = as.data.frame(t(grouping_values[1:group_numbers_to_plot,])) #we only look at 10 groups
    colnames(df) = 1:group_numbers_to_plot
    df = gather(df) #it needs to be in this form for the plot
    df[,1] = as.integer(df[,1]) #to show 2 after 1 and not 10 (don't treat it like a string)
    df = df %>% arrange(key)
    
    p = gather(df) %>% ggplot(aes(x = key, y = value, fill = key)) + ggtitle("Chr17 distance boxplots")+xlim(0,group_numbers_to_plot+1)+
      geom_boxplot(na.rm = TRUE)+theme_minimal()
    #+theme(axis.text.y = element_blank(),axis.ticks = element_blank())
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
    
