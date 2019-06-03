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
    set.seed(1024)

  #Data load:
  #we have saved it in tile_case_comparison.R file, we just load it from there
    tile_case = read_tsv('/home/noshadh/Codes/Germline_variation_detection/Tile_case_hl5.tsv',col_names = TRUE)
  #Add tile number as the first column
    tile_case = as.data.frame(tile_case) %>% add_column(tile = 0, .before = TRUE) #add index
    tile_case = tile_case %>% mutate(tile = row_number())
  
  #load the selected tiles from tile_Case_comparison.R file
    vol_tiles = read_tsv('/home/noshadh/Codes/Germline_variation_detection/Selected_volatile_tiles.tsv',col_names = TRUE)
  #join it with tile case to capture ncov:
    vol_tiles = vol_tiles %>% select(tile) #only select the ids
    vol_tiles = vol_tiles %>% left_join(tile_case,by = "tile")
    vol_tiles = vol_tiles %>% select(-tile)
#CAPTURING THE CORRELATION COEFFICIENT----------------------------------------------------------------------------------------------------------
  #---FIrst we will try our approach on the selected tiles from tile_Case_comparison.R file---
  #predefinning the data.frame
    tile_tile_corr_selected = data.frame(matrix(nrow = dim(vol_tiles)[1],ncol = dim(vol_tiles)[1]))
  #calculating the correlation (this function calculate between columns, so we transpose the data.frame)
    tile_tile_corr_selected = cor(t(vol_tiles))
    tile_tile_corr_selected = as.data.frame(tile_tile_corr_selected)
  #write the file
    write_tsv(tile_tile_corr_selected, "/home/noshadh/Codes/Germline_variation_detection/selected_tiles_pairwise_correlation.tsv")

#10MB BOX-PLOTS ANALYSIS-------------------------------------------------------------------------------------------------------------------------
#IN THIS PART, WE MEASURE EACH SELECTED TILE CORRELATION WITH:
#GROUP 1: ALL THE TILES IN THE DISTANCE OF <10MB
#GROUP 2: ALL THE TILES IN THE DISTANCE OF >10MB AND <20MB
#...
#DO THIS FOR ALL SELECTED TILES
#KEEP THE RESULTS IN A DATA.FRAME
#annotating where tiles belong (which box group)
#every 10mb (1000 tile) go into one group