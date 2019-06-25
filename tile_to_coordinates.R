tile_to_coordinates = function(tile,file){ #file is a one of the .RDS files of patients, don't matter which one
  range = as.data.frame(file$tile[tile]) %>% select(seqnames,start,end,blacklist)
  return(range)
}


coordinates_to_tile = function(chr,start_coord,end_coord,file){
  range = as.data.frame(file$tile)
  range = range %>% mutate(tile = row_number())
  start_tile = (range %>% filter(start < start_coord & end > start_coord & seqnames == as.character(chr)))$tile
  end_tile = (range %>% filter(start < end_coord & end > end_coord & seqnames == as.character(chr)))$tile
  output = as.data.frame(matrix(nrow = 1,ncol = 2))
  output[,1] = start_tile
  output[,2] = end_tile
  colnames(output) = c('start_tile','end_tile')
  return(output)
}

