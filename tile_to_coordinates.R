tile_to_coordinates = function(tile,RDSfile){ #RDSfile is a one of the .RDS files of patients, don't matter which one
  range = as.data.frame(file$tile[tile]) %>% select(seqnames,start,end,blacklist)
  return(range)
}


coordinates_to_tile = function(chr,start_coord,end_coord,RDSfile){
  range = as.data.frame(RDSfile$tile)
  range = range %>% mutate(tile = row_number())
  start_tile = (range %>% filter(start < start_coord & end > start_coord & seqnames == as.character(chr)))$tile
  end_tile = (range %>% filter(start < end_coord & end > end_coord & seqnames == as.character(chr)))$tile
  #output = as.data.frame(matrix(nrow = 1,ncol = 2))
  #output[,1] = start_tile
  #output[,2] = end_tile
  #colnames(output) = c('start_tile','end_tile')
  output = c(start_tile,end_tile)
  return(output)
}

bed_to_tile = function(file_address,col_names = FALSE, RDSfile){
  bed_file = read_tsv(file_address,col_names = col_names)
  if(col_names == FALSE){
    colnames(bed_file) = c('chr','start','end')  
  }
  tiles = list()
  for (i in 1:dim(bed_file)[1]){
    if (nchar(bed_file$chr[i]) == 4){
      tiles[i] = as.data.frame((coordinates_to_tile(bed_file$chr[i], bed_file$start[i],bed_file$end[i],file)))
      print(paste("Region has proccessed:",i))
    }
  }
  tiles = as.data.frame(do.call(rbind,tiles))
  colnames(tiles) = c("start_tile","end_tiles")
  return(tiles)
}
