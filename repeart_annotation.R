repeat_file = read_tsv("~/Codes/Germline_variation_detection/genome_repeat_5.bed",col_names = TRUE)

#generate a file for bed_to_tile function:
repeat_file_short = repeat_file %>% select(genoName,genoStart,genoEnd,repClass)
repeat_file_short = repeat_file_short %>% filter(nchar(genoName) == 4)
colnames(repeat_file_short) = c('chr','start','end','repeat_class')  



bed_to_tile_modified = function(bed_file,col_names = FALSE, RDSfile){
  if(col_names == FALSE){
    colnames(bed_file) = c('chr','start','end','repeat_class')  
  }
  tiles = list()
  
  PTIME = system.time({
    tiles = foreach(i = 1:5000) %dopar% {
      return(c((coordinates_to_tile(bed_file$chr[i], bed_file$start[i],bed_file$end[i],file))))
  }
  })
  tiles = as.data.frame(do.call(rbind,tiles))
  colnames(tiles) = c("start_tile","end_tiles")
  tiles = tiles %>% mutate(repeat_class = repeat_file_short$repeat_class[1:499])
  #if a region contains more than one tile, make it to a couple of entry, each containing one tile
  more_than_one_tile = tiles %>% filter(start_tile != end_tiles)
  tiles_single = tiles %>% anti_join(more_than_one_tile)
  
  tile_list = (tiles_single$start_tile)
  tile_anno = tiles_single$repeat_class
  for (i in 1:dim(more_than_one_tile)[1]){
    tile_list = c(tile_list,(more_than_one_tile$start_tile[i]):(more_than_one_tile$end_tiles[i]))
    for (i in (more_than_one_tile$start_tile[i]):(more_than_one_tile$end_tiles[i])) {tile_anno = c(tile_anno, more_than_one_tile$repeat_class[i])}
    print(i)
  }
  
  return(tile_list)
}

repeat_tilies = bed_to_tile_modified(repeat_file_short,RDSfile = file)
