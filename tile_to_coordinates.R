tile_to_coordinates = function(tile,file){ #file is a one of the .RDS files of patients, don't matter which one
  range = as.data.frame(file$tile[tile]) %>% select(seqnames,start,end,blacklist)
  return(range)
}


coordinates_to_tile = function(chr,start,end,file){
    
}
