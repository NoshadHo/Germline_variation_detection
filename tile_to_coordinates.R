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
  
  #if a region contains more than one tile, make it to a couple of entry, each containing one tile
  more_than_one_tile = tiles %>% filter(start_tile != end_tiles)
  tiles_single = tiles %>% anti_join(more_than_one_tile)
  
  tile_list = tiles_single$start_tile
  for (i in 1:dim(more_than_one_tile)[1]){
    tile_list = c(tile_list,(more_than_one_tile$start_tile[i]):(more_than_one_tile$end_tiles[i]))
    print(i)
  }
  
  return(tile_list)
}

file_address = "/home/noshadh/Codes/Germline_variation_detection/excluded_regions_cnvnator.bed"
ext_blacklist = bed_to_tile(file_address,RDSfile = file)

new_blacklist_tiles = c(1:90,260:280,12000:12341,1278:1340,1650:1697,24841:24871,
                        22861:22865,12341:14341,16144:16148,14341:15000,
                        257422:257422,33628:33640,25072,28183,28188,
                        34287:34367,34477:34505,34600:34660,35890:35965,36247:36262,36032:36052,35816,37887:37987,36032:36052,
                        56667:56695,58167:58207,
                        68657:68727,68807:68947,
                        73850:73947,73947:74147,
                        87910:87969,
                        68947:68957,
                        69855:69905,87969:88090,
                        92580:92849,
                        89710:89730,
                        90099:90119,92849:93020,
                        94820:95150,
                        105549:106123,106123:106173,112103:112303,
                        111953:112103,123103:123204,
                        109343:109423,122165:122200,123204:123254,129214:129564,
                        128830:129214,129710:129814,
                        130464:130564,139139:139250,143659:143800,153653:153703,157953:160500,171293:171473,
                        139800:140000,153553:153653,157503:157953,167453:167493,167493:167550,
                        140339:140439,147690:147730,
                        130654:130754,
                        133434:133504,
                        137550:137700,171473:171753,180873:180920,194300:194382,194382:194420,197932:198200,218875:218895,220867:221150,231752:232252,241510:241590,249040:249086,249086:249186,251596:251800,257412:257442,265300:265450,267870:268070,268070:268230,273850:274122,274122:274460,277757:278957,278957:279157,282428:283928,283928:284600,
                        180673:180873,189820:189900,207650:207710,209480:209580,229650:229852,239952:240052,241650:241730,251200:251596,252700:252750,258900:259262,265450:265520,265450:265520,277700:277757,282000:282250,287480:287510,
                        172050:172370,238220:238332,241852:241932,250900:250960,252850:252940,269075:269081,281150:281250,
                        232652:232752,242172:242332,253630:253760,
                        232820:232952,242902:242952,253405:253420,
                        233000:233152,243202:243732)
blacklist_new = blacklist %>% mutate(blacklist = case_when(blacklist == 1 ~ 1,
                                           tile %in% ext_blacklist ~ 1,
                                           TRUE ~ 0))
#make it sex exclusive
blacklist_new = blacklist_new[1:287509,]
sex_removed_tile_cov_gc = sex_removed_tile_cov_gc %>% mutate(blacklist = blacklist_new$blacklist)
sex_removed_tile_cov_gc_blacklist_newMask= sex_removed_tile_cov_gc %>% filter(blacklist  == 0)

#in_blacklist = (blacklist_new %>% filter(blacklist > 0))$tile

#length(intersect(in_blacklist,ext_blacklist))
