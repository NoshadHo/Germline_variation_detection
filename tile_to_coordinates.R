tile_to_coordinates = function(tile,file){ #RDSfile is a one of the .RDS files of patients, don't matter which one
  range = as.data.frame(file$tile[tile]) %>% select(seqnames,start,end,blacklist)
  return(range)
}


tile_to_bed(tiles,output_Address){ #example for tiles: tiles = blacklist_new_2 %>% filter(blacklist > 0)
  output = as.data.frame(matrix(nrow = dim(tiles)[1],ncol = 3))
  for (i in 1:dim(tiles)[1]){
    output[i,] = as.data.frame(tile_to_coordinates(tiles[i,1],file)[,1:3])
    print(i)
  }
  output = output %>% mutate(V1 = paste("chr",as.character(V1),sep = ""))
  output = output %>% mutate(V1 = if_else(V1=='chr23',"chrX",V1))
  output = output %>% mutate(V1 = if_else(V1=='chr24',"chrY",V1))
  colnames(output) = c("seqnames", "start", "end")
  GRoutput = makeGRangesFromDataFrame(output)
  GRoutput = GenomicRanges::reduce(GRoutput)
  output = as.data.frame(GRoutput) %>% select(seqnames,start,end)
  #write file
  write_tsv(output,output_Address,col_names = FALSE)
}
as.data.frame(do.call(rbind,output))

coordinates_to_tile = function(chr,start_coord,end_coord,file){ #file_tile is the list of tiles fro file = file$tile
  range = as.data.frame(file$tile)
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
      tiles[i] = as.data.frame((coordinates_to_tile(bed_file$chr[i], bed_file$start[i],bed_file$end[i],RDSfile)))
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

# new_blacklist_tiles = c(1:90,260:280,12000:12341,1278:1340,1650:1697,24841:24871,
#                         22861:22865,12341:14341,16144:16148,14341:15000,
#                         257422:257422,33628:33640,25072,28183,28188,
#                         34287:34367,34477:34505,34600:34660,35890:35965,36247:36262,36032:36052,35816,37887:37987,36032:36052,
#                         56667:56695,58167:58207,
#                         68657:68727,68807:68947,
#                         73850:73947,73947:74147,
#                         87910:87969,
#                         68947:68957,
#                         69855:69905,87969:88090,
#                         92580:92849,
#                         89710:89730,
#                         90099:90119,92849:93020,
#                         94820:95150,
#                         105549:106123,106123:106173,112103:112303,
#                         111953:112103,123103:123204,
#                         109343:109423,122165:122200,123204:123254,129214:129564,
#                         128830:129214,129710:129814,
#                         130464:130564,139139:139250,143659:143800,153653:153703,157953:160500,171293:171473,
#                         139800:140000,153553:153653,157503:157953,167453:167493,167493:167550,
#                         140339:140439,147690:147730,
#                         130654:130754,
#                         133434:133504,
#                         137550:137700,171473:171753,180873:180920,194300:194382,194382:194420,197932:198200,218875:218895,220867:221150,231752:232252,241510:241590,249040:249086,249086:249186,251596:251800,257412:257442,265300:265450,267870:268070,268070:268230,273850:274122,274122:274460,277757:278957,278957:279157,282428:283928,283928:284600,
#                         180673:180873,189820:189900,207650:207710,209480:209580,229650:229852,239952:240052,241650:241730,251200:251596,252700:252750,258900:259262,265450:265520,265450:265520,277700:277757,282000:282250,287480:287510,
#                         172050:172370,238220:238332,241852:241932,250900:250960,252850:252940,269075:269081,281150:281250,
#                         232652:232752,242172:242332,253630:253760,
#                         232820:232952,242902:242952,253405:253420,
#                         233000:233152,243202:243732)


new_blacklist_tiles = c(1:180,24952, 56650:56700, 69800:69920, 89710:89750, 111960:112120, 128840:129600, 140330:140420, 157500:158230, 171320:171760, 185950:186460, 197850:198170, 209300:209600, 220750:221430, 231500:232350, 241510:241750, 249050:249150, 258910:259520, 266320:266380, 273880:274470, 278210:278870, 283470:283750,
                        2679:2681,25187,58170:58500, 73850:74140, 92920:93000,120396,133410:133490, 143070:143110, 159680:160520, 172060:172320, 180292:180293, 200820:200825,216264,	229700:229860, 232670:233320, 241980:241990, 250730:251810, 262080:262130, 267880:267960, 279040:279170, 283910:285010,
                        10350:10390,28188,61190:61195,72551,92560:92710,120447,137550:137650, 143490:143750,161737, 174807,189040:189045, 204425:204427, 240010:240080, 242190:242340, 252680:252750, 265410:265510, 268160:268230,
                        12680:12690,44825,65210:65215,73510,94920:95120, 138570:138573, 147700:147740, 161673:161675,175750,188880:188885, 242880:243090, 252850:252920, 269070:269080,
                        14310:14980, 33580:34670, 65390:65410, 78415:78420,96148,139830:139990, 150047:150050, 180830:180910, 193975:193980, 243240:243910, 253380:253430, 269170:269180,
                        24820:24900, 49080:49130, 65810:65815, 87850:88070,96353,244650:244720, 253630:253790, 270450:270480,
                        22850:22870, 68660:68730, 106070:106220, 257390:257460,
                        24025:24030, 68920:68970,11960:12530)
# blacklist_new = blacklist %>% mutate(blacklist = case_when(blacklist > 0 ~ 'blacklist',
#                                            tile %in% new_blacklist_tiles ~ 'blacklist',
#                                            TRUE ~ 'normal'))
blacklist_new = blacklist %>% select(tile) %>% mutate(blacklist = if_else(tile %in% new_blacklist_tiles,'blacklist','normal'))

#make it sex exclusive
blacklist_new = blacklist_new[1:287509,]
sex_removed_tile_cov_gc = sex_removed_tile_cov_gc %>% mutate(blacklist = blacklist_new$blacklist)
blacklist_removed_tile_list_newMask = sex_removed_tile_cov_gc %>% filter(blacklist == 0) %>% select(tile)
blacklist_removed_tile_list_newMask = blacklist_removed_tile_list_newMask %>% mutate(new_row = row_number())
sex_removed_tile_cov_gc_blacklist_newMask= sex_removed_tile_cov_gc %>% filter(blacklist  == 0)
sex_removed_tile_cov_gc_blacklist_newMask = sex_removed_tile_cov_gc_blacklist_newMask %>% select(-tile,-blacklist)
#in_blacklist = (blacklist_new %>% filter(blacklist > 0))$tile

#length(intersect(in_blacklist,ext_blacklist))

