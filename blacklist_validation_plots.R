chromosome_standard_devation_plot = function(data, chr = 0, range = NULL, write_address = '~/Codes/Germline_variation_detection/chr_std_plot'){
  # data format should look like: [1] 308837    227
  #variance_sex function is in Variance_analysis
  variance_df = variance_sex(as.data.frame(t(data)))
  variance_df = variance_df %>% dplyr::mutate(blacklist = blacklist_new$blacklist)
  colnames(variance_df) = c("variance", "blacklist")
  
  if (chr == 0){
    
    for (chr in 1:22){
      start = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                           ungroup() %>% slice(chr) %>% select(tile))
      end = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                         ungroup() %>% slice(chr+1) %>% select(tile))
    
      pdf(paste(write_address,"_",chr,".pdf",sep = ""))
      ggplot()+geom_point(aes(x = start:end,y = variance_df$variance[start:end], color = variance_df$blacklist[start:end]),size = 0.4)+theme_linedraw()+
        labs(title = paste("Chr",chr),x = "Genome tiles", y = "Variance")+ylim(0,4)
      dev.off()
    }
  }else{
    start = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
      ungroup() %>% slice(chr) %>% select(tile))
    end = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                       ungroup() %>% slice(chr+1) %>% select(tile))
    ggplot()+geom_point(aes(x = start:end,y = variance_df$variance[start:end], color = variance_df$blacklist[start:end]),size = 0.4)+theme_linedraw()+
      labs(title = paste("PC removed",sc_num),x = "Genome tiles", y = "Variance")
  }

}

data = tile_cov_gc_normalized_227

#for ploting 6 in one plot----------------------------
chr = 1
start = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                     ungroup() %>% slice(chr) %>% select(tile))
end = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                   ungroup() %>% slice(chr+1) %>% select(tile))
p1 = ggplot()+geom_point(aes(x = start:end,y = variance_df$variance[start:end], color = variance_df$blacklist[start:end]),size = 0.4)+theme_linedraw()+
  labs(title = paste("Chr",chr),x = "Genome tiles", y = "Variance")+theme(legend.title = element_blank())

start_2 = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                     ungroup() %>% slice(chr+1) %>% select(tile))
end_2 = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                   ungroup() %>% slice(chr+2) %>% select(tile))
p2 = ggplot()+geom_point(aes(x = start_2:end_2,y = variance_df$variance[start_2:end_2], color = variance_df$blacklist[start_2:end_2]),size = 0.4)+theme_linedraw()+
  labs(title = paste("Chr",chr+1),x = "Genome tiles", y = "Variance")+theme(legend.title = element_blank())

start_3 = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                       ungroup() %>% slice(chr+2) %>% select(tile))
end_3 = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                     ungroup() %>% slice(chr+3) %>% select(tile))
p3 = ggplot()+geom_point(aes(x = start_3:end_3,y = variance_df$variance[start_3:end_3], color = variance_df$blacklist[start_3:end_3]),size = 0.4)+theme_linedraw()+
  labs(title = paste("Chr",chr+2),x = "Genome tiles", y = "Variance")+theme(legend.title = element_blank())

start_4 = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                       ungroup() %>% slice(chr+3) %>% select(tile))
end_4 = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                     ungroup() %>% slice(chr+4) %>% select(tile))
p4 = ggplot()+geom_point(aes(x = start_4:end_4,y = variance_df$variance[start_4:end_4], color = variance_df$blacklist[start_4:end_4]),size = 0.4)+theme_linedraw()+
  labs(title = paste("Chr",chr+3),x = "Genome tiles", y = "Variance")+theme(legend.title = element_blank())

start_5 = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                       ungroup() %>% slice(chr+4) %>% select(tile))
end_5 = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                     ungroup() %>% slice(chr+5) %>% select(tile))
p5 = ggplot()+geom_point(aes(x = start_5:end_5,y = variance_df$variance[start_5:end_5], color = variance_df$blacklist[start_5:end_5]),size = 0.4)+theme_linedraw()+
  labs(title = paste("Chr",chr+4),x = "Genome tiles", y = "Variance")+theme(legend.title = element_blank())

start_6 = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                       ungroup() %>% slice(chr+5) %>% select(tile))
end_6 = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                     ungroup() %>% slice(chr+6) %>% select(tile))
p6 = ggplot()+geom_point(aes(x = start_6:end_6,y = variance_df$variance[start_6:end_6], color = variance_df$blacklist[start_6:end_6]),size = 0.4)+theme_linedraw()+
  labs(title = paste("Chr",chr+5),x = "Genome tiles", y = "Variance")+theme(legend.title = element_blank())

grid.arrange(p1,p2,p3,p4,p5,p6)
#--------------------------------------------------

#box plot for blacklisted regions vs. none blacklisted regions
ggplot()+geom_boxplot(aes(x = variance_df$blacklist,y = variance_df$variance),size = 0.4)+theme_linedraw()+
  labs(title = paste("PC removed",sc_num),x = "Genome tiles", y = "Variance")+ylim(0,0.1)

#frequency of different sizes of blacklist (2 tile consequnsive, 3tile ...)
#read the bed file made earlier of blacklist:
blacklist_bed = read_tsv("~/Codes/Germline_variation_detection/blacklist3_withXY.bed",col_names = FALSE)
colnames(blacklist_bed) = c('seqnames','start','end')
blacklist_bed = GenomicRanges::makeGRangesFromDataFrame(blacklist_bed)
blacklist_bed = as.data.frame(GenomicRanges::reduce(blacklist_bed))
blacklist_bed_tile_count = blacklist_bed %>% select(width) %>% mutate(tile_count = floor(width/10000)+1)
blacklist_bed_tile_count = blacklist_bed_tile_count %>% group_by(tile_count) %>% summarise(freq = n())
blacklist_bed_tile_count %>% ggplot(aes(x = tile_count,y = freq))+geom_point()+geom_smooth()


