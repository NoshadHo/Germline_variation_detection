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


#unmasked plot of data:
tile_umasked = tile_unmasked %>% select(!!1)
tile_umasked = tile_umasked %>% mutate(blacklist_new$blacklist)
var = variance_sex(as.data.frame(t(tile_cov_gc_normalized_227)))
tile_umasked = tile_umasked %>% mutate(var = var[,1])
colnames(tile_umasked) = c('unmasked','blacklist','var')

var_unmasked_plot = function(tile_umasked){
scatterPlot = tile_umasked %>% ggplot()+geom_point(aes(x = unmasked,y = var, color = blacklist))+ylim(0,7)
xdensity = tile_umasked %>% ggplot()+geom_density(aes(x = unmasked, fill = blacklist))
ydensity = tile_umasked %>% ggplot()+geom_density(aes(x = var, fill = blacklist))+xlim(1,7.5)+coord_flip()
blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )

grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
}

#remove sex chromosomes
tile_umasked = tile_umasked %>% slice(1:287509)
tile_umasked = tile_umasked %>% filter(blacklist == "normal")
var_unmasked_plot(tile_umasked)

#remove the first pc:
svd_data = tile_cov_gc_normalized_227 %>% slice(1:287509) %>% mutate(blacklist = blacklist_new[1:287509,2]) %>% 
  filter(blacklist == "normal") %>% select(-blacklist)
svd = svd(svd_data)
sc_num = 1

svd$d[1:sc_num] = 0
svd$d = diag(svd$d)
purified_tile_cov_gc_normalized = svd$u %*% tcrossprod(svd$d,svd$v)
tile_umasked = tile_unmasked %>% select(!!1)
tile_umasked = tile_umasked %>% mutate(blacklist = blacklist_new$blacklist)
var_2 = variance_sex(as.data.frame(t(purified_tile_cov_gc_normalized)))
tile_umasked = tile_umasked %>% slice(1:287509) %>% filter(blacklist == "normal") %>% mutate(var = var_2[,1])
colnames(tile_umasked) = c('unmasked','blacklist','var')
var_unmasked_plot(tile_umasked)

#lr-tile plot with blacklisted tiles
patient_lr = tile_lr %>% select(!!45)
patient_lr = patient_lr %>% mutate(blacklist = blacklist_new$blacklist)
colnames(patient_lr) = c('lr','blacklist')
patient_lr %>% slice(153652:167492) %>% ggplot(aes(x = 1:(167492-153652+1),y = lr,color = blacklist))+geom_point()
167492153652


#lr with variance coloring
  variance_df = variance_sex(as.data.frame(t()),"m")
  variance_df = variance_df2
  variance_df = variance_df %>% mutate(lr = tile_lr[,106])

  #plot for chr1
  variance_df = variance_df %>% filter(V1 < 1)
  chr_start = 106122
  chr_end = 123203
  ggplot()+geom_point(aes(x = chr_start:chr_end, y = variance_df$lr[chr_start:chr_end], color = variance_df$V1[chr_start:chr_end]))
#range analysis
  variance_df = variance_df %>% mutate(range = cov_ranges$range)
  variance_df = variance_df %>% mutate(unmasked = tile_unmasked[,1])  
  variance_df %>% ggplot(aes(x = range,y = V1,color = unmasked))+geom_point()+theme_linedraw()
  variance_df %>% ggplot(aes(x = unmasked))+geom_histogram(binwidth = 0.01)+theme_linedraw()


  
#coverage + variance plot for blacklists:-------------------------------------------
  coverage_raw = tile_cov_gc_normalized_227 %>% select(!!101)
  variance_raw = variance_sex(as.data.frame(t(tile_cov_gc_normalized_227)))
  
  start = 45001
  end = 50000
  49116
  coverage = coverage_raw %>% slice(start:end) %>% mutate(blacklist = 0)
  variance = variance_raw %>% slice(start:end) %>% mutate(blacklist = 0)
  colnames(coverage) = c("coverage", "blacklist")
  colnames(variance) = c("variance", "blacklist")
  
  

#using 10 tile variance aggregated
  start_10 = floor(start/10)+1
  end_10 = floor(end/10)
  variance.temp = variance %>% mutate(group = floor(seq(from = 1,length.out = (end-start+1),by = 0.1)))
  variance.temp = variance.temp %>% group_by(group) %>% summarise(mean(variance))
  variance.temp = variance.temp %>% select(-group)%>% mutate(blacklist = 0)
  #variance.temp = variance.temp %>% slice(start_10:end_10) %>% mutate(blacklist = 0)
  colnames(variance.temp) = c("variance", "blacklist")
  
  p1 = coverage %>% ggplot(aes(x = start:end, y = coverage, color = as.factor(blacklist))) + geom_point(size = 0.5)
  p2 = variance %>% ggplot(aes(x = start:end, y = variance, color = as.factor(blacklist))) + geom_point(size = 0.5)
  
  
  p3 = variance.temp %>% ggplot(aes(x = start_10:end_10, y = variance, color = as.factor(blacklist))) + geom_point(size = 0.5)+
    coord_cartesian(xlim = c(start_10,end_10))+geom_hline(yintercept = 0.06)+coord_cartesian(ylim = c(0,1))
  grid.arrange(p1,p2,p3)
  
  variance.temp[408:413,2] = 1
  coverage[4080:4130,2] = 1
  variance[4080:4130,2] = 1
  
  