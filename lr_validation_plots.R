#heatmap of all normals
mat = tile_cov_gc_normalized_227[1:5000,]
mat = tile_lr[240051:249085,]
mat = tile_lr[1:1000,]
mat2 = tile_lr[1001:2000,]
ht = Heatmap(mat,cluster_rows = F,show_row_dend = F,show_column_dend = F)
draw(ht)


#plot lr vs coordinates:
lr_plot = function(patient_num = 10,file1){
  file_num = patient_num
  #file = readRDS(paste("./",files[file_num],"/",files[file_num],".rds",sep = ""))   #use read_rds from readr next time
  file = file1
  print(1)
  file$tile$n.cov.variance = variance_raw$variance
  file$tile$n.cov.range = variance_raw$range
  file$tile$n.IQR = coverage_raw$IQR
  file$tile$n.var.zscore = variance_raw$var.sigma.dist
  file$tile$n.cov.zscore = coverage_raw$cov.sigma.dist
  file$tile$n.peak.dist = kmedian_dist$peak_dist
  file$tile$blacklist.2 = as.factor(variance_raw$blacklist)
  
  print(2)
  data = as.data.frame(file$tile) %>% select(seqnames,start,end,t.cov,n.cov,lr,n.cov.variance,n.peak.dist,seg,blacklist.2)
  data_seg = as.data.frame(file$seg) %>% mutate(seg = row_number()) %>% select(-strand)
  data = data %>% left_join(data_seg,by = "seg")
  data = data %>% group_by(seg) %>% mutate(seg_lr = mean(lr,na.rm = T)) %>% ungroup()
  
  print(3)
  p = data %>% filter(n.cov.variance > 0.025 & n.cov.variance < 10.1 & seqnames.x == 'chr19') %>% ggplot()+
    facet_grid(.~seqnames.x, scales="free_x")+
    geom_point(size = 0.3,aes(x = start.x,y = lr,color = as.factor(blacklist.2)))+
    geom_segment(aes(x = start.y,xend = end.y,y = seg_lr,yend = seg_lr),colour = 'red')+
    theme_linedraw()+ggtitle("sample36-chr19-coord:whole Chr, No Weight")#+coord_cartesian(ylim = c(-1,1),xlim = c(1.03e8,1.035e8))+ggtitle("sample41-chr6-coord:1.03e8,1.035e8, No Weight")
  
  return(p)
}
file_num = 102
file102 = readRDS(paste("./",files[file_num],"/",files[file_num],".rds",sep = ""))   #use read_rds from readr next time
#14,25,101,86
#45,36,41,49
patient_num = 25
p1 = lr_plot(patient_num,file14)
p2 = lr_plot(patient_num,file14.weighted)

p3 = lr_plot(patient_num,file25)
p4 = lr_plot(patient_num,file25.weighted)

p5 = lr_plot(patient_num,file25)
p6 = lr_plot(patient_num,file25.weighted)

p7 = lr_plot(patient_num,file36)
p8 = lr_plot(patient_num,file36.weighted)

p9 = lr_plot(patient_num,file41)
p10 = lr_plot(patient_num,file41.weighted)

p11 = lr_plot(patient_num,file36)
p12 = lr_plot(patient_num,file36.weighted) 

grid.arrange(p1,p2,p3,p4,p5,p6,nrow = 6) 
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol = 2) 


#change the segmentation:
cnv = file41
cnv$tile$n.cov.variance = variance_raw$variance
cnv$tile$n.cov.range = variance_raw$range
cnv$tile$n.IQR = coverage_raw$IQR
cnv$tile$n.var.zscore = variance_raw$var.sigma.dist
cnv$tile$n.cov.zscore = coverage_raw$cov.sigma.dist
cnv$tile$n.peak.dist = kmedian_dist$peak_dist
cnv$tile$blacklist.2 = as.factor(variance_raw$blacklist)
cnv$tile$lr.weight = (1/variance_raw$variance)
cnv$tile$seg = NULL
cnv$seg = NULL
file41.weighted = addJointSegment(cnv,opts)
    
#calculate the distance of each tile to the nearest breakpoint