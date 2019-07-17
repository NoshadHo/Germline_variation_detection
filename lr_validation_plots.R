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
  p = data %>% filter(n.cov.variance > 0.025 & n.cov.variance < 10.1 & seqnames.x == 'chr10') %>% ggplot()+
    facet_grid(.~seqnames.x, scales="free_x")+
    geom_point(size = 0.3,aes(x = start.x,y = lr))+
    geom_segment(aes(x = start.y,xend = end.y,y = seg_lr,yend = seg_lr),colour = 'red')+
    theme_linedraw()#+ggtitle("sample36-chr19-coord:whole Chr, No Weight")+coord_cartesian(ylim = c(-1,1),xlim = c(3e7,4e7))
  
  return(p)
}
file_num = 102
file102 = readRDS(paste("./",files[file_num],"/",files[file_num],".rds",sep = ""))   #use read_rds from readr next time
#14,25,101,86
#45,36,41,49

p11 = lr_plot(patient_num,file67)
p22 = lr_plot(patient_num,file102)
p33 = lr_plot(patient_num,file36)
p44 = lr_plot(patient_num,file45)
p55 = lr_plot(patient_num,file49)
p66 = lr_plot(patient_num,file25)
p66_w = lr_plot(patient_num,file25.weighted)
grid.arrange(p66,p66_w)
grid.arrange(p11,p22,p33,p44,p55,p66,nrow = 6)
  
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
cnv = file36
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
file25.weighted = addJointSegment(cnv,opts)
    
#calculate the distance of each tile to the nearest breakpoint

#defining new weights using kernel.
#kernel1: linear
variance.withkernel = function(cnv = file25,kernel = "Gaussian",k = 5, intial_seg = TRUE){
  if (initial_seg == TRUE){
    cnv = addJointSegment(cnv,opts)
  }
  output = numeric()
  segments = split(cnv$tile,cnv$tile$seg)
  if (kernel == "Linear"){
    time = system.time({foreach(i = 1:length(segments)) %dopar%{
      segment = segments[[i]]
      variance_kernel = list()
      if (length(segment) > k){
        for (j in 1:length(segment)){
          window = as.data.frame(segment[((max(1,j-k)):min(length(segment),j+k))])
          weights = (c(1:k,k+1,(k):1)/(k+1))
          if (j < k+1){
            window = window %>% mutate(weights = weights[(ceiling(length(weights)/2)-j+1):length(weights)])
          }else if (j > length(segment)-k){
            window = window %>% mutate(weights = weights[1: min((length(weights)+(length(segment)-k)-j),length(segment))])
          }else{
            window = window %>% mutate(weights = weights)
          }
          variance_kernel[j] = window %>% transmute(weighted_variance = weights * n.cov.variance) %>% summarise(var = mean(weighted_variance))
        }
      }else{
        for (j in 1:length(segment)){
          variance_kernel[j] = segment$n.cov.variance
        }
      }
      return(unlist(variance_kernel))
    }
    })
  }else if(kernel == "Gaussian"){
    
  }
}
