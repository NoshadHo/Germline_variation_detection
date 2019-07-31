#heatmap of all normals
mat = tile_cov_gc_normalized_227[1:5000,]
mat = tile_lr[240051:249085,]
mat = tile_lr[1:1000,]
mat2 = tile_lr[1001:2000,]
ht = Heatmap(mat,cluster_rows = F,show_row_dend = F,show_column_dend = F)
draw(ht)


#plot lr vs coordinates:
lr_plot = function(patient_num = 10,file1,start,end,chr){
  file_num = patient_num
  #file = readRDS(paste("./",files[file_num],"/",files[file_num],".rds",sep = ""))   #use read_rds from readr next time
  file = file1
  print(1)
  file$tile$n.cov.variance = variance_raw$variance
  # file$tile$n.cov.range = variance_raw$range
  # file$tile$n.IQR = coverage_raw$IQR
  # file$tile$n.var.zscore = variance_raw$var.sigma.dist
  # file$tile$n.cov.zscore = coverage_raw$cov.sigma.dist
  # file$tile$n.peak.dist = kmedian_dist$peak_dist
  file$tile$blacklist.2 = as.factor(variance_raw$blacklist)
  
  print(2)
  # data = as.data.frame(file$tile) %>% select(seqnames,start,end,t.cov,n.cov,lr,n.cov.variance,n.peak.dist,seg,blacklist.2,arm)
  data = as.data.frame(file$tile) %>% dplyr::select(seqnames,start,end,t.cov,n.cov,lr,n.cov.variance,seg,arm,blacklist.2)
  data_seg = as.data.frame(file$seg) %>% mutate(seg = row_number()) %>% dplyr::select(-strand)
  data = data %>% left_join(data_seg,by = "seg")
  data = data %>% group_by(seg) %>% mutate(seg_lr = mean(lr,na.rm = T)) %>% ungroup()
  data = data %>% group_by(seg) %>% mutate(seg_lr_med = median(lr,na.rm = T)) %>% ungroup()
  
  print(3)
  p = data %>% filter(n.cov.variance > 0.00025 & n.cov.variance < 100.1 & seqnames.x == chr) %>% ggplot()+
    facet_grid(.~seqnames.x, scales="free_x")+
    geom_point(size = 0.3,aes(x = start.x,y = lr, color = n.cov.variance > 0.020))+
    geom_segment(aes(x = start.y,xend = end.y,y = seg_lr,yend = seg_lr),colour = 'red')+
    theme_linedraw()+
    #ggtitle("sample25-chr5-, weight and No Weight, pval_diff = 5.815507e-10")+
    coord_cartesian(ylim = c(-2,2),xlim = c(start,end))+
    geom_vline(xintercept = 1)
  
  return(p)
}
lr_plot(patient_num,new_cnv,7.5e7,8e7,'chr3')

##------TEMP------------------
##to select points with norm variance of more than 0.02 from no pool, rest frompool
#then re segment,
# then plot
new_cnv = cnv_ica
new_cnv$tile$lr = if_else(cnv_org$tile$n.cov.variance < 0.02, cnv_ica$tile$lr, cnv_org$tile$lr)
new_cnv = addJointSegment(new_cnv, opts, 0)
##

segm = 2
as.data.frame(cnv$seg)[segm:(segm+1),]
as.data.frame(data) %>% group_by(seg) %>% dplyr::slice(1) %>% filter(seg %in% (segm):(segm+1)) %>% select(-seqnames.y) %>% mutate(width = width/1e4)

grid.arrange(p14_1_tn,p14_1_pool,p14_2_tn,p14_2_pool,p14_3_tn,p14_3_pool, ncol = 2)

file_num = 14
cnv_org = readRDS(paste("./",files[file_num],"/",files[file_num],".rds",sep = ""))   #use read_rds from readr next time
#14,25,101,86
#45,36,41,49

p11 = lr_plot(patient_num,file25.weighted.nokernel)
p22 = lr_plot(patient_num,file25.weighted.nokernel)
p33 = lr_plot(patient_num,file25.weighted.nokernel)
p44 = lr_plot(patient_num,file25.weighted.nokernel)
p55 = lr_plot(patient_num,file25.weighted.nokernel)
p66 = lr_plot(patient_num,file25.weighted.nokernel)

p66_w = lr_plot(patient_num,file25.weighted)
grid.arrange(p11,p22,p33,nrow = 3)
grid.arrange(p11,p22,p33,p44,p66,nrow = 5)
  
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
cnv = file1
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
cnv = addJointSegment(cnv,opts)
    
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
    time = system.time({new_var = foreach(i = 1:length(segments)) %dopar%{
    #for (i in 1:length(segments)){
      segment = segments[[i]]
      variance_kernel = list()
      if (length(segment) > 2*k + 1){
        for (j in 1:length(segment)){
          window = as.data.frame(segment[((max(1,j-k)):min(length(segment),j+k))])
          weights = (c(1:k,k+1,(k):1)/(k+1))
          if (j < k+1){
            window = window %>% mutate(weights = weights[(ceiling(length(weights)/2)-j+1):min(length(weights),length(segment))])
          }else if (j > length(segment)-k){
            window = window %>% mutate(weights = weights[1: min((length(weights)+(length(segment)-k)-j),length(segment))])
          }else{
            window = window %>% mutate(weights = weights)
          }
          variance_kernel[j] = window %>% transmute(weighted_variance = weights * n.cov.variance) %>% summarise(var = mean(weighted_variance))
          #kprint(j)        
          }
      }else{
        for (j in 1:length(segment)){
          variance_kernel[j] = segment$n.cov.variance[j]
        }
      }
      return(unlist(variance_kernel))
    }
    })
  }else if(kernel == "Gaussian"){
    
  }
  output = unlist(new_var)
}

cnv$tile$lr.weight = (1/variance_raw$variance)
file25.weighted.nokernel = addJointSegment(cnv,opts)

cnv$tile$lr.weight = (1/output)
file25.weighted.kernel = addJointSegment(cnv,opts)


#measure:
data = as.data.frame(cnv$tile) %>% select(seqnames,start,end,t.cov,n.cov,lr,n.cov.variance,n.peak.dist,seg,blacklist.2)
data_seg = as.data.frame(cnv$seg) %>% mutate(seg = row_number()) %>% select(-strand)
data = data %>% left_join(data_seg,by = "seg")
data = data %>% group_by(seg) %>% mutate(seg_lr = mean(lr,na.rm = T)) %>% ungroup()

data = data %>% filter(n.cov.variance > 0.0025) %>% rowwise() %>% 
  mutate(breakpoint.dist = min((start.x - start.y),(end.y - end.x)))

a = (data$breakpoint.dist)
a[a<0]

data %>% filter(breakpoint.dist < 0) %>% select(-lr,-seg_lr,-blacklist.2,-seqnames.y) %>% group_by(seqnames.x,seg) %>% summarise(n())



#calculate the pvalues of every 2 adjacent segments:---------------------------------------
seg_count = length(cnv$seg)
pvals.t = numeric()
pvals.wc = numeric()
for (i in 1:(seg_count-1)){
  pvals.t[i] = t.test((cnv$tile$lr[cnv$tile$seg == i]),cnv$tile$lr[cnv$tile$seg == i+1])$p.value
  pvals.wc[i] = wilcox.test((cnv$tile$lr[cnv$tile$seg == i]),cnv$tile$lr[cnv$tile$seg == i+1])$p.value
  print(i)  
}
ggplot()+geom_histogram(aes(x = pvals.t),bins = 100)

#index them
pvals.t = as.data.frame(pvals.t) %>% mutate(seg = row_number())
pvals.wc = as.data.frame(pvals.wc) %>% mutate(seg = row_number())

pvals.t %>% filter(pvals.t > 0.05) %>% arrange(desc(pvals.t))

segm = 120
as.data.frame(cnv$seg)[segm:(segm+1),]

as.data.frame(data) %>% group_by(seg) %>% dplyr::slice(1) %>% filter(seg %in% segm:(segm+1)) %>% select(-seqnames.y)

#calculate the p-value of each segment, comparing with the null model of lr:------------------------

seg_count = length(cnv$seg)
null.pvals.t = numeric()
for (i in 1:(seg_count)){
  null.pvals.t[i] = t.test((cnv$tile$lr[cnv$tile$seg == i]),lr_null_regions$lr)$p.value
  print(i)  
}

#now merge every 2 segment and measure the new p-value:
null.pvals.t.2 = numeric()
null.pvals.t.2.change = numeric()
for (i in 1:(seg_count-1)){
  new_seg = cnv$tile$lr[cnv$tile$seg == i | cnv$tile$seg == i+1]
  null.pvals.t.2[i] = t.test(new_seg,lr_null_regions$lr)$p.value
  null.pvals.t.2.change[i] = min(null.pvals.t[i+1],null.pvals.t[i]) - null.pvals.t.2[i]
  null.pvals.t.2[i] = ifelse(null.pvals.t.2[i] < null.pvals.t[i] & null.pvals.t.2[i] < null.pvals.t[i+1],1,0)
  print(i)  
}

#index them
null.pvals.t = as.data.frame(null.pvals.t) %>% mutate(seg = row_number())
null.pvals.t.2 = as.data.frame(null.pvals.t.2) %>% mutate(null.pvals.t.2.change) %>% mutate(seg = row_number())
null.pvals.t.2 %>% arrange(desc(null.pvals.t.2.change))

#compare with previous method, see how many of the selected segments are mutual: -----
seg1 = pvals.t %>% filter(pvals.t > 0.05)
seg1 = seg1$seg
seg2 = null.pvals.t.2 %>% filter(null.pvals.t.2 == 1)
seg2 = seg2$seg
intersect(seg1,seg2)
setdiff(seg2,seg1)

##--------------------------Measuring variance around a couple of selected segments for pca ica and matched tumor normal--------------
as.data.frame(cnv_pca$seg) %>% mutate(seg = row_number()) %>% filter(seqnames == "chr12") #selecting using this code
#selecting the seg:
cnv_org$tile$n.cov.variance = variance_raw$variance
seg1_tn = cnv_org$tile$lr[cnv_org$tile$seg == 260 & cnv_org$tile$n.cov.variance < 0.02] #seg15 on original cnv
seg1_ica = cnv_ica$tile$lr[cnv_ica$tile$seg == 428 & cnv_ica$tile$n.cov.variance < 0.02] #seg25 on ica cnv
seg1_pca = cnv_pca$tile$lr[cnv_pca$tile$seg == 418 & cnv_pca$tile$n.cov.variance < 0.02] #seg25 on pca cnv
sd(seg1_tn,na.rm = TRUE)
sd(seg1_ica,na.rm = TRUE)
sd(seg1_pca,na.rm = TRUE)

#second seg
seg2_tn = cnv_org$tile$lr[cnv_org$tile$seg == 58] #seg15 on original cnv
seg2_ica = cnv_ica$tile$lr[cnv_ica$tile$seg == 79] #seg25 on ica cnv
seg2_pca = cnv_pca$tile$lr[cnv_pca$tile$seg == 198] #seg25 on pca cnv
sd(seg2_tn,na.rm = TRUE)
sd(seg2_ica,na.rm = TRUE)
sd(seg2_pca,na.rm = TRUE)
