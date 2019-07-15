setwd("~/Data/cnvex_normals/cnvex/")
#get a list of folders
files = list.files()

output = as.data.frame(matrix(nrow = 0,ncol = 8))
for (i in 171:227){
  print(i)
  #finding small segments (less than 5 tiles)
  file_num = i
  file = read_rds(paste("./",files[file_num],"/",files[file_num],".rds",sep = ""))   #use read_rds from readr next time  
  file$tile$n.cov.variance = variance_raw$variance
  file$tile$n.cov.range = variance_raw$range
  file$tile$n.IQR = coverage_raw$IQR
  file$tile$n.var.zscore = variance_raw$var.sigma.dist
  file$tile$n.cov.zscore = coverage_raw$cov.sigmaa.dist
  file$tile$n.peak.dist = kmedian_dist$peak_dist
  file$tile$blacklist.2 = as.factor(variance_raw$blacklist)
  
  data = as.data.frame(file$tile) %>% 
    select(seqnames,start,end,t.cov,n.cov,lr,n.cov.variance,n.cov.range,n.IQR,n.peak.dist,n.var.zscore,n.cov.zscore,seg,blacklist.2)
  data_seg = as.data.frame(file$seg) %>% mutate(seg = row_number()) %>% select(-strand)
  data = data %>% left_join(data_seg,by = "seg")
  data = data %>% group_by(seg) %>% mutate(seg_lr = mean(lr,na.rm = T)) %>% ungroup()
  
  small_segments = (data_seg %>% filter(width < 50001))$seg
  data_small_segments = data %>% filter(seg %in% small_segments)
  #define it's features
  x = (data_small_segments %>% select(n.cov.variance,n.cov.range,n.IQR,n.peak.dist,n.var.zscore,n.cov.zscore,blacklist.2,seg_lr))
  #y = as.matrix(data_small_segments %>% select(seg_lr))
  
  boruta = Boruta(seg_lr ~ .,data = x, doTrace = 0,maxRuns = 200)
  roughFixMod <- TentativeRoughFix(boruta)
  one_sample_results = attStats(roughFixMod)
  one_sample_results = one_sample_results %>% mutate(sample = file_num) %>% mutate(features = rownames(one_sample_results))
  output = rbind(output,one_sample_results)
}
save.image("~/Codes/Germline_variation_detection/data_july11.RData")

#analysis of random forest results:

#mean importance of each feature:
output %>% group_by(features) %>% summarise(imp = mean(meanImp)) %>% arrange(desc(imp))

#add confirmed as an factor for decision:
output = output %>% mutate(confirmed = if_else(decision == "Confirmed",1,0))
output %>% group_by(features) %>% mutate(wconf = confirmed*meanImp) %>% summarize(res = mean(wconf)) %>% arrange(desc(res))

small_segments = (data_seg %>% filter(width < 50001))$seg
data_small_segments = data %>% filter(seg %in% small_segments)


#find the variability of each tile (its segment value) and then correlate it with other stuff
tile_seg_lr = list()
for (file_num in 1:(length(files))){ #can be potentially multithreat
  file = read_rds(paste("./",files[file_num],"/",files[file_num],".rds",sep = ""))   #use read_rds from readr next time
  tiles = as.data.frame(file$tile)
  
  data = tiles %>% 
    select(seqnames,start,end,t.cov,n.cov,lr,seg)
  data_seg = as.data.frame(file$seg) %>% mutate(seg = row_number()) %>% select(-strand)
  data = data %>% left_join(data_seg,by = "seg")
  data = data %>% group_by(seg) %>% mutate(seg_lr = mean(lr,na.rm = T)) %>% ungroup()
  
  tile_seg_lr[file_num] = data %>% select(seg_lr)
  colnames_list[file_num] = files[file_num]
  #sex[file_num,1] = file_num
  #sex[file_num,2] = detect.sex(file$var,file$tile)
  print(paste("File processed: ",file_num,"/",length(files),sep = ""))
}
tile_seg_lr = as.data.frame(matrix(unlist(tile_seg_lr),ncol = length(tile_seg_lr), byrow = FALSE))
colnames(tile_seg_lr) = colnames_list

variance_lr = variance_sex(as.data.frame(t(tile_seg_lr)))
rfInput = variance_lr %>% mutate(n.cov.variance = variance_raw$variance)%>% mutate(n.cov.range = variance_raw$range) %>% 
  mutate(n.IQR = coverage_raw$IQR) %>% mutate(n.var.zscore = variance_raw$var.sigma.dist) %>% mutate(n.cov.zscore = coverage_raw$cov.sigmaa.dist) %>% 
  mutate(n.peak.dist = kmedian_dist$peak_dist) %>% mutate(blacklist.2 = variance_raw$blacklist)
colnames(rfInput)[1] = "lr.seg.var"

save.image("~/Codes/Germline_variation_detection/data_july11.RData")
boruta = Boruta(lr.seg.var ~ .,data = rfInput, doTrace = 1, maxRuns = 40)
save.image("~/Codes/Germline_variation_detection/data_july11.RData")
names(boruta)
getSelectedAttributes(boruta, withTentative = TRUE)
roughFixMod <- TentativeRoughFix(boruta)
imps <- attStats(roughFixMod)
imps2 = imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
plot(boruta, cex.axis=.7, las=2, xlab="", main="Variable Importance")  


a = as.data.frame(unique(variance_lr$V1))

colnames(a) = 'V1'
#####
a %>% ggplot() +geom_density(aes(x = V1),binwidth = 0.0001)+ggtitle("Seg_lr_variance")
#####

# Fit random forrest
boruta = Boruta(seg_lr ~ .,data = x, doTrace = 1)
names(boruta)
getSelectedAttributes(boruta, withTentative = TRUE)
roughFixMod <- TentativeRoughFix(boruta)
imps <- attStats(roughFixMod)
imps2 = imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
plot(boruta, cex.axis=.7, las=2, xlab="", main="Variable Importance")  
