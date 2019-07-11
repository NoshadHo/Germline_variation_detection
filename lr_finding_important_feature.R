setwd("~/Data/cnvex_normals/cnvex/")
#get a list of folders
files = list.files()

output = as.data.frame(matrix(nrow = 0,ncol = 8))
for (i in 1:227){
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

# Fit the LASSO model (Lasso: Alpha = 1)
cv.lasso <- cv.glmnet(x, y, alpha=1, parallel=TRUE, standardize=TRUE, type.measure='auc')
coef(cv.lasso)
# Fit random forrest
boruta = Boruta(seg_lr ~ .,data = x, doTrace = 1)
names(boruta)
getSelectedAttributes(boruta, withTentative = TRUE)
roughFixMod <- TentativeRoughFix(boruta)
imps <- attStats(roughFixMod)
imps2 = imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
plot(boruta, cex.axis=.7, las=2, xlab="", main="Variable Importance")  
