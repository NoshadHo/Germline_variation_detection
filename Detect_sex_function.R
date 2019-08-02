detect_sex_new = function(var, tile) {
  tile$n.cov = tile$t.cov/(2^tile$lr)
  n.cov.a = (as.data.frame(tile) %>% filter(seqnames != "chrX" & (seqnames != "chrY")))$n.cov
  n.cov.x = (as.data.frame(tile) %>% filter(seqnames == "chrX"))$n.cov
  n.cov.y = (as.data.frame(tile) %>% filter(seqnames == "chrY"))$n.cov
  THRESHOLD = 0.05 
  if (findInterval(x=(median(n.cov.x,na.rm = T)/median(n.cov.a, na.rm = T)), c(1-THRESHOLD, 1+THRESHOLD)) == 1L) {
    return("female")
  }else{
    if (findInterval(x=(median(n.cov.x,na.rm = T)/median(n.cov.a, na.rm = T)), c(0.5-THRESHOLD, 0.5+THRESHOLD)) == 1L){
      return("male")
    }
    if (findInterval(x=abs(median(n.cov.y,na.rm = T)/median(n.cov.a, na.rm = T)), c(0, THRESHOLD)) == 1L){
      return("female")
    }
  }
}


sex_results = foreach(i = 1:227) %dopar% {
  cnv = readRDS(paste("./",files[i],"/",files[i],".rds",sep = ""))   #use read_rds from readr next time
  cat(as.character(i),file="progress.txt",sep="\n",append=TRUE)
  return(list(bafBased = detect.sex(var = cnv$var,tile = cnv$tile), tlrBased = detect_sex_new(var = cnv$var,tile = cnv$tile)))
}
baf = numeric()
tlr = numeric()
for (i in 1:227) {
  if(!is.null(sex_results[[i]])){
    baf[i] = sex_results[[i]][[1]]  
    tlr[i] = sex_results[[i]][[2]]
  }
}
sex_results_df = as.data.frame(cbind(baf, tlr))

sex_stat = foreach(i = 1:227) %dopar% {
  cnv = readRDS(paste("./",files[i],"/",files[i],".rds",sep = ""))   #use read_rds from readr next time
  cnv$tile$n.cov = cnv$tile$t.cov/(2^cnv$tile$lr)
  n.cov.auto = (as.data.frame(cnv$tile) %>% filter(seqnames != "chrX") %>% filter(seqnames != "chrY"))$n.cov
  n.cov.x = (as.data.frame(cnv$tile) %>% filter(seqnames == "chrX"))$n.cov
  n.cov.y = (as.data.frame(cnv$tile) %>% filter(seqnames == "chrY"))$n.cov
  sex = detect.sex(var = cnv$var,tile = cnv$tile)
  cat(as.character(i),file="progress.txt",sep="\n",append=TRUE)
  return(list(sex = sex,
              x.median = median(n.cov.x, na.rm = T),
              y.median = median(n.cov.y, na.rm = T),
              auto.median = median(n.cov.auto, na.rm = T)))
}

sex = character()
x.median = numeric()
y.median = numeric()
auto.median = numeric()
for (i in 1:227) {
  if(!is.null(sex_stat[[i]])){
    sex[i] = sex_stat[[i]][[1]]  
    x.median[i] = sex_stat[[i]][[2]]
    y.median[i] = sex_stat[[i]][[3]]
    auto.median[i] = sex_stat[[i]][[4]]
  }
}
sex_stat = as.data.frame(cbind(sex, x.median,y.median,auto.median))

max((sex_stat %>% mutate(xTOa = (as.numeric(as.character(y.median))/as.numeric(as.character(x.median)))) %>% filter(sex == "female"))$xTOa)


#plot lr vs coordinates:
lr_plot_2 = function(patient_num = 10,file1,start,end,chr){
  file = file1
  file$tile$n.cov.variance = variance_raw$variance
  file$tile$blacklist.2 = as.factor(variance_raw$blacklist)
  
  data = as.data.frame(file$tile) %>% dplyr::select(seqnames,start,end,t.cov,n.cov,lr,n.cov.variance,seg,arm,blacklist.2,baf)
  data_seg = as.data.frame(file$seg) %>% mutate(seg = row_number()) %>% dplyr::select(-strand)
  data = data %>% left_join(data_seg,by = "seg")
  data = data %>% group_by(seg) %>% mutate(seg_lr = mean(lr,na.rm = T)) %>% ungroup()
  data = data %>% group_by(seg) %>% mutate(seg_lr_med = median(lr,na.rm = T)) %>% ungroup()
  #data$n.cov = n.cov = 1/((2^data$lr)/data$t.cov)
  p = data %>% dplyr::filter(n.cov.variance > 0.00025 & n.cov.variance < 100.1 & seqnames.x == chr) %>% ggplot()+
    facet_grid(.~seqnames.x, scales="free_x")+
    geom_point(size = 0.3,aes(x = start.x,y = n.cov, color = n.cov.variance > 0.020))+
    geom_hline(yintercept = mean(data$n.cov, na.rm = T))+
    coord_cartesian(ylim = c(-1,1))+
    geom_vline(xintercept = 1)
  
  return(p)
}
lr_plot_2(patient_num,cnv_92,1,8e7,c('chr22','chrX',"chrY"))
