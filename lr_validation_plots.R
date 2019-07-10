#heatmap of all normals
mat = tile_cov_gc_normalized_227[1:5000,]
mat = tile_lr[240051:249085,]

ht = Heatmap(mat,cluster_rows = F,show_row_dend = F,show_column_dend = F)
draw(ht)


#plot lr vs coordinates:
lr_plot = function(patient_num = 10){
  file_num = patient_num
  file = readRDS(paste("./",files[file_num],"/",files[file_num],".rds",sep = ""))   #use read_rds from readr next time
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
  p = data %>% filter(n.cov.variance > 0.025) %>% ggplot()+
    facet_grid(.~seqnames.x, scales="free_x")+
    geom_point(size = 0.3,aes(x = start.x,y = lr,color = n.cov.variance))+
    geom_segment(aes(x = start.y,xend = end.y,y = seg_lr,yend = seg_lr),colour = 'red')+
    theme_linedraw()
  
  return(p)
}

patient_num = 86
p4 = lr_plot(patient_num)

grid.arrange(p1,p2,p3,p4,nrow = 4)

getCovPlot <- function(stats, data, seg=FALSE, col=FALSE) {
  cov.plt <- with(stats, {
    plt <- ggplot(data)
    if (seg) {
      plt <- plt + geom_segment(size=2)
    } else {
      plt <- plt + geom_point(size=0.5)
    }
    if (col) {
      plt <- plt +
        aes(x = start, y = lr, xend = end, yend = lr, color = factor(seg %% 3)) +
        geom_hline(aes(yintercept=lr), data=Clr, size=1) +
        scale_color_manual(guide=FALSE, values=c("red", "black", "blue"))
    } else {
      plt <- plt +
        aes(x = start, y = lr, xend = end, yend = lr) +
        geom_hline(aes(yintercept=lr, color=factor(C)), data=Clr, size=1) +
        scale_color_manual(guide=FALSE, values=STRING_COL)
    }
    plt <- plt +
      facet_grid(.~seqnames, scales="free_x") +
      coord_cartesian(ylim=c(ymin, ymax)) +
      scale_y_continuous(breaks=Clr$lr, labels=Clr$C) +
      ylab("Absolute Copy Number") +
      theme_pubclean(base_size=14) +
      theme(
        panel.spacing = unit(0, "lines"),
        axis.title.x=element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank()
      )
    return(plt)
  })
  return(cov.plt)
}