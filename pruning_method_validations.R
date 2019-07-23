#overRiding requiered functions-----------------
  .mergeBreakpoints <- function(seg0, sd.lr, sd.baf, arm, algo) {
    best1 <- integer()
    ## got at least one breakpoint
    if (length(seg0) > 0) {
      ## compute breakpoint stats
      beg0 <- c(1, seg0+1)
      end0 <- c(seg0, length(arm))
      len0 <- diff(c(0, end0))
      idx0 <- rep(seq_along(beg0), len0)
      lr0 <- split(arm$lr, idx0)
      baf0 <- split(arm$baf, idx0)
      var0 = split(arm$n.cov.variance, idx0)
      stat0 <- data.table(
        seg=seg0,
        lr.diff =sapply(2:length(beg0), function(i) .absMedDiff( lr0[[i]],  lr0[[i-1]])),
        baf.diff=sapply(2:length(beg0), function(i) .absMedDiff(baf0[[i]], baf0[[i-1]])),
        min.len =sapply(2:length(beg0), function(i) min(length(lr0[[i]]), length(lr0[[i-1]]))),
        max.lr.var = sapply(2:length(beg0), function(i) max(mean(var0[[i]]), mean(var0[[i-1]])))
      )
      stat0[is.na(lr.diff), lr.diff:=0]
      stat0[is.na(baf.diff), baf.diff:=0]
      stat0[,len.penalty:=.len.penalty(min.len)]
      ## pick weakest breakpoint to merge
      stat0 <- stat0[order(lr.diff, baf.diff)]
      
      if (algo == 1){
        stat0[,len.penalty:=.len.penalty(min.len-1)]
        seg1 <- stat0[(
          lr.diff <  sd.lr +  sd.lr * (len.penalty) &
            baf.diff < sd.baf + sd.baf * len.penalty
        ), seg]
      }else if(algo == 2){
        theta = 1
        seg1 <- stat0[(
          lr.diff <  (sd.lr +  sd.lr * len.penalty)*(1+theta*max.lr.var)&
            baf.diff < sd.baf + sd.baf * len.penalty
        ), seg]
      }else if(algo == 3){
        #adjust the weights
        a = 0.02
        b = 0.08
        stat0[max.lr.var < a, max.lr.var:= 0]
        stat0[max.lr.var > b, max.lr.var:= 1]
        stat0[max.lr.var < b & max.lr.var > a, max.lr.var:= (max.lr.var-a)/(b-a)]
        
        seg1 <- stat0[(
          lr.diff <  (sd.lr +  sd.lr * len.penalty + sd.lr*max.lr.var) &
            baf.diff < sd.baf + sd.baf * len.penalty
        ), seg]
      }else if (algo == 4){
        a = 0.02
        b = 0.08
        stat0[max.lr.var < a, max.lr.var:= 0]
        stat0[max.lr.var > b, max.lr.var:= 1]
        stat0[max.lr.var < b & max.lr.var > a, max.lr.var:= (exp((max.lr.var-a)/(b-a))-1)/exp(1)]
        
        seg1 <- stat0[(
          lr.diff <  (sd.lr +  sd.lr * len.penalty + sd.lr*max.lr.var) &
            baf.diff < sd.baf + sd.baf * len.penalty
        ), seg]
      }else{ #this is the default
        seg1 <- stat0[(
          lr.diff <  sd.lr +  sd.lr * len.penalty &
            baf.diff < sd.baf + sd.baf * len.penalty
        ), seg]
      }
      
      best1 <- head(seg1, 1)
    }
    return(best1)
  }

  .jointSegArm <- function(arm, sd.lr, sd.baf, method, tile.width, len.min, cbs.lr, cbs.baf, rbs.selection, sd.prune, len.prune,algo) {
    ## initial segmentation
    pruned_BP_count = 0
    seg1 <- .runJointSeg(arm, method, tile.width, len.min, cbs.lr, cbs.baf, rbs.selection, FALSE)
    ## iteratively remove spurious breakpoints
    if (sd.prune) {
      repeat({
        best1 <- .mergeBreakpoints(seg1, sd.lr, sd.baf, arm, algo)
        seg1 <- setdiff(seg1, best1)
        if (length(best1) == 0) {break}
        pruned_BP_count = pruned_BP_count + 1
      })
    }
    ## skip short segments
    if (len.prune) {
      end1 <- c(seg1, length(arm))
      len1 <- diff(c(0, end1))
      seg1 <- head(end1[len1>=len.min], -1)
    }
    ## append segmentation
    bpt1 <- c(1, seg1 + 1)
    len1 <- diff(c(bpt1, length(arm)+1))
    idx1 <- rep(seq_along(bpt1), len1)
    arm$seg <- idx1
    arm$pruned_BP_count = pruned_BP_count
    return(arm)
  }
  
  
  .addJointSeg <- function(gt, ...) {
    #gt$rep = 0
    gt <- sort(unname(unlist(endoapply(split(gt, gt$arm), .jointSegArm, ...))))
    ## provide globally unique ids
    tmp <- paste(gt$arm, gt$seg)
    gt$seg <- as.integer(factor(tmp, levels=unique(tmp)))
    return(gt)
  }
  
  gt <- sort(unname(unlist(endoapply(split(gt, gt$arm), .jointSegArm, sd.lr, sd.baf, method, tile.width, len.min, cbs.lr, cbs.baf, rbs.selection, sd.prune, len.prune,algo))))
  
  addJointSegment <- function(cnv, opts, algo) {
    ## estimated standard-deviation
    hq.lr <- cnv$tile[cnv$tile$unmasked > 0.25 & cnv$tile$blacklist==0.0]$lr
    sd.lr <- estimateSd(hq.lr) * opts$seg.sd.lr.penalty
    sd.baf <- estimateSd(cnv$tile$baf) * opts$seg.sd.baf.penalty
    sd.baf <- ifelse(!is.finite(sd.baf), .Machine$double.eps, sd.baf)
    ## create segmentation
    cnv$tile <- .addJointSeg(cnv$tile, sd.lr, sd.baf, opts$seg.method,
                             opts$tile.width, opts$seg.len.min, opts$seg.cbs.lr,
                             opts$seg.cbs.baf, opts$seg.rbs.selection, opts$seg.sd.prune, opts$seg.len.prune,algo)    

    cnv$seg <- unname(unlist(range(split(cnv$tile, cnv$tile$seg))))
    ## seg -> var
    cnv$var$seg <- findOverlaps(cnv$var, cnv$seg, select="first", maxgap = opts$tile.shoulder-1)
    return(cnv)
  }

#run the functions and calculate the percentage of removed breakpoints
  pruned_percent = function(cnv,algo){
    cnv = addJointSegment(cnv,opts,algo)
    pruned_BP_count = sum((as.data.frame(cnv$tile) %>% group_by(arm) %>% dplyr::slice(1))$pruned_BP_count)
    pruned_BP_percent = pruned_BP_count/(pruned_BP_count + max(cnv$tile$seg))
    
    data = as.data.frame(cnv$tile) %>% select(seqnames,start,end,t.cov,n.cov,lr,n.cov.variance,seg)
    data_seg = as.data.frame(cnv$seg) %>% mutate(seg = row_number()) %>% select(-strand)
    data = data %>% left_join(data_seg,by = "seg")
    data = data %>% group_by(seg) %>% mutate(seg_lr = mean(lr,na.rm = T)) %>% ungroup()
    
    data = data %>% filter(n.cov.variance > 0.02) %>% rowwise() %>% 
      mutate(breakpoint.dist = min((start.x - start.y),(end.y - end.x)))%>% ungroup()
    
    mean_seg_dist = (data %>% group_by(seg) %>% summarise(seg_dist = mean(breakpoint.dist)) %>% ungroup() %>% summarise(dist = mean(seg_dist)))$dist
    
    
    return(c(pruned_BP_percent,mean_seg_dist))
  }
  (as.data.frame(cnv$tile) %>% group_by(arm) %>% dplyr::slice(1))$pruned_BP_count
  
  algo0_percent = numeric()
  algo1_percent = numeric()
  algo2_percent = numeric()
  algo3_percent = numeric()
  algo4_percent = numeric()
  
  algo0_dist = numeric()
  algo1_dist = numeric()
  algo2_dist = numeric()
  algo3_dist = numeric()
  algo4_dist = numeric()
  
  highly_variable_points = (variance_raw %>% filter(variance > 0.02))$tile
  for (file_num in 1:(length(files))){
    print(paste("[1.1] Start Reading:", file_num))
    time1 = system.time({
      file = readRDS(paste("./",files[file_num],"/",files[file_num],".rds",sep = ""))   #use read_rds from readr next time
    })
    print(paste("[1.2] Reading time:", as.numeric(time1)[3]))
    print(paste("[2.1] Start Processing:", file_num))
    time2 = system.time({
    file$tile$n.cov.variance = variance_raw$variance
    
    algo0 = pruned_percent(file,0)
    algo1 = pruned_percent(file,1)
    algo2 = pruned_percent(file,2)
    algo3 = pruned_percent(file,3)
    algo4 = pruned_percent(file,4)
    
    algo0_percent[file_num] = algo0[1]
    algo1_percent[file_num] = algo1[1]
    algo2_percent[file_num] = algo2[1]
    algo3_percent[file_num] = algo3[1]
    algo4_percent[file_num] = algo4[1]
    
    algo0_dist[file_num] = algo0[2]
    algo1_dist[file_num] = algo1[2]
    algo2_dist[file_num] = algo2[2]
    algo3_dist[file_num] = algo3[2]
    algo4_dist[file_num] = algo4[2]
    
    rm(algo0,algo1,algo2,algo3,algo4)
    })
    print(paste("[2.2] Processing time:", as.numeric(time2)[3],'\n'))
  }
  
  
  

  algo_performance = foreach (file_num = 1:(length(files))) %dopar% {
    #print(paste("[1.1] Start Reading:", file_num))
    time1 = system.time({
      file = readRDS(paste("./",files[file_num],"/",files[file_num],".rds",sep = ""))   #use read_rds from readr next time
    })
    #print(paste("[1.2] Reading time:", as.numeric(time1)[3]))
    #print(paste("[2.1] Start Processing:", file_num))
    time2 = system.time({
      file$tile$n.cov.variance = variance_raw$variance
      
      algo0 = pruned_percent(file,0)
      algo1 = pruned_percent(file,1)
      algo2 = pruned_percent(file,2)
      algo3 = pruned_percent(file,3)
      algo4 = pruned_percent(file,4)
      
      algo0_percent = algo0[1]
      algo1_percent = algo1[1]
      algo2_percent = algo2[1]
      algo3_percent = algo3[1]
      algo4_percent = algo4[1]
      
      algo0_dist = algo0[2]
      algo1_dist = algo1[2]
      algo2_dist = algo2[2]
      algo3_dist = algo3[2]
      algo4_dist = algo4[2]
      
      cat(as.character(file_num),file="progress.txt",sep="\n",append=TRUE)
      
      return(c(algo0_percent,algo1_percent,algo2_percent,algo3_percent,algo4_percent,
               algo0_dist,algo1_dist,algo2_dist,algo3_dist,algo4_dist))
    })
    #print(paste("[2.2] Processing time:", as.numeric(time2)[3],'\n'))
  }
  algo_performance = as.data.frame(t(do.call(cbind,algo_performance)))
  