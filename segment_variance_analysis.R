#adding segment to tile
seg_tiles <- queryHits(findOverlaps(cnv_old$tile,cnv_old$seg))
seg_ids <- subjectHits(findOverlaps(cnv_old$tile,cnv_old$seg))
cnv_old$tile$seg <- seg_ids

seg_tiles <- queryHits(findOverlaps(cnv_new$tile,cnv_new$seg))
seg_ids <- subjectHits(findOverlaps(cnv_new$tile,cnv_new$seg))
cnv_new$tile$seg <- seg_ids


#calculate the sd around each segment
old_sd = as.data.frame(cnv_old$tile) %>% group_by(seg) %>% summarise(sd = sd(lr,na.rm = T))
new_sd = as.data.frame(cnv_new$tile) %>% group_by(seg) %>% summarise(sd = sd(lr,na.rm = T))
mean(new_sd$sd,na.rm = T)
mean(old_sd$sd,na.rm = T)

#for all the samples:
agi.v4.files = read_delim(delim = "\n", file = "/data/cohorts/agi.v4.samples.txt", col_names = F)
olds = agi.v4.files$X1
news = paste0("/data/WXS/rerun1", substr(olds, start = 18, stop = nchar(olds)))

list = foreach(i = 1:length(olds)) %dopar% {
  cnv_old = read_rds(olds[i])
  cnv_new = read_rds(news[i])
  
  seg_tiles <- queryHits(findOverlaps(cnv_old$tile,cnv_old$seg))
  seg_ids <- subjectHits(findOverlaps(cnv_old$tile,cnv_old$seg))
  cnv_old$tile$seg <- seg_ids
  
  seg_tiles <- queryHits(findOverlaps(cnv_new$tile,cnv_new$seg))
  seg_ids <- subjectHits(findOverlaps(cnv_new$tile,cnv_new$seg))
  cnv_new$tile$seg <- seg_ids
  
  old_sd = as.data.frame(cnv_old$tile) %>% group_by(seg) %>% summarise(sd = sd(lr,na.rm = T))
  new_sd = as.data.frame(cnv_new$tile) %>% group_by(seg) %>% summarise(sd = sd(lr,na.rm = T))
  new = median(new_sd$sd,na.rm = T)
  old = median(old_sd$sd,na.rm = T)
  cat(as.character(i),file="/data/progress.txt",sep="\n",append=TRUE)
  return(c(old,new))
}
seg_medians = do.call(rbind,list)
colnames(seg_means) = c("old","new")
# these id's, sd has been increased
# [1]  49 118 127 130 133 134 136 159 160 179 185 191 193 200 207 213 216 217 218 219 220 221 222 223 224 225 226 227 228 230 231 239 240 241 244 246 249 250 252 253
# [41] 254 262 264 265 266 267 268 269 271 272 273 274 277 279 280 284 285 287 291 297 298 300 313 317 336 337 339 341 342 363

#Box plot
box_data = as.data.frame(seg_means) %>% gather()
box = box_data %>% ggplot() + geom_boxplot(aes(x = key, y = value), fill = c("#E69F00", "#56B4E9"))+scale_fill_brewer(palette="Dark2")


## PC_NUMBERS EFFECT -----------------------------------------
# load the pool data
pd = read_rds("/data/pools/agi.v4.pca.V2.pd.rds")
pool = read_rds("/data/pools/agi.v4.pca.V2.pool.rds")
PC_SIZES = seq(5,300,5)
opts = OPTS #(set opts accordingly)
opts$pool.hi.nvar = 10

#copying createpool code, so we can run it faster by some modifications
filter.female <- pool$female$filter
cov.female <- .outlierMask(pd$female$cov, filter.female, opts$pool.sd.out)
lr.female <- .medianLogRatio(cov.female, pd$target, filter.female, pd$sex.chr)
mc.female <- pool$female$cov.med

## male
filter.male <- .poolFilter(pd$male$cov, pd$male$nvar, pd$target, opts$pool.lo.cov, opts$pool.hi.cov, opts$pool.hi.zero, opts$pool.hi.nvar, opts$lr.method)
cov.male <- .outlierMask(pd$male$cov, filter.male, opts$pool.sd.out)
lr.male <- .medianLogRatio(cov.male, pd$target, filter.male, pd$sex.chr)
mc.male <- .medianCoverage(cov.male, pd$sex.chr)

#for all the samples:
agi.v4.files = read_delim(delim = "\n", file = "/data/cohorts/agi.v4.samples.txt", col_names = F)
olds = agi.v4.files$X1
news = paste0("/data/WXS/rerun1", substr(olds, start = 18, stop = nchar(olds)))

#modified segmentation
addJointSegmentModif <- function(cnv, opts) {
  tile = cnv$tile
  var = cnv$var
  ## add mirrored BAF to tiles
  tile <- .addBafTile(tile, var, opts)
  ## segment each arm
  arms <- split(tile, tile$arm)
  tmp <- lapply(arms,.jointSegArm,opts)
  ## compute genome-wide unique indexes
  tmp <- paste(tile$arm, unlist(tmp))
  tmp <- as.integer(factor(tmp, levels=unique(tmp)))
  ## create segmentation ranges
  cnv$seg <- unname(unlist(range(split(tile, tmp))))
  return(cnv)
}
#modified pool
modifyPool <- function(pool,opts) { #not an indipendent function, codes before it must be ran before running ti
  ## projection
  if (opts$pool.method=="pca") {
    p.female <- svd(lr.female[filter.female & pd$target,])$u
    p.male <- svd(lr.male[filter.male & pd$target,])$u
  } else if (opts$pool.method=="ica") {
    p.female <- fastICA(lr.female[filter.female & pd$target,], opts$pool.n.comp, method="C")$S
    p.male <- fastICA(lr.male[filter.male & pd$target,], opts$pool.n.comp, method="C")$S
  } else {
    p.female <- NULL
    p.male <- NULL
  }
  pool$female$projection <- p.female
  pool$male$projection <- p.male
  
  return(pool)
}

reAdjustCNV <- function(cnv,pool,opts) {
  cnv <- addProp(cnv, NULL, opts)
  ## log-ratio
  if (opts$lr.method=="nvar") {
    cnv <- addNvarLogRatio(cnv, pool, opts)
  } else if (opts$lr.method=="pool") {
    cnv <- addPoolLogRatio(cnv, pool, opts)
  } else if (opts$lr.method=="pair") {
    cnv <- addPairLogRatio(cnv, opts)
  } else if (opts$lr.method=="mean") {
    cnv <- addMeanLogRatio(cnv, opts)
  }
  cnv <- addScaleLogRatio(cnv, opts)
  cnv <- addGcLogRatio(cnv, opts)
  cnv <- addSmoothLogRatio(cnv, opts)
  
  ## segmentation
  registerDoParallel(cores = 10)
  cnv <- addJointSegmentModif(cnv, opts)
  stopImplicitCluster()
  return(cnv)
}

calcSd <- function(pool,opts) {
  list = foreach(i = 1:length(olds)) %dopar% {
    cnv_old = read_rds(olds[i])
    cnv_new = read_rds(news[i])
    # cnv_new = reAdjustCNV(cnv_old,pool,opts)
    
    seg_tiles <- queryHits(findOverlaps(cnv_old$tile,cnv_old$seg))
    seg_ids <- subjectHits(findOverlaps(cnv_old$tile,cnv_old$seg))
    cnv_old$tile$seg <- seg_ids
    
    seg_tiles <- queryHits(findOverlaps(cnv_new$tile,cnv_new$seg))
    seg_ids <- subjectHits(findOverlaps(cnv_new$tile,cnv_new$seg))
    cnv_new$tile$seg <- seg_ids
    
    old_sd = as.data.frame(cnv_old$tile) %>% group_by(seg) %>% summarise(sd = sd(lr,na.rm = T))
    new_sd = as.data.frame(cnv_new$tile) %>% group_by(seg) %>% summarise(sd = sd(lr,na.rm = T))
    new = mean(new_sd$sd,na.rm = T)
    old = mean(old_sd$sd,na.rm = T)
    return(c(old,new))
  }
  seg_means = do.call(rbind,list)
  colnames(seg_means) = c("old","new")
  
  return(seg_means)
}

run <- function(){
  # res = foreach(i = 1:3) %dopar% {
  #   opts$pool.n.comp = PC_SIZES[i]
  #   pool = modifyPool(pool,opts)
  #   seg_means = calcSd(pool,opts)
  #   cat(as.character(i),file="/data/progress.txt",sep="\n",append=TRUE)
  #   # return(c(median(seg_means[,1],na.rm = T),median(seg_means[,2],na.rm = T)))
  #   return(c(seg_means[,1],seg_means[,2]))
  # }
  res = lapply(X = 1:length(PC_SIZES), FUN = function(X){
    opts$pool.n.comp = PC_SIZES[X]
    pool = modifyPool(pool,opts)
    seg_means = calcSd(pool,opts)
    cat(as.character(i),file="/data/progress.txt",sep="\n",append=TRUE)
    # return(c(median(seg_means[,1],na.rm = T),median(seg_means[,2],na.rm = T)))
    return(c(seg_means[,1],seg_means[,2]))
  })
  #here seperate them into columns(each column have the old+new ones of each pc size)
  sds = do.call(cbind,res)
  #Then seperate them into two different data sets
  old_sds = sds[1:371,]
  colnames(old_sds) = paste(PC_SIZES," comp")
  new_sds = sds[372:742,]
  colnames(new_sds) = str_pad(as.character(PC_SIZES), 3, pad = "0")
  #gather
  new_sds_gather = gather(as.data.frame(new_sds))
  old_sds_gather = gather(as.data.frame(old_sds))
  #plot two boxplotos
  colnames(seg_means) = c("old","new")
  
  ##then plot a box plot for each PC_SIZE
  box = new_sds_gather %>% ggplot()+geom_boxplot(aes(x = key, y = value))
}