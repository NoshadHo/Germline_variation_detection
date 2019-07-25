#implementing pool

#we have to do importPoolData manually:
# .importPoolData <- function(cnv.fns, opts) {
  # cnvs <- lapply(cnv.fns, function(fn) {
  #   cnv <- readRDS(fn)
  # })
#   cov <- do.call(cbind, lapply(cnvs, function(cnv) cnv$tile$n.cov))
#   sex <- sapply(cnvs, function(cnv) .detect.sex(cnv$var, cnv$tile))
#   sex.chr <- as.logical(as.character(seqnames(cnvs[[1]]$tile)) %in% c("chrX", "chrY"))
#   target <- cnvs[[1]]$tile$target
#   pd <- list(cov=cov, sex=sex, sex.chr=sex.chr, target=target)
#   return(pd)
# }
files = list.files()

  imported_pool = foreach (file_num = 1:227) %dopar% {
    cnv = readRDS(paste("./",files[file_num],"/",files[file_num],".rds",sep = ""))
    cov = cnv$tile$n.cov
    sex = detect.sex(cnv$var,cnv$tile)
    sex.chr = as.character(seqnames(cnv$tile)) %in% c("chrX", "chrY")
    target = cnv$tile$target
    
    cat(as.character(file_num),file="progress.txt",sep="\n",append=TRUE)
    return(list(cov,sex,sex.chr,target))
  }
  cov = list()
  sex = character()
  sex.chr = list()
  target = list()
  
  for (i in 1:length(imported_pool)){
    l = imported_pool[[i]]
    cov[[i]] = l[[1]]  
    sex[i] = l[[2]]
    sex.chr[[i]] = l[[3]]
    target[[i]] = l[[4]]
  }
  cov = as.matrix(matrix(unlist(cov),ncol = length(cov), byrow = FALSE))
  colnames(cov) = files[1:227]
  
  sex.chr = sex.chr[[1]]
  #sex.chr = as.data.frame(matrix(unlist(sex.chr),ncol = length(sex.chr), byrow = FALSE))
  #colnames(sex.chr) = files[1:227]
  
  target = target[[1]]
  #target = as.data.frame(matrix(unlist(target),ncol = length(target), byrow = FALSE))
  #colnames(target) = files[1:227]

pd <- list(cov=cov, sex=as.list(sex), sex.chr=sex.chr, target=target)

pool <- .createPool(pd, opts)
cnv$tile$lr = .poolPcaDenoise(cnv, pool, opts)
###########
.createPool <- function(pd, opts) {
  filter <- .poolFilter(pd$cov, pd$sex, pd$target, opts$pool.lo.cov, opts$pool.hi.cov, opts$pool.hi.zero)
  cov <- .outlierMask(pd$cov, pd$sex, filter, opts$pool.sd.out)
  lr <- .medianLogRatio(cov, pd$target, filter, pd$sex, pd$sex.chr)
  mc <- .medianCoverage(cov, pd$sex, pd$sex.chr)
  ##
  if (opts$pool.method=="pca") {
    p <- svd(lr[filter,])$u
  } else if (opts$pool.method=="ica") {
    p <- fastICA(lr[filter,], opts$pool.n.comp, method="C")$S
  } else {
    p <- NULL
  }
  pool <- list(method=opts$pool.method, sex=pd$sex, target=pd$target, sex.chr=pd$sex.chr,
               cov=cov, cov.med=mc, projection=p, filter=filter)
  return(pool)
}

###########

algo_performance_df = as.data.frame(t(do.call(cbind,algo_performance)))
colnames(algo_performance_df) = c("algo0_percent","algo1_percent","algo2_percent","algo3_percent","algo4_percent",
                                  "algo0_dist","algo1_dist","algo2_dist","algo3_dist","algo4_dist")

