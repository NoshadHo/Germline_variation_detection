#vcf
.preFilterVar <- function(var, opts) {
  var <- var[
    lengths(fixed(var)$ALT) == 1 &
      isSNV(var) &
      !is.na(geno(var)$DP[,1]) &
      !is.na(qual(var))
    ]
  var <- var[
    geno(var)$DP[,1] > opts$baf.min.het.dp
    ]
  return(var)
}

.maskVar <- function(var, opts) {
  ## masked regions 1000G
  pilot.fn <- system.file("extdata/hg38/1000G-pilot-unmask.bed.gz", package="cnvex")
  pilot <- .robust.import(pilot.fn, opts$genome)
  strict.fn <- system.file("extdata/hg38/1000G-strict-unmask.bed.gz", package="cnvex")
  strict <- .robust.import(strict.fn, opts$genome)
  tmp <- DataFrame(Number=c("1","1"),
                   Type=c("Integer", "Integer"),
                   Description=c("loose mask", "strict mask"),
                   row.names=c("mask.loose", "mask.strict"))
  info(header(var)) <- rbind(info(header(var)), tmp)
  info(var)$mask.loose <- as.integer(!(var %over% pilot))
  info(var)$mask.strict <- as.integer(!(var %over% strict))
  return(var)
}

.importVcf <- function(vcf, param, opts) {
  var <- .robust.import.vcf(vcf, param, opts$genome)
  var <- .preFilterVar(var, opts)
  var <- .maskVar(var, opts)
  return(var)
}
#tile
.annotateTiles <- function(genome.tile, genome) {
  
  bl1.fn <- system.file("extdata/hg38/sv-blacklist-10x-ucsc.bed.gz", package="cnvex")
  bl2.fn <- system.file("extdata/hg38/blacklist.bed.gz", package="cnvex")
  cyto.fn <- system.file("extdata/hg38/cytoband.bed.gz", package="cnvex")
  strict.fn <- system.file("extdata/hg38/1000G-strict-unmask.bed.gz", package="cnvex")
  
  ## blacklist
  bl.1 <- .robust.import(bl1.fn, genome)
  bl.2 <- .robust.import(bl2.fn, genome)
  bl <- reduce(c(granges(bl.1), granges(bl.2)))
  tmp <- findOverlaps(genome.tile, bl)
  tmp <- data.table(
    tile=queryHits(tmp),
    blacklist=width(pintersect(genome.tile[queryHits(tmp)], bl[subjectHits(tmp)]))
  )
  setkey(tmp, tile)
  tmp <- tmp[J(seq_along(genome.tile))]
  tmp[is.na(blacklist), blacklist:=0]
  tmp <- tmp[,.(blacklist=sum(blacklist)),by=tile]
  genome.tile$blacklist <- tmp$blacklist / width(genome.tile)
  
  ## masking
  strict <- .robust.import(strict.fn, genome)
  tmp <- findOverlaps(genome.tile, strict)
  tmp <- data.table(
    tile=queryHits(tmp),
    unmasked=width(pintersect(genome.tile[queryHits(tmp)], strict[subjectHits(tmp)]))
  )
  setkey(tmp, tile)
  tmp <- tmp[J(seq_along(genome.tile))]
  tmp[is.na(unmasked), unmasked:=0]
  tmp <- tmp[,.(unmasked=sum(unmasked)),by=tile]
  genome.tile$unmasked <- tmp$unmasked / width(genome.tile)
  
  ## GC content
  tmp <- getSeq(BSgenome.Hsapiens.UCSC.hg38, genome.tile)
  genome.tile$gc <- letterFrequency(tmp, "GC", as.prob=TRUE)[,1]
  
  ## has assembly gaps
  genome.tile$gap <- letterFrequency(tmp, "N", as.prob=TRUE)[,1]
  
  ## cytobands
  cyto <- .robust.import(cyto.fn, genome)
  tmp <- findOverlaps(genome.tile, cyto, select="first")
  genome.tile$cytoband <- cyto[tmp]$name
  
  ## order arms
  genome.tile <- sort(genome.tile)
  genome.tile$arm <- paste0(seqnames(genome.tile), str_sub(genome.tile$cytoband, 1, 1))
  genome.tile$arm <- factor(genome.tile$arm, unique(genome.tile$arm), ordered=TRUE)
  
  ## define HQ tile
  genome.tile$hq <-
    genome.tile$gap < 0.005 &
    genome.tile$unmasked > 0.25 &
    genome.tile$blacklist < 0.001
  return(genome.tile)
}

getGenomeTiles <- function(tgt.fn, opts) {
  ## tile genome
  genome.tile <- .tilegenome(opts$genome, opts$tile.width)
  genome.tile$target <- TRUE
  genome.tile <- .annotateTiles(genome.tile, opts$genome)
  return(genome.tile)
}

getTargetTiles <- function(tgt.fn, opts) {
  ## tile targets
  seqi <- .seqinfo(opts$genome)
  tgt <- .robust.import(tgt.fn, opts$genome)
  tgt$name <- NULL
  tgt$score <- NULL
  tgt$target <- TRUE
  gap <- keepSeqlevels(gaps(tgt), seqlevels(seqi), pruning.mode="coarse")
  gap <- gap[strand(gap)=="*"]
  gap <- gap[width(gap)>2*opts$tile.shoulder+opts$tile.min.gap]
  gap <- gap-opts$tile.shoulder
  gap$target <- FALSE
  target.tile <- sort(c(tgt, gap))
  target.tile <- .annotateTiles(target.tile, opts$genome)
  return(target.tile)
}

#seg
.len.penalty <- function(x) {
  1/x
}

.runJointSeg <- function(arm, method, tile.width, len.min, cbs.lr, cbs.baf, rbs.selection, only.target) {
  arm.lr <- arm$lr
  arm.lr.weigth = arm$lr.weight
  arm.baf <- arm$baf
  arm.baf.weight <- arm$baf.weight
  arm.baf.weight[is.na(arm.baf.weight)] <- 1.0
  
  if (only.target) {
    arm.lr[!arm$target] <- NA_real_
    arm.baf[!arm$target] <- NA_real_
  }
  ## make sure we have enough points to segment
  seg0 <- integer()
  if (method %in% c("RBS", "DynamicProgramming")) {
    n.points <- sum(!is.na(arm.lr) & !is.na(arm.baf))
    if (n.points >= 2*len.min) {
      if (!is.na(tile.width)) {
        tile.per.mb <- 1e6 / tile.width
        opt.K <- ceiling(n.points / tile.per.mb)
      } else {
        opt.K <- max(2, floor(n.points / (2*len.min)))
      }
      seg0 <- suppressWarnings(sort(unique(jointSeg(
        cbind(arm.lr, arm.baf), method=method, modelSelectionMethod=rbs.selection, K=opt.K)$bestBkp)))
    }
  } else if (method=="CBS") {
    n.points.lr <- sum(!is.na(arm.lr))
    if (n.points.lr >= 2*len.min) {
      seg0.lr <- .runCBS(arm.lr,arm.lr.weigth, args=cbs.lr)
    } else {
      seg0.lr <- c()
    }
    n.points.baf <- sum(!is.na(arm.baf))
    if (n.points.baf >= 2*len.min) {
      seg0.baf <- .runCBS(arm.baf, arm.baf.weight, args=cbs.baf)
    } else {
      seg0.baf <- c()
    }
    seg0 <- unique(sort(c(seg0.lr, seg0.baf)))
  } else {
    stop("Joint segmentation method not supported.")
  }
  return(seg0)
}

.mergeBreakpoints <- function(seg0, sd.lr, sd.baf, arm) {
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
    stat0 <- data.table(
      seg=seg0,
      lr.diff =sapply(2:length(beg0), function(i) .absMedDiff( lr0[[i]],  lr0[[i-1]])),
      baf.diff=sapply(2:length(beg0), function(i) .absMedDiff(baf0[[i]], baf0[[i-1]])),
      min.len =sapply(2:length(beg0), function(i) min(length(lr0[[i]]), length(lr0[[i-1]])))
    )
    stat0[is.na(lr.diff), lr.diff:=0]
    stat0[is.na(baf.diff), baf.diff:=0]
    stat0[,len.penalty:=.len.penalty(min.len)]
    ## pick weakest breakpoint to merge
    stat0 <- stat0[order(lr.diff, baf.diff)]
    seg1 <- stat0[(
      lr.diff <  sd.lr +  sd.lr * len.penalty &
        baf.diff < sd.baf + sd.baf * len.penalty
    ), seg]
    best1 <- head(seg1, 1)
  }
  return(best1)
}

.jointSegArm <- function(arm, sd.lr, sd.baf, method, tile.width, len.min, cbs.lr, cbs.baf, rbs.selection, sd.prune, len.prune) {
  ## initial segmentation
  seg1 <- .runJointSeg(arm, method, tile.width, len.min, cbs.lr, cbs.baf, rbs.selection, FALSE)
  ## iteratively remove spurious breakpoints
  if (sd.prune) {
    repeat({
      best1 <- .mergeBreakpoints(seg1, sd.lr, sd.baf, arm)
      seg1 <- setdiff(seg1, best1)
      if (length(best1) == 0) {break}
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
  return(arm)
}

.addJointSeg <- function(gt, ...) {
  gt <- sort(unname(unlist(endoapply(split(gt, gt$arm), .jointSegArm, ...))))
  ## provide globally unique ids
  tmp <- paste(gt$arm, gt$seg)
  gt$seg <- as.integer(factor(tmp, levels=unique(tmp)))
  return(gt)
}

#rule
.detect.sex <- function(var, tile) {
  sel <- var$n.GT %in% c("0/1", "1/0") & var$n.DP>16
  chrx.snp <- length(var[sel & seqnames(var)=="chrX"])
  chr2.snp <- length(var[sel & seqnames(var)=="chr2"])
  sex.snp <-  if (chrx.snp / chr2.snp > 0.25) "female" else "male"
  return(sex.snp)
}

.detect.offset <- function(cnv, opts) {
  tmp <- as.data.table(mcols(cnv$tile)[,c("seg", "lr", "baf")])
  tmp <- tmp[,.(.N, lr=mean(lr, na.rm=TRUE), baf=mean(baf, na.rm=TRUE)), by=seg]
  tmp <- tmp[is.finite(lr) & is.finite(baf) & baf>0.38]
  tmp <- tmp[order(lr)][1:(nrow(tmp)%/%2)]
  tmp <- tmp[order(-N)][1:(nrow(tmp)%/%1.5)]
  offset <- median(tmp$lr)
  return(offset)
}

#pool
.coverageFilterSel <- function(xs, lo, hi) {
  uxs <- rowMeans(xs)
  muxs <- median(uxs)
  f <- lo * muxs < uxs & uxs < hi * muxs
  return(f)
}

.zeroFilterSel <- function(xs, hi) {
  f <- rowMeans(xs==0) < hi & rowMads(xs)!=0
  return(f)
}

.coverageFilter <- function(x, s, lo, hi) {
  f <- logical(length(s))
  f[ s] <- .coverageFilterSel(x[ s,], lo, hi)
  f[!s] <- .coverageFilterSel(x[!s,], lo, hi)
  return(f)
}

.zeroFilter <- function(x, s, hi) {
  f <- logical(length(s))
  f[ s] <- .zeroFilterSel(x[ s,], hi)
  f[!s] <- .zeroFilterSel(x[!s,], hi)
  return(f)
}

.poolFilter <- function(cov, sex, target, lo.cov, hi.cov, hi.zero) {
  cfm <- .coverageFilter(cov[,sex==  "male"], target, lo.cov, hi.cov)
  cff <- .coverageFilter(cov[,sex=="female"], target, lo.cov, hi.cov)
  zfm <- .zeroFilter(cov[,sex==  "male"], target, hi.zero)
  zff <- .zeroFilter(cov[,sex=="female"], target, hi.zero)
  f <- (cfm | cff) & (zfm | zff)
  return(f)
}

.outlierMaskSel <- function(xs, sd) {
  ## identify outlier tiles
  xs.med <- rowMedians(xs)
  xs.mad <- rowMads(xs)
  xs.out <- (xs - xs.med) / xs.mad
  ## lo-threshold
  xs[xs.out < -sd] <- NA_real_
  xs.min <- rowMins(xs, na.rm = TRUE)
  xs <- apply(xs, 2, function(col) {
    col[is.na(col)] <- xs.min[is.na(col)]
    col
  })
  ## hi-threshold
  xs[xs.out > sd] <- NA_real_
  xs.max <- rowMaxs(xs, na.rm = TRUE)
  xs <- apply(xs, 2, function(col) {
    col[is.na(col)] <- xs.max[is.na(col)]
    col
  })
  return(xs)
}

.outlierMask <- function(cov, sex, filter, sd) {
  cov[ filter,sex==  "male"] <- .outlierMaskSel(cov[filter,sex==  "male"], sd)
  cov[ filter,sex=="female"] <- .outlierMaskSel(cov[filter,sex=="female"], sd)
  cov[!filter,] <- NA_real_
  return(cov)
}

.fixLogRatioBounds <- function(xsr) {
  xsrm <- max(min(xsr[is.finite(xsr)]), -4)
  xsrx <- min(max(xsr[is.finite(xsr)]),  4)
  xsr[is.nan(xsr)] <- 0    ## 0 / 0
  xsr[xsr < xsrm] <-  xsrm ## 0 / x
  xsr[xsr > xsrx] <-  xsrx ## x / 0
  return(xsr)
}

.medianLogRatioSel <- function(xs, sex, sex.chr) {
  xsr <- xs
  ## autosomes
  xsa <- xs[!sex.chr,]
  xsam <- rowMedians(xsa)
  xsr[!sex.chr,] <- log2(xsa/xsam)
  ## sex chromosomes
  xsm <- xs[sex.chr, sex=="male"]
  xsmm <- rowMedians(xsm)
  xsr[sex.chr,sex=="male"] <- log2(xsm/xsmm)
  xsf <- xs[sex.chr, sex=="female"]
  xsfm <- rowMedians(xsf)
  xsr[sex.chr,sex=="female"] <- log2(xsf/xsfm)
  ## fix bounds
  xsr <- .fixLogRatioBounds(xsr)
  return(xsr)
}

.medianLogRatio <- function(x, s, f, sex, sex.chr) {
  x[( s & f),] <- .medianLogRatioSel(x[( s & f),], sex, sex.chr[( s & f)])
  x[(!s & f),] <- .medianLogRatioSel(x[(!s & f),], sex, sex.chr[(!s & f)])
  x[f,] <- t(t(x[f,]) - colMedians(x[f,]))
  return(x)
}

.medianCoverage <- function(cov, sex, sex.chr) {
  m <- rowMedians(cov)
  mm <- rowMedians(cov[,sex==  "male"])
  mf <- rowMedians(cov[,sex=="female"])
  tile.median <- cbind(
    male=ifelse(sex.chr, mm, m),
    female=ifelse(sex.chr, mf, m)
  )
  return(tile.median)
}

.importPoolData <- function(cnv.fns, opts) {
  cnvs <- lapply(cnv.fns, function(fn) {
    cnv <- readRDS(fn)
  })
  cov <- do.call(cbind, lapply(cnvs, function(cnv) cnv$tile$n.cov))
  sex <- sapply(cnvs, function(cnv) .detect.sex(cnv$var, cnv$tile))
  sex.chr <- as.logical(as.character(seqnames(cnvs[[1]]$tile)) %in% c("chrX", "chrY"))
  target <- cnvs[[1]]$tile$target
  pd <- list(cov=cov, sex=sex, sex.chr=sex.chr, target=target)
  return(pd)
}

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

.poolIcaDenoise <- function(cnv, pool, opts) {
  sex <- .detect.sex(cnv$var, cnv$tile)
  S <- pool$projection[,seq(1,opts$pool.n.comp)]
  lr <- log2(cnv$tile$t.cov / pool$cov.med[,sex])
  x <- .fixLogRatioBounds(lr[pool$filter])
  x <- x - as.vector(tcrossprod(x %*% t(ginv(S)), S))
  lr[ pool$filter] <- x
  lr[!pool$filter] <- NA_real_
  return(lr)
}

.poolPcaDenoise <- function(cnv, pool, opts) {
  sex <- .detect.sex(cnv$var, cnv$tile)
  P <- pool$projection[,seq(1,opts$pool.n.comp)]
  lr <- log2(cnv$tile$t.cov / pool$cov.med[,sex])
  x <- .fixLogRatioBounds(lr[pool$filter])
  x <- x - as.vector(tcrossprod(x %*% P, P))
  lr[ pool$filter] <- x
  lr[!pool$filter] <- NA_real_
  return(lr)
}

.poolKnnCoverage <- function(cnv, pool, opts) {
  sex <- .detect.sex(cnv$var, cnv$tile)
  pool.sex <- pool$cov[,pool$sex==sex]
  ## rank normals by variance in log2(tumor/normal)
  lr.pool.mx <- log2(cnv$tile$t.cov/pool.sex)
  lr.pool.mx <- ifelse(is.finite(lr.pool.mx), lr.pool.mx, NA_real_)
  lr.pool.sd <- apply(lr.pool.mx, 2, estimateSd)
  lr.pool.rk <- order(lr.pool.sd)
  ## greedy search for k-best normals to pool
  k <- which.min(sapply(seq_len(min(25, ncol(pool.sex))), function(i) {
    n <- lr.pool.rk[1:i]
    n.pool <- pool.sex[,n,drop=FALSE]
    n.cov <- .normCoverage(rowMedians(n.pool), cnv$tile, opts)
    lr.pool <- log2(cnv$tile$t.cov/n.cov)
    lr.pool <- ifelse(is.finite(lr.pool), lr.pool, NA_real_)
    sd.pool <- estimateSd(lr.pool)
  }))
  ## compute final log-ratio
  n <- lr.pool.rk[1:k]
  n.pool <- pool.sex[,n,drop=FALSE]
  n.cov <- .normCoverage(rowMedians(n.pool), cnv$tile, opts)
  return(n.cov)
}

#plot
plotGC <- function(cnv) {
  tmp <- as.data.table(mcols(cnv$tile))
  if (!all(is.na(tmp$n.cov))) {
    tmp[,lr.raw:=log2(t.cov/n.cov)]
  } else {
    m.cov <- rep(mean(cnv$tile$t.cov, 0.01), length(cnv$tile))
    tmp[,lr.raw:=log2(t.cov/m.cov)]
  }
  tmp <- melt(tmp[,.(gc, blacklist, target, lr.raw, lr, lr.off=lr.raw-lr)], id.vars=c("gc", "blacklist", "target"))
  tmp[,delta:=FALSE]
  tmp[variable=="lr.off", delta:=TRUE]
  tmp[variable=="lr.off", variable:="lr.raw"]
  tmp[,target.lab:=ifelse(target, "on-target", "off-target")]
  tmp[,lr.lab:=ifelse(variable=="lr.raw", "raw", "adjusted")]
  plt <- ggplot(tmp) + aes(x=gc, y=value, color=delta) +
    facet_grid(lr.lab~target.lab) +
    geom_point(alpha=0.05, size=0.1) +
    scale_color_manual(values=c("black", "red"), guide=FALSE) +
    scale_x_continuous(labels=scales::percent) +
    coord_cartesian(xlim=c(0.2, 0.8), ylim=c(-3, 3)) +
    ylab("log2(tumor/normal)") +
    xlab("GC [%]") +
    theme_pubr()
  return(plt)
}

plotRatio <- function(cnv, opts, sel.chr=NULL) {
  if (is.null(sel.chr)) {
    sel.chr <- paste("chr", c(1:22, "X", "Y"), sep="")
  }
  cov <- cnv$tile[cnv$tile$unmasked>0.25]
  cov <- cov[seqnames(cov) %in% sel.chr]
  
  cov.dt <- data.table(
    chr=as.character(seqnames(cov)),
    pos=floor((start(cov)+end(cov))/2),
    val=mcols(cov)[["lr"]],
    type="COV",
    seg=cov$seg,
    tgt=mcols(cov)[["target"]]
  )[!is.na(val)]
  cov.dt[,chr:=factor(chr, sel.chr, ordered=TRUE)]
  cov.dt <- cov.dt[-c(1,nrow(cov.dt))][order(-tgt)]
  ## COV
  cov.plt <- ggplot(cov.dt) +
    facet_grid(.~chr, scales="free_x") +
    aes(x=pos, y=val, color=factor(ifelse(!tgt,4,as.integer(seg%%3)))) +
    coord_cartesian(ylim=c(-3,3)) +
    scale_color_igv(guide=FALSE) +
    geom_point(size=0.75, alpha=0.5, shape=16) +
    theme_pubr() +
    ylab("log2(tumor/normal)") +
    theme(
      panel.spacing = unit(0, "lines"),
      axis.title.x=element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )
  
  return(cov.plt)
}

plotSegPoint <- function(cnv, opts, sel.chr=NULL) {
  if (is.null(sel.chr)) {
    sel.chr <- paste("chr", c(1:22, "X", "Y"), sep="")
  }
  pass <- .filterVar(cnv$tile, cnv$var, opts)
  snp <- cnv$var[pass]
  snp <- snp[seqnames(snp) %in% sel.chr]
  cov <- cnv$tile
  cov <- cov[seqnames(cov) %in% sel.chr]
  
  if (length(unique(snp$SOURCE))==1) {
    snp$SOURCE <- "genome"
  }
  
  snp.dt <- data.table(
    chr=as.character(seqnames(snp)),
    pos=start(snp),
    val=snp$t.AF,
    src=snp$SOURCE,
    seg=snp$seg,
    type="BAF"
  )[!is.na(val)]
  snp.dt[,chr:=factor(chr, sel.chr, ordered=TRUE)]
  
  cov.dt <- data.table(
    chr=as.character(seqnames(cov)),
    pos=floor((start(cov)+end(cov))/2),
    val=mcols(cov)[["lr"]],
    type="COV",
    seg=cov$seg,
    tgt=mcols(cov)[["target"]]
  )[!is.na(val)]
  cov.dt[,chr:=factor(chr, sel.chr, ordered=TRUE)]
  cov.dt <- cov.dt[-c(1,nrow(cov.dt))][order(-tgt)]
  ## fix
  snp.dt <- rbind(snp.dt, cov.dt[,.SD[1],by=chr][,.(chr, pos, val=0.5, src="genome", seg=0, type="BAF")])
  cov.dt <- rbind(cov.dt, snp.dt[,.SD[1],by=chr][,.(chr, pos, val=0.0, seg=0, type="COV", tgt=TRUE)])
  
  ## COV
  offset <- .detect.offset(cnv, opts)
  cov.plt <- ggplot(cov.dt) +
    facet_grid(.~chr, scales="free_x") +
    aes(x=pos, y=val, color=factor(ifelse(!tgt,4,as.integer(seg%%3)))) +
    coord_cartesian(ylim=c(-3,3)) +
    scale_color_igv(guide=FALSE) +
    geom_point(size=0.75, alpha=0.5, shape=16) +
    ## geom_hline(yintercept = 0, color="blue") +
    ## geom_hline(yintercept = offset, color="red") +
    theme_pubr() +
    ylab("log2(tumor/normal)") +
    theme(
      panel.spacing = unit(0, "lines"),
      axis.title.x=element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )
  
  ## SNP
  if (nrow(snp.dt)>1e5) {
    p.size=0.2
    p.alpha=0.5
  } else {
    p.size=0.5
    p.alpha=0.75
  }
  snp.plt <- ggplot(snp.dt) +
    facet_grid(.~chr, scales="free_x") +
    aes(x=pos, y=val, color=factor(as.integer(seg%%3))) +
    coord_cartesian(ylim=c(0,1)) +
    geom_point(data=snp.dt[src=="genome"], size=p.size, shape=16) +
    geom_point(data=snp.dt[src=="target"], color="black", size=p.size, shape=16) +
    scale_color_igv(guide=FALSE) +
    theme_pubr() +
    ylab("BAF") +
    theme(
      panel.spacing = unit(0, "lines"),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )
  capture.output({
    pdf(NULL)
    plt <- rbind(ggplotGrob(cov.plt), ggplotGrob(snp.plt), size = "last")
    dev.off()
  })
  return(plt)
}

plotGrid <- function(grid, cand, opts, var="Llr", better="higher") {
  grid$X <- grid[[var]]
  cand$X <- cand[[var]]
  if (better=="higher") {
    midp <- min(cand$X)
    delta <- (max(cand$X) - midp)
    lowp <- midp - delta
    grid[X<lowp, X:=lowp]
    labelcol <- "blue"
  } else if (better=="lower") {
    midp <- max(cand$X)
    delta <- (midp-min(cand$X))
    lowp <- midp + delta
    grid[X>lowp, X:=lowp]
    labelcol <- "red"
  }
  plt <- ggplot(grid) +
    aes(x=p0, y=P0, fill=X) +
    geom_tile() +
    scale_x_continuous(breaks=seq(0, 1, 0.1)) +
    scale_y_continuous(breaks=seq(opts$opt.P.lo, opts$opt.P.hi, 1)) +
    scale_fill_gradient2(low="blue", mid="white", high="red", name=var, midpoint=midp) +
    geom_point(data=cand, size=2) +
    theme_pubr(legend="right")
  if (!is.null(cand$cand)) {
    plt <- plt + geom_text(aes(label=cand, x=p0+0.01), data=cand,
                           size=6, fontface = "bold", color=labelcol, hjust = 0)
  }
  return(plt)
}


#opt
.get.nC <- function(cnv, opts) {
  if (is.null(opts$sex)) {
    sex <- .detect.sex(cnv$var, cnv$tile)
  } else {
    sex <- opts$sex
  }
  copy <- rep(2L, length(cnv$tile))
  if (sex=="male") {
    x.copy <- 1L
    y.copy <- 1L
  } else if (sex=="female"){
    x.copy <- 2L
    y.copy <- 0L
  } else {
    stop("wrong sex.")
  }
  copy[as.logical(seqnames(cnv$tile) %in% c("chrX"))] <- x.copy
  copy[as.logical(seqnames(cnv$tile) %in% c("chrY"))] <- y.copy
  return(copy)
}

## data generation
.lr.opt.data <- function(cnv, opts) {
  tmp <- data.table(
    seg=cnv$tile$seg,
    lr=cnv$tile$lr,
    nC=.get.nC(cnv, opts)
  )
  ##
  global.sd <- estimateSd(tmp$lr)
  if (opts$opt.local.sd) {
    tmp[,sd := sd(lr, na.rm=TRUE), by=seg]
    tmp[,sd := ifelse(nlr>30, sd, global.sd)]
  } else {
    tmp[,sd := global.sd]
  }
  return(tmp[])
}

.af.opt.data <- function(cnv, opts) {
  ##
  cols <- c("t.AF", "t.DP", "SOURCE", "mask.loose", "mask.strict", "seg")
  snpt <- as.data.table(mcols(cnv$var)[,cols])
  snpt[,idx:=.I]
  snpt[,pass:=.filterVar(cnv$tile, cnv$var, opts)]
  snpt <- snpt[(pass)]
  if (any(cnv$var$SOURCE=="target")) {
    snpt <- snpt[SOURCE=="target"]
  } else {
    snpt <- snpt[order(!mask.strict, !mask.loose, t.DP, decreasing=TRUE),
                 head(.SD, opts$opt.max.snp.per.segment), seg]
  }
  return(snpt)
}

.seg.opt.data <- function(cnv, opts) {
  tmp <- as.data.table(mcols(cnv$tile))
  tmp$len <- width(cnv$tile)
  seg.len <- width(cnv$seg)
  segt <- tmp[,.(
    len=seg.len[seg],
    nlr=sum(!is.na(lr)),
    gap=weighted.mean(gap, len),
    ntot=.N
  ),by=seg]
  setkey(segt, seg)
  return(segt)
}

.opt.llik <- function(lrC, af, seg, p, P, opts) {
  ## LR llik
  lrl.pick <- .llik.lrC.outer(lrC, p, P, opts$opt.max.sC, opts$opt.p.lr.anom, TRUE)[seg]
  setkey(lrl.pick, seg, C)
  tmp <- lrl.pick
  y1 <- tmp[(len / nlr < opts$opt.max.len.per.probe),.(
    p0 = ..p, P0 = ..P,
    hP = weighted.mean(ifelse(sC * nC/2 < opts$opt.max.C + 0.5, C  * nC/2, sC * nC/2), len), ## heuristic ploidy
    hL = sum((1 - ((2 * pmin(d, 0.5))**1.5)) * llik), ## heuristic likelihood
    tL = sum(llik), ## total segment log-likelihood
    cP = weighted.mean(C  * nC/2, len), ## clonal ploidy
    sP = weighted.mean(sC * nC/2, len), ## sub-clonal ploidy
    aD = weighted.mean(d, len), ## average sub-clonal deviation
    pH = weighted.mean(h, len), ## proportion homozygous
    pN = weighted.mean(n, len), ## proportion negative copy-number
    nlr = sum(nlr) # total number of TILEs
  )]
  ## AF llik
  if (!is.null(af)) {
    afl.pick <- .llik.afC.outer(af, p, P, opts$opt.max.C, opts$opt.p.af.anom, opts$opt.dp.af.max, TRUE)
    setkey(afl.pick, seg, C)
    tmp <- afl.pick[lrl.pick[,.(seg, C)], nomatch=0]
    y2 <- tmp[,.(
      aL = sum(llik), ## total af likelihood
      pA = sum(anom*naf)/sum(naf), ## proportion anomalous
      mse = weighted.mean(mse, naf), ## mean-squered error
      naf = sum(naf) ## total number of SNPs
    )]
    y <- c(y1, y2)
    r <- min((y$naf / y$nlr), 1)
    y$L <- r * y$hL + 1 * y$aL
  } else {
    y <- y1
    y$L <- y$hL
  }
  return(y)
}

.opt.llik.p <- function(lrC, af, seg, pi, Pi, opts) {
  ## optimize p
  p.lox <- max(pi-opts$opt.fine.p.off, opts$opt.p.lo)
  p.hix <- min(pi+opts$opt.fine.p.off, opts$opt.p.hi)
  pP <- cbind(p=seq(p.lox, p.hix, opts$opt.fine.p.res), P=Pi)
  tmp <- rbindlist(lapply(seq_len(nrow(pP)), function(i) {
    .opt.llik(lrC, af, seg, pP[[i,1]], pP[[i,2]], opts)
  }))
  opt.p <- tmp[order(-L), .SD[1]]
  return(opt.p)
}

.opt.seg <- function(cnv, data, p, P, opts) {
  lrC <- .lr.grid.lrC(data$lr, opts$opt.max.C)
  lrl.pick <- .llik.lrC.outer(lrC, p, P, opts$opt.max.sC, opts$opt.p.lr.anom, TRUE)[data$seg]
  setkey(lrl.pick, seg, C)
  if (!is.null(data$af)) {
    afl.pick <- .llik.afC.outer(data$af, p, P, opts$opt.max.C, opts$opt.p.af.anom, opts$opt.dp.af.max, TRUE)
    setkey(afl.pick, seg, C)
    tmp <- afl.pick[lrl.pick]
    tmp[nC==1, ":="(C=1, K=0)]
    tmp[nC==0, ":="(C=0, K=0)]
    tmp <- tmp[,.(seg, C, K, tL=i.llik, aL=llik, d, anom, mse, nlr, naf, len, sC=sC * nC/2)]
  } else {
    tmp <- lrl.pick[,.(seg, C, K=NA_integer_, tL=llik, aL=NA_real_, d, anom=NA_real_, mse=NA_real_, nlr, naf=NA_real_, len, sC=sC * nC/2)]
  }
  seg <- cnv$seg[tmp$seg]
  mcols(seg) <- tmp
  return(seg)
}

optData <- function(cnv, opts) {
  if (length(cnv$var)>0) {
    af <- .af.opt.data(cnv, opts)
  } else {
    af <- NULL
  }
  lr <- .lr.opt.data(cnv, opts)
  seg <- .seg.opt.data(cnv, opts)
  data <- list(lr=lr, af=af, seg=seg)
  return(data)
}

## hP (heuristic ploidy):
## Ploidy (trimmed) if subclonal copy (sC) is above max.C,
## C is not a good estimate of copy-number and we use sC
## hL (heuristic log-likelihood):
## down-weights sub-clonal fragments, this rewards models for which the majority
## of segments are clonal (low d), because d is correlated with P, this also
## rewards low ploidy models
optGrid <- function(data, opts) {
  ## grid over all p and P combinations
  pP <- .lr.grid.pP(opts$opt.grid.p.res, opts$opt.p.lo, opts$opt.p.hi, opts$opt.P.lo, opts$opt.P.hi)
  ## grid over all C's for each segment
  lrC <- .lr.grid.lrC(data$lr, opts$opt.max.C)
  ## calculate likelihood
  grid <- rbindlist(mclapply(seq_len(nrow(pP)), function(i) {
    x <- .llik.lrC.outer(lrC, pP[[i,1]], pP[[i,2]], opts$opt.max.sC, opts$opt.p.lr.anom, TRUE)[data$seg]
    ## pick segments with enough probe density
    y <- x[(len / nlr < opts$opt.max.len.per.probe),.(
      p0 = pP[[i,1]], P0 = pP[[i,2]],
      hP = weighted.mean(ifelse(sC * nC/2 < opts$opt.max.C + 0.5, C  * nC/2, sC * nC/2), len), ## heuristic ploidy
      hL = sum((1 - 2 * pmin(d, 0.5)) * llik), ## heuristic likelihood
      tL = sum(llik), ## total segment log-likelihood
      cP = weighted.mean(C  * nC/2, len), ## clonal ploidy
      sP = weighted.mean(sC * nC/2, len), ## sub-clonal ploidy
      aD = weighted.mean(d, len), ## average sub-clonal deviation
      pH = weighted.mean(h, len), ## proportion homozygous
      pN = weighted.mean(n, len), ## proportion negative copy-number
      nlr = sum(nlr) # total number of TILEs
    )]
  }, mc.cores=opts$cores))
  return(grid)
}

optFine <- function(data, grid, opts) {
  cand <- .grid.localmax(grid, res=opts$opt.cand.res, x="p0", y="P0", var="hL")
  lrC <- .lr.grid.lrC(data$lr, opts$opt.max.C)
  fine <- rbindlist(mclapply(seq_len(nrow(cand)), function(i) {
    ## setup variables
    iter <- 0
    pi <- cand[i,p0]
    Pi <- cand[i,P0]
    cPi <- cand[i,cP]
    trace <- .opt.llik(lrC, data$af, data$seg, pi, Pi, opts)
    while(TRUE) {
      iter <- iter + 1
      opt.p <- .opt.llik.p(lrC, data$af, data$seg, pi, Pi, opts)
      pj <- opt.p[,p0] # correct
      Pj <- opt.p[,hP]
      cPj <- opt.p[,cP]
      trace <- rbind(trace, opt.p)
      if (
        (iter==opts$opt.cand.max.iter) ||
        (cPi==cPj)
      ) {
        break
      }
      pi <- pj
      Pi <- Pj
      cPi <- cPj
    }
    trace[,":="(cand=i, iter=.I-1)]
    return(trace)
  }, mc.cores=opts$cores))
  return(fine)
}

optSegs <- function(cnv, data, fine, opts) {
  tmp <- fine[order(-iter),.SD[1], cand]
  segs <- mclapply(seq_len(nrow(tmp)), function(i) {
    .opt.seg(cnv, data, tmp[i,p0], tmp[i, P0], opts)
  }, mc.cores=opts$cores)
  return(segs)
}

#misc
.seqinfo <- function(genome) {
  genome <- tolower(genome)
  if (genome %in% c("hg38", "grch38")) {
    seqi <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
    seqi <- keepStandardChromosomes(seqi)
    seqi <- dropSeqlevels(seqi, "chrM")
  } else if (genome %in% c("hg19", "grch37")) {
    seqi <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
    seqi <- keepStandardChromosomes(seqi)
    seqi <- dropSeqlevels(seqi, "chrM")
  } else if (genome %in% c("mm10", "grcm38")) {
    seqi <- seqinfo(BSgenome.Hsapiens.UCSC.mm10)
    seqi <- keepStandardChromosomes(seqi)
    seqi <- dropSeqlevels(seqi, "chrM")
  } else {
    stop(sprintf("Genome: %s not supported", genome))
  }
  return(seqi)
}

.tilegenome <- function(genome, tile.width) {
  seqi <- .seqinfo(genome)
  tile <- tileGenome(seqi, tilewidth=tile.width, cut.last.tile.in.chrom=TRUE)
  return(tile)
}

.robust.import <- function(fn, genome) {
  tmp <- import(fn)
  seqi <- .seqinfo(genome)
  shared.levels <- intersect(seqlevels(tmp), seqlevels(seqi))
  tmp <- keepSeqlevels(tmp, shared.levels, pruning.mode="coarse")
  seqlevels(tmp) <- seqlevels(seqi)
  seqinfo(tmp) <- seqi
  return(tmp)
}

.robust.import.vcf <- function(fn, param, genome) {
  var <- readVcf(fn, param=param)
  seqi <- .seqinfo(genome)
  shared.levels <- intersect(seqlevels(var), seqlevels(seqi))
  var <- keepSeqlevels(var, shared.levels, pruning.mode="coarse")
  seqlevels(var) <- seqlevels(seqi)
  return(var)
}

.log.sum.exp <- function(x) {
  offset <- max(x)
  log(sum(exp(x - offset))) + offset
}

.dtRound <- function(dt) {
  for (i in names(dt)) {
    if (class(dt[,get(i)])=="numeric") {
      dt[,":="((i), round(get(i), 6))]
    }
  }
  return(dt)
}

max.na.rm <- function(x) max(x, na.rm=TRUE)

merge.list <- function (x, y) {
  if (length(x) == 0)
    return(y)
  if (length(y) == 0)
    return(x)
  i = match(names(y), names(x))
  i = is.na(i)
  if (any(i))
    x[names(y)[which(i)]] = y[which(i)]
  return(x)
}

.absMedDiff <- function(x, y) {
  abd <- abs(median(x, na.rm=TRUE)-median(y, na.rm=TRUE))
  return(abd)
}

plog_sum_exp <- function(u, v) {
  m <- pmax(u, v)
  m + log(exp(u - m) + exp(v - m))
}

.grid.localmax <- function(grid, res, x, y, var) {
  mx <- as.matrix(dcast.data.table(grid, reformulate(x, response=y), value.var=var)[,-1])
  r <- raster(mx)
  ## nearest odd integer >= to 3
  rres <- max(2*floor((nrow(r)*res)/2)+1, 3)
  cres <- max(2*floor((ncol(r)*res)/2)+1, 3)
  wind <- matrix(1, nrow=rres, ncol=cres)
  localmax <- focal(r, fun = max.na.rm, w = wind, pad=TRUE, padValue=NA)
  cand <- grid[Which(localmax==r, cells=TRUE)]
  return(cand)
}

#lr
.smoothOutliers <- function(y, tile, opts) {
  ## remove gross outliers
  y.tar.n <- sum(!is.na(y[ tile$target]))
  if (y.tar.n>1) {
    y[ tile$target] <- .runSmooth(y[ tile$target], tile[ tile$target])
  }
  y.off.n <- sum(!is.na(y[!tile$target]))
  if (y.off.n>1) {
    y[!tile$target] <- .runSmooth(y[!tile$target], tile[!tile$target])
  }
  return(y)
}

.rawLogRatio <- function(t.cov, n.cov, opts) {
  lr.raw <- log2(t.cov / n.cov)
  lr.raw <- ifelse(is.finite(lr.raw), lr.raw, NA_real_)
  return(lr.raw)
}

.scaleLogRatio <- function(lr, gt, opts) {
  weight <- width(gt)
  lr.mean.tgt <- weighted.mean(lr[ gt$target & gt$hq], weight[ gt$target & gt$hq], na.rm=TRUE)
  lr.mean.off <- weighted.mean(lr[!gt$target & gt$hq], weight[!gt$target & gt$hq], na.rm=TRUE)
  lr[ gt$target] <- lr[ gt$target] - lr.mean.tgt
  lr[!gt$target] <- lr[!gt$target] - lr.mean.off
  return(lr)
}

.gcLogRatio <- function(lr, gt, opts) {
  ## normalize, smooth, and gc-correct
  if (opts$gc.adjust.trend) {
    for (sel in unique(gt$target)) {
      lr.sel <- lr[gt$target==sel]
      gt.sel <- gt[gt$target==sel]
      weight.sel <- ifelse(gt.sel$hq, 1, 0)
      span.sel <- ifelse(sel, opts$gc.adjust.span.on, opts$gc.adjust.span.off)
      gc.residuals.sel <- limma::loessFit(y=lr.sel, x=gt.sel$gc,
                                          weight=weight.sel, min.weight=1e-9,
                                          span=span.sel)$residuals
      if (opts$gc.adjust.offset) {
        lr.offset.sel <- lm(gc.residuals.sel~lr.sel, weights=weight.sel)$coefficients[1]
        gc.residuals.sel <- gc.residuals.sel - lr.offset.sel
      }
      if (sel) {
        gc.range <- gt$gc > opts$gc.adjust.on[1] & gt$gc < opts$gc.adjust.on[2]
        lr[gt$target & gc.range] <- gc.residuals.sel[gc.range[gt$target]]
      } else {
        gc.range <- gt$gc > opts$gc.adjust.off[1] & gt$gc < opts$gc.adjust.off[2]
        lr[!gt$target & gc.range] <- gc.residuals.sel[gc.range[!gt$target]]
      }
    }
  }
  lr[!is.finite(lr)] <- NA_real_
  return(lr)
}

.smoothLogRatio <- function(lr, tile, opts) {
  lr.smooth <- unname(unlist(lapply(split(lr, tile$arm), function(y) {
    if (length(y) >= opts$lr.smooth.window) {
      y.smooth <- suppressWarnings(
        hybrid.filter(y, opts$lr.smooth.window, method=c("MEAN", "MED"))$level[["MEAN"]]
      )
    } else {
      y.smooth <- y
    }
    return(y.smooth)
  })))
  return(lr.smooth)
}

#.lr.grid.pP <- function(grid.p.res, p.lo, p.hi, P.lo, P.hi) {
p <- seq(p.lo, p.hi, grid.p.res)
P <- seq(P.lo, P.hi, length.out=length(p))
pP <- as.matrix(expand.grid(p=p, P=P))
return(pP)
}

## lr data grid generation
.lr.grid.lrC <- function(lr, max.C) {
  ## only keep tiles with coverage for llik calculations
  lr <- lr[!is.na(lr)]
  lrC <- expand.grid(lr=lr$lr, C=0:max.C)
  lrC$seg <- lr$seg
  lrC$sd <- lr$sd
  lrC$nC <- lr$nC
  setDT(lrC)
  setkey(lrC, C, seg)
  return(lrC[])
}

## llik_lr
.llik.lrC.inner <- function(lrC, p, P, p.lr.anom) {
  lrC[,":="(
    ## normal distribution to model clonal variants
    norm = dnorm(lr, mean=log2((p * C + (1-p) * 2) / ((p * P) + (1-p) * 2)),
                 sd=sd, log=TRUE),
    ## uniform distribution to model anomalies
    unif = dunif(0, min=-6, max=6, log=TRUE)
  )]
  lrC[,":="(
    llik=plog_sum_exp(norm + log(1-p.lr.anom), unif + log(p.lr.anom))
  )]
  x <- lrC[,.( ## sum by C and seg
    p=p, P=P,
    lr = mean(lr),
    llik = sum(llik),
    nC = nC[1]
  ), by=.(C, seg)]
  return(x)
}

.llik.lrC.outer <- function(lrC, p, P, max.sC, p.lr.anom, collapse) {
  ## compute likelihood for each segment and each C
  x <- .llik.lrC.inner(lrC, p, P, p.lr.anom)
  ## ML sub-clonal copy-number
  raw.sC <- x[,(2^(lr) * ((p * P) + (1-p) * 2) - ((1-p) * 2)) / p]
  x[,sC:=pmax(pmin((raw.sC), max.sC), 0)]
  ## sub-clonal deviation
  x[,d:=abs(C - sC) * nC/2]
  ## homozygous
  x[,h:=(C==0)]
  ## negative copy-number
  x[,n:=(raw.sC < 0)]
  ## for each segment pick C with highest llik
  if (collapse) {
    x <- x[order(-llik),.SD[1],by=seg]
  }
  setkey(x, seg)
  return(x[])
}

#llilk_af
.af.grid.MC <- function(p, P, max.C) {
## M - number of chromosome with variant
## C - total number of chromosomes
## Ef - expected frequency of variant
MC <- as.data.table(expand.grid(M=0:max.C,C=0:max.C))[M<=C]
MC[,":="(
  K=pmin(M,C-M),
  Ef=(p * M + 1 * (1-p)) / (p * C + 2 * (1-p)),
  prior.M=ifelse(M==C-M, log(1), log(1/2))
)]
return(MC[])
}

.llik.afC.inner <- function(af, MC, p.af.anom, dp.af.max) {
  ## af data grid generation
  tmp1 <- rbindlist(rep(list(af), nrow(MC)))[order(idx)]
  tmp2 <- rbindlist(rep(list(MC), nrow(af)))
  afC <- cbind(tmp1, tmp2)
  afC[,":="(
    beta=dbeta(
      Ef,
      shape1=   t.AF  * pmin(t.DP, dp.af.max) + 1,
      shape2=(1-t.AF) * pmin(t.DP, dp.af.max) + 1,
      log=TRUE
    ),
    unif=1 # dbeta(x, 1, 1) or dunif(x, 0, 1)
  )]
  bllik <- afC[,beta + prior.M + log(1-p.af.anom)]
  ullik <- afC[,unif + prior.M + log(p.af.anom)]
  afC[,":="(
    llik=plog_sum_exp(bllik, ullik),
    anom=bllik < ullik,
    se=(Ef-t.AF)**2
  )]
  ## pick best M per K C
  x <- afC[order(-llik),.SD[1],.(idx, K, C)]
  ## total/average by K C
  x <- x[,.(
    llik=sum(llik),
    anom=mean(anom),
    mse=mean(se),
    naf=.N
  ), .(seg, K, C)]
  return(x[])
}

.llik.afC.outer <- function(af, p, P, max.C, p.af.anom, dp.af.max, collapse) {
  ## get af
  MC <- .af.grid.MC(p, P, max.C)
  afs <- split(af, af$seg)
  x <- rbindlist(lapply(afs, function(af) {
    .llik.afC.inner(af, MC, p.af.anom, dp.af.max)
  }))
  if (collapse) {
    x <- x[order(-llik),.SD[1],.(seg, C)]
  }
  return(x[])
}

#gc
.correctGC <- function(lr, gt, opts) {
  ## normalize, smooth, and gc-correct
  if (opts$gc.adjust.trend) {
    for (sel in unique(gt$target)) {
      lr.sel <- lr[gt$target==sel]
      gt.sel <- gt[gt$target==sel]
      weight.sel <- ifelse(gt.sel$blacklist==0 & gt.sel$unmasked>0.9, 1, 0)
      span.sel <- ifelse(sel, opts$gc.adjust.span.on, opts$gc.adjust.span.off)
      gc.residuals.sel <- limma::loessFit(y=lr.sel, x=gt.sel$gc,
                                          weight=weight.sel, min.weight=1e-9,
                                          span=span.sel)$residuals
      if (opts$gc.adjust.offset) {
        lr.offset.sel <- lm(gc.residuals.sel~lr.sel, weights=weight.sel)$coefficients[1]
        gc.residuals.sel <- gc.residuals.sel - lr.offset.sel
      }
      if (sel) {
        gc.range <- gt$gc > opts$gc.adjust.on[1] & gt$gc < opts$gc.adjust.on[2]
        lr[gt$target & gc.range] <- gc.residuals.sel[gc.range[gt$target]]
      } else {
        gc.range <- gt$gc > opts$gc.adjust.off[1] & gt$gc < opts$gc.adjust.off[2]
        lr[!gt$target & gc.range] <- gc.residuals.sel[gc.range[!gt$target]]
      }
    }
  }
  lr[!is.finite(lr)] <- NA_real_
  return(lr)
}

#filter
.filterGenomeGermlineHets <- function(var, opts) {
  ## genome germline hets
  germline <-
    ## germline SNV
    !var$SOMATIC &
    var$TYPE=="SNV" &
    ## heterozygous in normal
    var$n.GT %in% c("0/1", "1/0") &
    var$n.DP > opts$baf.min.het.dp &
    (
      (!var$mask.strict & !var$mask.loose &
         var$n.AF > opts$baf.het.range[1]+0.00 & var$n.AF < opts$baf.het.range[2]-0.00) |
        (var$mask.strict & !var$mask.loose &
           var$n.AF > opts$baf.het.range[1]+0.05 & var$n.AF < opts$baf.het.range[2]-0.05) |
        (var$mask.strict &  var$mask.loose &
           var$n.AF > opts$baf.het.range[1]+0.10 & var$n.AF < opts$baf.het.range[2]-0.10)
    )
  return(germline)
}

.filterTargetGermlineHets <- function(var, opts) {
  ## target germline hets
  germline <-
    ## germline SNV
    !var$SOMATIC &
    var$TYPE=="SNV" &
    ## heterozygous in normal
    var$n.GT %in% c("0/1", "1/0") &
    var$n.DP > opts$baf.min.het.dp &
    (
      (!var$mask.strict & !var$mask.loose &
         var$n.AF > opts$baf.het.range[1]+0.00 & var$n.AF < opts$baf.het.range[2]-0.00) |
        (var$mask.strict & !var$mask.loose &
           var$n.AF > opts$baf.het.range[1]+0.02 & var$n.AF < opts$baf.het.range[2]-0.02) |
        (var$mask.strict &  var$mask.loose &
           var$n.AF > opts$baf.het.range[1]+0.04 & var$n.AF < opts$baf.het.range[2]-0.04)
    )
  return(germline)
}

.filterGenomeHets <- function(var, opts) {
  ## genome germline hets
  het <-
    ## germline SNV
    !var$SOMATIC &
    var$TYPE=="SNV" &
    ## heterozygous in normal
    var$t.GT %in% c("0/1", "1/0") &
    var$t.DP > opts$baf.min.het.dp &
    var$t.AF > 0.05 & var$n.AF < 0.95 &
    !var$mask.strict
  return(het)
}

.filterTargetHets <- function(var, opts) {
  ## target germline hets
  het <-
    ## germline SNV
    !var$SOMATIC &
    var$TYPE=="SNV" &
    ## heterozygous
    var$t.GT %in% c("0/1", "1/0") &
    var$t.DP > opts$baf.min.het.dp &
    var$t.AF > 0.05 & var$t.AF < 0.95 &
    !var$mask.strict
  return(het)
}

.filterGenomeCoverage <- function(var, opts) {
  t.hi.dp <- quantile(var$t.DP, 0.95, na.rm=TRUE)
  t.lo.dp <- quantile(var$t.DP, 0.05, na.rm=TRUE)
  n.hi.dp <- quantile(var$n.DP, 0.95, na.rm=TRUE)
  n.lo.dp <- quantile(var$n.DP, 0.05, na.rm=TRUE)
  ## high-quality
  covered <-
    var$t.DP > t.lo.dp & var$t.DP < t.hi.dp &
    var$n.DP > n.lo.dp & var$n.DP < n.hi.dp &
    var$t.DP > opts$baf.min.genome.dp
  return(covered)
}

.filterTargetCoverage <- function(var, opts) {
  t.hi.dp <- quantile(var$t.DP, 0.99, na.rm=TRUE)
  t.lo.dp <- quantile(var$t.DP, 0.01, na.rm=TRUE)
  n.hi.dp <- quantile(var$n.DP, 0.99, na.rm=TRUE)
  n.lo.dp <- quantile(var$n.DP, 0.01, na.rm=TRUE)
  ## coverage
  covered <-
    var$t.DP > t.lo.dp & var$t.DP < t.hi.dp &
    var$n.DP > n.lo.dp & var$n.DP < n.hi.dp &
    var$t.DP > opts$baf.min.target.dp
  return(covered)
}

.filterOnTarget <- function(tile, var, opts) {
  splash <- tile[tile$target] + opts$tile.shoulder
  ontarget <- var %over% splash
  return(ontarget)
}

.filterTumorNormal <- function(tile, var, opts) {
  pass <- logical(length(var))
  if (any(var$SOURCE=="genome")) {
    var.genome <- var[var$SOURCE=="genome"]
    germline.genome <- .filterGenomeGermlineHets(var.genome, opts)
    covered.genome <- .filterGenomeCoverage(var.genome, opts)
    pass[var$SOURCE=="genome"] <- germline.genome & covered.genome
  }
  if (any(var$SOURCE=="target")) {
    var.target <- var[var$SOURCE=="target"]
    germline.target <- .filterTargetGermlineHets(var.target, opts)
    covered.target <- .filterTargetCoverage(var.target, opts)
    ontarget <- .filterOnTarget(tile, var.target, opts)
    pass[var$SOURCE=="target"] <- germline.target & covered.target & ontarget
  }
  pass[is.na(pass)] <- FALSE
  return(pass)
}

.filterTumorOnly <- function(tile, var, opts) {
  pass <- logical(length(var))
  if (any(var$SOURCE=="genome")) {
    var.genome <- var[var$SOURCE=="genome"]
    het.genome <- .filterGenomeHets(var.genome, opts)
    covered.genome <- .filterGenomeCoverage(var.genome, opts)
    pass[var$SOURCE=="genome"] <- het.genome & covered.genome
  }
  if (any(var$SOURCE=="target")) {
    var.target <- var[var$SOURCE=="target"]
    het.target <- .filterTargetHets(var.target, opts)
    covered.target <- .filterTargetCoverage(var.target, opts)
    ontarget <- .filterOnTarget(tile, var.target, opts)
    pass[var$SOURCE=="target"] <- het.target & covered.target & ontarget
  }
  pass[is.na(pass)] <- FALSE
  return(pass)
}

.filterVar <- function(tile, var, opts) {
  tumor.normal <- "n.GT" %in% names(mcols(var)) ## FIXME
  if (tumor.normal) {
    pass <- .filterTumorNormal(tile, var, opts)
  } else {
    pass <- .filterTumorOnly(tile, var, opts)
  }
  return(pass)
}

#external
.runMosdepth <- function(bam, by, fasta, cores) {
  if (is.character(by) && file.exists(by)) {
    out.col <- "V5"
  } else if (!is.na(suppressWarnings(as.integer(by))) ) {
    out.col <- "V4"
  } else {
    stop("by should be a file or window-size")
  }
  prefix <- tempfile("mosdepth_")
  ret <- system2("mosdepth", sprintf("-x -f %s -F 1796 -n -t%s -b %s %s %s", fasta, cores, by, prefix, bam))
  out.fn <- list.files(dirname(prefix), paste0(basename(prefix), ".regions.bed.gz$"),
                       full.names = TRUE)
  if (!ret & file.exists(out.fn)) {
    ## this makes the assumption that the input bed
    ## output and regions.bed.gz file have the same
    ## order... for sure the input bed has to be
    ## sorted.
    out <- fread(out.fn)[[out.col]]
  } else {
    stop(sprintf("mosdepth run failed: %s", ret))
  }
  out.fns <- list.files(dirname(prefix), paste0(basename(prefix)), full.names = TRUE)
  rets <- sapply(out.fns, unlink)
  if (any(rets)) {
    stop("could not remove temp. files")
  }
  return(out)
}

.runMosdepthMulti <- function(bams, ...) {
  out <- rowSums(mapply(FUN=function(bam) {
    out1 <- .runMosdepth(bam, ...)
  }, bams))
  return(out)
}

.runMosdepthTile <- function(bams, gt, fasta, cores) {
  bed <- tempfile("mosdepth_", fileext=".bed")
  export(granges(gt), bed)
  out <- .runMosdepthMulti(bams, bed, fasta, cores)
  unlink(bed)
  return(out)
}

.runCBS <- function(y, y.weights=NULL, args) {
  n <- length(y)
  chrom <- rep(1, n)
  maploc <- 1:n
  genomdat <- y
  cna <- DNAcopy::CNA(genomdat, chrom, maploc)
  ## this is crazy
  if (!is.null(y.weights)) {
    inp <- c(list(cna, weights=y.weights), args)
  } else {
    inp <- c(list(cna), args)
  }
  capture.output(
    res <- do.call(DNAcopy::segment, inp)
  )
  bkp <- res$output$loc.end[-length(res$output$loc.end)]
  if (length(bkp)==0) {
    bkp <- integer()
  }
  return(bkp)
}

.runSmooth <- function(lr, gr) {
  obj <- CNA(lr, as.integer(seqnames(gr)), floor((start(gr)+end(gr))/2),
             data.type = "logratio",
             sampleid = "sample")
  adj <- smooth.CNA(obj)$sample
  return(adj)
}

#cov
.normCoverage <- function(cov, gt, opts) {
  ## normalized coverage
  cov.n <- cov*width(gt)
  ## on-target
  cov.n.tar.total <- sum(cov.n[gt$target & gt$hq])
  cov.n.tar <- cov.n[gt$target]
  cov.n.tar <- cov.n.tar/(cov.n.tar.total/1e6)
  cov.n.tar <- 1000*cov.n.tar/width(gt[ gt$target])
  cov.n[ gt$target] <- cov.n.tar
  ## off-target
  cov.n.off.total <- sum(cov.n[!gt$target & gt$hq])
  cov.n.off <- cov.n[!gt$target]
  cov.n.off <- cov.n.off/(cov.n.off.total/1e6)
  cov.n.off <- 1000*cov.n.off/width(gt[!gt$target])
  cov.n[!gt$target] <- cov.n.off
  return(cov.n)
}

#caller
.importSentieon <- function(vcf, opts) {
  param <-ScanVcfParam(
    info=c("AF"),
    geno=c("GT", "AD", "DP", "GQ")
  )
  if (!is.null(opts$which)) {
    vcfWhich(param) <- which
  }
  var <- .importVcf(vcf, param, opts)
  var.gr <- rowRanges(var)
  var.gr$mask.loose <- info(var)$mask.loose
  var.gr$mask.strict <- info(var)$mask.strict
  var.gr$SOMATIC <- FALSE
  var.gr$TYPE <- "SNV"
  t.AD <- sapply(geno(var)$AD[,1], "[", 2)
  t.DP <- geno(var)$DP[,1]
  var.gr$t.GT <- geno(var)$GT[,1]
  var.gr$t.AF <- t.AD / t.DP
  var.gr$t.DP <- t.DP
  if (ncol(var)==2) {
    n.AD <- sapply(geno(var)$AD[,2], "[", 2)
    n.DP <- geno(var)$DP[,2]
    var.gr$n.GT <- geno(var)$GT[,2]
    var.gr$n.AF <- n.AD / n.DP
    var.gr$n.DP <- n.DP
  }
  return(var.gr)
}

getVcfImporter <- function(caller) {
  ## import
  if (is.null(caller)) {
    func <- NULL
  } else if (caller=="sentieon") {
    func <- .importSentieon
  } else {
    stop("Variant caller not supported.")
  }
  return(func)
}

#baf

.getTargetHits <- function(gt, snp, opts) {
  hits <- findOverlaps(snp, gt, maxgap = opts$tile.shoulder-1)
  hits <- hits[gt[subjectHits(hits)]$target] # prefer assignment to target
  hits <- hits[!duplicated(queryHits(hits))] # if snp close to two targets pick one
  return(hits)
}

.getGenomeHits <- function(gt, snp, opts) {
  hits <- findOverlaps(snp, gt, maxgap = opts$tile.shoulder-1)
  return(hits)
}

.getHits <- function(gt, snp, opts) {
  if (length(opts$capture)==0) {
    hits <- .getGenomeHits(gt, snp, opts)
  } else {
    hits <- .getTargetHits(gt, snp, opts)
  }
  return(hits)
}

.addBafTile <- function(gt, var, opts) {
  ## filter variants to SNPs
  pass <- .filterVar(gt, var, opts)
  snp <- var[pass]
  hits <- .getHits(gt, snp, opts)
  if (length(hits)>0) {
    bad=ifelse(snp[queryHits(hits)]$t.AF < 0.5,
               round(   snp[queryHits(hits)]$t.AF  * snp[queryHits(hits)]$t.DP),
               round((1-snp[queryHits(hits)]$t.AF) * snp[queryHits(hits)]$t.DP))
    tmp <- data.table(
      idx=subjectHits(hits),
      bad=bad,
      depth=snp[queryHits(hits)]$t.DP
    )
    setkey(tmp, idx)
    tmp <- tmp[J(1:length(gt))]
    tmp <- tmp[,.(
      baf=sum(bad)/sum(depth),
      depth=sum(depth)
    ), by=idx]
    tmp <- tmp[,":="(
      ## weight proportional to inverse varince ~ 1/sqrt(n)
      baf.weight = 1/sqrt(0.25/pmin(pmax(1, depth), opts$baf.max.eff.dp))
    )]
    gt$baf <- tmp$baf
    gt$baf.weight <- tmp$baf.weight
    gt$baf <- .smoothOutliers(gt$baf, gt, opts)
  } else {
    gt$baf <- NA_real_
  }
  return(gt)
}

#api
importFeat <- function(gtf, opts) {
  ann <- .robust.import(gtf, opts$genome)
  feat <- ann[ann$type==opts$feature]
  return(feat)
}

#' @export
importVcf <- function(vcf, opts) {
  func <- getVcfImporter(opts$caller)
  var <- unlist(GRangesList(unname(lapply(names(vcf), function(name) {
    fn <- vcf[[name]]
    v <- func(fn, opts)
    mcols(v)$SOURCE <- name
    return(v)
  }))))
  return(var)
}

#' @export
importTile <- function(opts) {
  if (is.na(opts$tile.width)) {
    target.fun <- getTargetTiles
  } else {
    target.fun <- getGenomeTiles
  }
  tile <- target.fun(opts$capture, opts)
  return(tile)
}

#' @export
createCnv <- function(vcf, opts) {
  var <- importVcf(vcf, opts)
  tile <- importTile(opts)
  cnv <- list(tile=tile, var=var)
  return(cnv)
}

#' @export
createPool <- function(cnv.fns, opts) {
  pool.data <- .importPoolData(cnv.fns, opts)
  pool <- .createPool(pool.data, opts)
  return(pool)
}

#' @export
addCoverage <- function(cnv, t.bam, n.bam, opts) {
  ## normalize by sequencing depth
  t.cov.raw <- .runMosdepthTile(t.bam, cnv$tile, opts$fasta, opts$cores)
  t.cov <- .normCoverage(t.cov.raw, cnv$tile, opts)
  if (!is.null(n.bam)) {
    n.cov.raw <- .runMosdepthTile(n.bam, cnv$tile, opts$fasta, opts$cores)
    n.cov <- .normCoverage(n.cov.raw, cnv$tile, opts)
  } else {
    n.cov.raw <- NA_real_
    n.cov <- NA_real_
  }
  mcols(cnv$tile) <- cbind(mcols(cnv$tile), cbind(t.cov.raw, n.cov.raw, t.cov, n.cov))
  return(cnv)
}

#' @export
addKnnLogRatio <- function(cnv, pool, opts) {
  p.cov <- .poolKnnCoverage(cnv, pool, opts)
  cnv$tile$lr <- .rawLogRatio(cnv$tile$t.cov, p.cov, opts)
  return(cnv)
}

#' @export
addPcaLogRatio <- function(cnv, pool, opts) {
  cnv$tile$lr <- .poolPcaDenoise(cnv, pool, opts)
  return(cnv)
}

#' @export
addIcaLogRatio <- function(cnv, pool, opts) {
  cnv$tile$lr <- .poolIcaDenoise(cnv, pool, opts)
  return(cnv)
}

#' @export
addMeanLogRatio <- function(cnv, pool, opts) {
  cnv$tile$lr <- .rawLogRatio(cnv$tile$t.cov, 1, opts)
  return(cnv)
}

#' @export
addPairLogRatio <- function(cnv, opts) {
  cnv$tile$lr <- .rawLogRatio(cnv$tile$t.cov, cnv$tile$n.cov, opts)
  return(cnv)
}

#' @export
addScaleLogRatio <- function(cnv, opts) {
  cnv$tile$lr <- .scaleLogRatio(cnv$tile$lr, cnv$tile, opts)
  return(cnv)
}

#' @export
addGcLogRatio <- function(cnv, opts) {
  cnv$tile$lr <- .gcLogRatio(cnv$tile$lr, cnv$tile, opts)
  return(cnv)
}

#' @export
addSmoothLogRatio <- function(cnv, opts) {
  lr.smooth <- cnv$tile$lr
  if (opts$lr.smooth=="hybrid") {
    ## hybrid smooth only on targeted
    lr.smooth[cnv$tile$target] <-
      .smoothLogRatio(lr.smooth[cnv$tile$target], cnv$tile[cnv$tile$target], opts)
  } else if (opts$lr.smooth=="outlier") {
    lr.smooth <- .smoothOutliers(lr.smooth, cnv$tile, opts)
  }
  cnv$tile$lr <- lr.smooth
  return(cnv)
}

#' @export
addBaf <- function(cnv, opts) {
  ## add BAF to segmentation tile
  cnv$tile <- .addBafTile(cnv$tile, cnv$var, opts)
  return(cnv)
}

#' @export
addJointSegment <- function(cnv, opts) {
  ## estimated standard-deviation
  hq.lr <- cnv$tile[cnv$tile$unmasked > 0.25 & cnv$tile$blacklist==0.0]$lr
  sd.lr <- estimateSd(hq.lr) * opts$seg.sd.lr.penalty
  sd.baf <- estimateSd(cnv$tile$baf) * opts$seg.sd.baf.penalty
  sd.baf <- ifelse(!is.finite(sd.baf), .Machine$double.eps, sd.baf)
  ## create segmentation
  print("tile")
  cnv$tile <- .addJointSeg(cnv$tile, sd.lr, sd.baf, opts$seg.method,
                           opts$tile.width, opts$seg.len.min, opts$seg.cbs.lr,
                           opts$seg.cbs.baf, opts$seg.rbs.selection, opts$seg.sd.prune, opts$seg.len.prune)
  print("seg")
  cnv$seg <- unname(unlist(range(split(cnv$tile, cnv$tile$seg))))
  ## seg -> var
  print("last step")
  cnv$var$seg <- findOverlaps(cnv$var, cnv$seg, select="first", maxgap = opts$tile.shoulder-1)
  return(cnv)
}

#' @export
addLogRatio <- function(cnv, pool, opts) {
  if (is.null(pool)) {
    if (all(is.na(cnv$tile$n.cov))) {
      ## no matched normal
      cnv <- addMeanLogRatio(cnv, opts)
    } else {
      ## has matched normal
      cnv <- addPairLogRatio(cnv, opts)
    }
  } else {
    if (pool$method=="fake") {
      cnv <- addFakeLogRatio(cnv, pool, opts)
    }
    else if (pool$method=="pca") {
      cnv <- addPcaLogRatio(cnv, pool, opts)
    }
    else if (pool$method=="ica") {
      cnv <- addIcaLogRatio(cnv, pool, opts)
    }
  }
  cnv <- addScaleLogRatio(cnv, opts)
  cnv <- addGcLogRatio(cnv, opts)
  cnv <- addSmoothLogRatio(cnv, opts)
  return(cnv)
}

#' @export
addSegment <- function(cnv, opts) {
  cnv <- addBaf(cnv, opts)
  cnv <- addJointSegment(cnv, opts)
  return(cnv)
}

#' @export
getOpts <- function(settings, capture, opts=list()) {
  ENV <- new.env(parent = .BaseNamespaceEnv)
  source(settings, local=ENV)
  opts$capture <- capture
  opts <- merge.list(opts, ENV$OPTS)
  return(opts)
}

#' @export
segFeat <- function(feat, cnv, opts) {
  ## seg -> feat
  hit.prom <- as.data.table(findOverlaps(promoters(feat, 0, 1), cnv$seg))
  hit.feat <- as.data.table(findOverlaps(feat, cnv$seg))
  hit.feat <- hit.feat[!(queryHits %in% hit.prom$queryHits)]
  hit.feat <- hit.feat[,.SD[1],by=queryHits]
  hit <- rbind(hit.prom, hit.feat)
  feat$seg <- NA_integer_
  feat$seg[hit$queryHits] <- hit$subjectHits
  ## tile -> feat
  tile.data <- as.data.table(mcols(cnv$tile))
  tmp.seg <- tile.data[,.(n.seg=.N, lr.seg=mean(lr, na.rm=TRUE), mzd.seg=mean(baf, na.rm=TRUE)), by=seg]
  mcols(feat) <- cbind(mcols(feat), .dtRound(tmp.seg[feat$seg, -1]))
  return(feat)
}

#' @export
segGtf <- function(cnv) {
  tmp1 <- as.data.table(mcols(cnv$tile))
  tmp1 <- tmp1[,.(lr=mean(lr, na.rm=TRUE), mzd=mean(baf, na.rm=TRUE), .N),by=seg]
  tmp2 <- as.data.table(cnv$seg)[,.(chr=seqnames, start, end, seg=.I)]
  setkey(tmp1, seg)
  setkey(tmp2, seg)
  seg.stat <- tmp1[tmp2,.(chr, start, end, N, lr, mzd)]
  seg.stat.gr <- with(seg.stat, GRanges(chr, IRanges(start, end), "*", N, lr, mzd))
  return(seg.stat.gr)
}

opts <- list(
  
  ## tiling
  tile.width=10000,
  tile.min.gap=NA_integer_,
  tile.shoulder=0,
  
  ## log-ratio
  lr.smooth="outlier",
  lr.smooth.window=21,
  
  ## BAF
  baf.max.eff.dp=300,
  baf.het.range=c(0.3, 0.7),
  baf.min.het.dp=6,
  baf.min.target.dp=24,
  baf.min.genome.dp=12,
  
  ## segmentation
  seg.strategy="joint",
  seg.method="CBS",
  seg.sd.prune=TRUE,
  seg.sd.lr.penalty=1,
  seg.sd.baf.penalty=1,
  seg.len.prune=TRUE,
  seg.len.min=2,
  seg.cbs.baf=list(alpha=0.01, trim=0.025, min.width=2),
  seg.cbs.lr=list(alpha=0.01, trim=0.025, min.width=2),
  seg.rbs.selection="Lebarbier",
  
  ## GC-content
  gc.adjust.trend=TRUE,
  gc.adjust.offset=TRUE,
  gc.adjust.span.on=0.5,
  gc.adjust.span.off=NA_real_,
  gc.adjust.on=c(0.3, 0.7),
  gc.adjust.off=c(NA_real_, NA_real_),
  
  ## pool
  pool.method="pca",
  pool.lo.cov=0.25,
  pool.hi.cov=4,
  pool.hi.zero=0.25,
  pool.sd.out=3,
  pool.n.comp=20,
  
  ## llik
  opt.local.sd=FALSE,
  opt.max.C=12,
  opt.max.sC=20,
  opt.max.len.per.probe=1e6,
  opt.max.snp.per.segment=200,
  opt.p.af.anom=0.001,
  opt.p.lr.anom=0.001,
  opt.dp.af.max=300,
  
  ## grid
  opt.p.lo=0.10,
  opt.p.hi=0.99,
  opt.P.lo=1,
  opt.P.hi=9,
  opt.grid.p.res=0.0125,
  opt.fine.p.res=0.0050,
  opt.fine.p.off=0.0250,
  opt.cand.res=0.075,
  opt.cand.max.iter=5
  
)
