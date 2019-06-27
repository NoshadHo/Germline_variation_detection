library(GenomicRanges)
library(cnvex) #this is a local library, might not work on other machines
library(tidyverse)
library(rlist)
library(gplots)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(mclust)
library(gridExtra)
library(corrr)
library(plotly)
library(rlang)
library(limma) #bioconductor lib
library(gridExtra)
library(ggfortify)
library(mclust)
library(foreach)
library(doParallel)
library(Rcpp)
library(fpc)
library(dbscan)
library(factoextra)
set.seed(1024)
numCores = detectCores()
registerDoParallel(numCores-1)
set.seed(1024)