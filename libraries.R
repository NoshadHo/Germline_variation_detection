library(GenomicRanges)
library(cnvex) #this is a local library, might not work on other machines
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
library(ggplot2)
library(rbenchmark)
library(Gmedian)
library(ComplexHeatmap)
library(glmnet)
library(Boruta)
library(caret)
library(tidyverse)
numCores = detectCores()
registerDoParallel(numCores-1)
registerDoParallel(8)
set.seed(1024)

list_obj_sizes <- function(list_obj=ls(envir=.GlobalEnv)){
  sizes <- sapply(list_obj, function(n) object.size(get(n)), simplify = FALSE)  
  print(sapply(sizes[order(-as.integer(sizes))], function(s) format(s, unit = 'auto')))
}

as.data.frame(list_obj_sizes())
rm(list = ls(pattern = "variance*"))
base::load("~/Codes/Germline_variation_detection/data_june26.RData") 

rm(repeat_file,repeat_file_short,svd)
rm(tile_modal,ttile_modal,subdata,tile_coverage,right_form_of_data,sample)
