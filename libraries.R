library(cnvex) #this is a local library, might not work on other machines
library(GenomicRanges)
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
library(jointseg)
library(data.table)
library(tidyverse)
library(sigmoid)
library(matrixStats)
library(fastICA)
library(MASS)
numCores = detectCores()
registerDoParallel(numCores-1)
registerDoParallel(10)
set.seed(1024)

list_obj_sizes <- function(list_obj=ls(envir=.GlobalEnv)){
  sizes <- sapply(list_obj, function(n) object.size(get(n)), simplify = FALSE)  
  print(sapply(sizes[order(-as.integer(sizes))], function(s) format(s, unit = 'auto')))
}

as.data.frame(list_obj_sizes())
rm(list = ls(pattern = "variance_df*"))
base::load("~/Codes/Germline_variation_detection/data_june26.RData") 

rm(repeat_file,repeat_file_short,svd)
rm(tile_modal,ttile_modal,subdata,tile_coverage,right_form_of_data,sample)


rm(p14_pool,p14_tn,p36_tn,p36_pool,new_cnv,cnv_pca,cnv_org,cnv_ica)
