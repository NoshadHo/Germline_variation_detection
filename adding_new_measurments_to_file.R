#variance add to variance_raw------------------------------------------------------------
variance_df = variance_sex(as.data.frame(t(tile_cov_gc_normalized_227)))
variance_raw = variance_raw %>% mutate(variance = variance_df$V1)
variance_raw = variance_raw %>% mutate(range = cov_ranges$range)

# Calculating IQR----------------------------------------------------------------------------
IQR_coverage = function(data){ #data here is tile_coverage matrix, with 110 rows
  time1 = system.time({ #we break data set to data sets with the size equal to step_size, then run them parallel.
    #processing big data self will be too time consuming.
    STEP_SIZE = 5000
    coverage_list = foreach(i = 0:(floor(dim(data)[2]/STEP_SIZE))) %dopar% {
      subdata = data[,(i*STEP_SIZE+1):(min((i+1)*STEP_SIZE,dim(data)[2]))]
      subdata %>% summarise_all(IQR)
    }
  })
  coverage_df = as.data.frame(t(do.call(cbind,coverage_list)))
  print(time1)
  return(coverage_df)
}


coverage_raw = IQR_coverage(as.data.frame(t(tile_cov_gc_normalized_227)))
colnames(coverage_raw) = "IQR"
#calculate the sigma distanse of each tile variance from the variance null model---------------------------------

variance_raw = variance_raw %>% mutate(var.sigma.dist = 
                                         ((variance_raw$variance - mean(variance_region_variances$variance))/sd(variance_region_variances$variance)))

#calculating the sigma dist of each tile coverage from the coverage null model--------------------------------------
#average over all coverage:
average_coverage = function(data){ #data here is tile_coverage matrix, with 110 rows
  time1 = system.time({ #we break data set to data sets with the size equal to step_size, then run them parallel.
    #processing big data self will be too time consuming.
    STEP_SIZE = 5000
    coverage_list = foreach(i = 0:(floor(dim(data)[2]/STEP_SIZE))) %dopar% {
      subdata = data[,(i*STEP_SIZE+1):(min((i+1)*STEP_SIZE,dim(data)[2]))]
      subdata %>% summarise_all(mean)
    }
  })
  coverage_df = as.data.frame(t(do.call(cbind,coverage_list)))
  print(time1)
  return(coverage_df)
}
a = average_coverage(as.data.frame(t(tile_cov_gc_normalized_227)))
coverage_raw = coverage_raw %>% mutate(average = (average_coverage(as.data.frame(t(tile_cov_gc_normalized_227)))))
coverage_raw = coverage_raw %>% mutate(cov.sigma.dist = 
                                        ((coverage_raw$V1 - mean(variance_region_coverage$coverage))/sd(variance_region_coverage$coverage)))

#KMedian distance----------------------------------------------
kmedian_distance_coverage = function(data){ #data here is tile_coverage matrix, with 110 rows
  time1 = system.time({ #we break data set to data sets with the size equal to step_size, then run them parallel.
    #processing big data self will be too time consuming.
    coverage_list = foreach(i = 1:length(data[i,])) %dopar% {
      subdata = (data[,i])
      if (length(unique(subdata)) > 2){
        k = kmeans(data[,i],centers = 2,nstart = 100)
        subdata = as.data.frame(data[,i]) %>% mutate(clusters = k$cluster)
        colnames(subdata) = c('cov',"clusters")
        if(min(k$size) > 1){
          coverage = (subdata %>% filter(clusters == 1))$cov
          peak1 = density(coverage)$x[which(density(coverage)$y == max(density(coverage)$y))]
          coverage = (subdata %>% filter(clusters == 2))$cov
          peak2 = density(coverage)$x[which(density(coverage)$y == max(density(coverage)$y))]
          return(c(peak1,peak2))
        }else{
          coverage = (subdata %>% 
            filter(clusters == which(table(subdata$clusters) == max(table(subdata$clusters)))))$cov
          peak1 = density(coverage)$x[which(density(coverage)$y == max(density(coverage)$y))]
          peak2 = (subdata %>% 
                     filter(clusters == which(table(subdata$clusters) == min(table(subdata$clusters)))))$cov
          return(c(peak1,peak2))
      }
      }else(
        return(c(0,0))
      )
    }
  })
  coverage_df = as.data.frame(t(do.call(cbind,coverage_list)))
  colnames(coverage_df) = c('peak1',"peak2")
  coverage_df = coverage_df %>% transmute(peak_dist = abs(peak2 - peak1))
  print(time1)
  return(coverage_df)
}

save.image("~/Codes/Germline_variation_detection/data_july9.RData")
kmedian_dist = kmedian_distance_coverage(as.data.frame(t(tile_cov_gc_normalized_227)))
save.image("~/Codes/Germline_variation_detection/data_july9.RData")
###################################################
a = rnorm(40000)
ggplot()+geom_density(aes(x = a))
k = kmeans(a,centers = 2, nstart = 500)
a = as.data.frame(a) %>% mutate(clusters = k$cluster)
ggplot(data = a)+geom_density(aes(x = cov,fill = as.factor(clusters)))
geom
#adding the new modal
c = rnorm(100000,mean = 2.5)
a = c(a,c)
###################################################
for (i in 1:200) {
  subdata = (data[,i])
  if (length(unique(subdata)) > 2){
    k = kmeans(data[,i],centers = 2,nstart = 100)
    subdata = as.data.frame(data[,i]) %>% mutate(clusters = k$cluster)
    colnames(subdata) = c('cov',"clusters")
    coverage = (subdata %>% filter(clusters == 1))$cov
    if(min(k$size) > 2){
    peak1 = density(coverage)$x[which(density(coverage)$y == max(density(coverage)$y))]
    coverage = (subdata %>% filter(clusters == 2))$cov
    peak2 = density(coverage)$x[which(density(coverage)$y == max(density(coverage)$y))]
    print(i)
    }
    else{
      
    }
  }else(
    print("here")
  )
}

###################################################

file$tile$n.cov.variance = variance_raw$variance
file$tile$n.cov.range = variance_raw$range
file$tile$n.IQR = coverage_raw$IQR
file$tile$n.var.sigma.dist = variance_raw$var.sigma.dist
file$tile$n.cov.sigma.dist = coverage_raw$cov.sigma.dist

