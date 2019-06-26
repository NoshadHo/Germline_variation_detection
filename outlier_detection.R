cppFunction('
  double median2(std::vector<double> x){
    double median;
    size_t size = x.size();
    sort(x.begin(), x.end());
    if (size  % 2 == 0){
        median = (x[size / 2 - 1] + x[size / 2]) / 2.0;
    }
    else {
        median = x[size / 2];
    }
    return median;
}
')

df = as.numeric(purified_tile_cov_gc_normalized[186436,])

#mad
mad = median(abs(df-median(df)))
sum((df > (median(input_vector)+mad*1.5)))
sum((df < (median(input_vector)-mad*1.5)))

#dbscan
db = fpc::dbscan(data = df,eps = (0.12),MinPts = 10,method = 'raw')
a = as.data.frame(df) %>% mutate(clust = db$cluster)
a
a %>% ggplot()+geom_point(aes(y = 0, x = df, color = as.factor(clust)))+geom_vline(xintercept = median(input_vector)+mad*1)+
  geom_vline(xintercept = -(median(input_vector)+mad*1))+theme_linedraw()+xlim(-1,1)


#finding the max and min for all the tiles:
TIME = system.time({
  cov_ranges = foreach(i = 1:dim(purified_tile_cov_gc_normalized)[1]) %dopar%{
    row = purified_tile_cov_gc_normalized[i,]
    return(c(min(row),max(row)))
  }
})
cov_ranges = as.data.frame(do.call(rbind,cov_ranges))
colnames(cov_ranges) = c('min','max')
cov_ranges = cov_ranges %>% mutate(range = max - min)
cov_ranges %>% ggplot()+geom_histogram(aes(x = range),binwidth = 0.001)+xlim(0,1.5)+theme_linedraw()
cov_ranges %>% filter(range >7) %>% summarise(n())


db = fpc::dbscan(data = cov_ranges$range,eps = (0.1),MinPts = 10,method = 'raw')
