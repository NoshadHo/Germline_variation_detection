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
db = fpc::dbscan(data = df,eps = (0.1),MinPts = 3,method = 'raw')
a = as.data.frame(df) %>% mutate(clust = db$cluster)
a
a %>% ggplot()+geom_point(aes(y = 0, x = df, color = as.factor(clust)))+geom_vline(xintercept = median(input_vector)+mad*1)+
  geom_vline(xintercept = -(median(input_vector)+mad*1))+theme_linedraw()+xlim(-1,1)

