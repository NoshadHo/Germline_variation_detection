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

df = as.numeric(purified_tile_cov_gc_normalized[216564,])

#mad
mad = median(abs(df-median(df)))
sum((df > (median(input_vector)+mad*1.5)))
sum((df < (median(input_vector)-mad*1.5)))

#dbscan
db = fpc::dbscan(data = df,eps = (0.1),MinPts = 4,method = 'raw')
a = as.data.frame(df) %>% mutate(clust = db$cluster)
a %>% ggplot()+geom_point(aes(y = 0, x = df, color = as.factor(clust)))+geom_vline(xintercept = median(input_vector)+mad*1)+
  geom_vline(xintercept = -(median(input_vector)+mad*1))+theme_linedraw()+xlim(-1,1)

a %>% ggplot()+geom_point(aes(y = 0, x = df, color = as.factor(clust)))+geom_vline(xintercept = median(input_vector)+mad*1)+
  geom_vline(xintercept = -(median(input_vector)+mad*1))+theme_linedraw()

dbscan_tiles = foreach(i = 1:dim(purified_tile_cov_gc_normalized)[1]) %dopar% {
  db = fpc::dbscan(data = purified_tile_cov_gc_normalized[i,],eps = (0.05),MinPts = 5,method = 'raw')
  if(length(unique(db$cluster)) > 1){
    return(c(i,length(unique(db$cluster))))
  }
}
dbscan_tiles = as.data.frame(do.call(rbind,dbscan_tiles))
colnames(dbscan_tiles) = c('tile','clust_num')


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
cov_ranges %>% ggplot()+geom_histogram(aes(x = range),binwidth = 0.001)+xlim(0,5)+theme_linedraw()
cov_ranges %>% filter(range >13) %>% summarise(n())
cov_ranges = cov_ranges %>% mutate(tile = row_number())
cov_ranges %>% filter(range > 5)

#tiles ranges multimodal
k = Mclust(cov_ranges[,3],verbose = FALSE,G = 3)
summary(k)
cov_ranges_clust = (cov_ranges %>% select(range)) %>% mutate(cluster = k$classification) #add clusters to data.frame
cov_ranges_clust$cluster = as.factor(cov_ranges_clust$cluster)
cov_ranges_clust %>% ggplot(aes(x = range, y = cluster))+geom_point()  #look at the position of points in clusters

cov_ranges_clust %>% ggplot(aes(x = range))+geom_density(alpha = 0.4)+theme_linedraw() #look at the distribution of coverage in a specific tile as a whole

cov_ranges_clust %>% ggplot(aes(x = range, fill = cluster))+geom_density(alpha = 0.4)+theme_linedraw()+xlim(-1,5)+ylim(0,200)   #look at the distribution of coverage in a specific tile for each modal
range_dist_tiles = cov_ranges_clust %>% mutate(tile = row_number()) %>% filter(cluster == 9) %>% select(tile)
range_dist_tiles_2 = cov_ranges_clust %>% mutate(tile = row_number()) %>% filter(cluster == 9 | cluster == 8) %>% select(tile)


#Lets see how many of significant tiles are captured in other methods
multimodal_tiles = significant_tiles %>% select(tile)                   #tile multimodal
range_tiles = cov_ranges %>% filter(range > 0.5) %>% select(tile)       #tile range
range_dist_tiles
dbscan_tiles_2 = dbscan_tiles %>% filter(clust_num == 2) %>% select(tile)
dbscan_tiles_3 = dbscan_tiles %>% filter(clust_num == 3) %>% select(tile)
dbscan_tiles_all = dbscan_tiles %>% select(tile)

dim(intersect(range_dist_tiles_3,range_dist_tiles_2))

#remove intersect of tile_ranges and tile_Ranges_multimodal, see what are the rest of the selected
excluded_trmultimodal_tr = range_dist_tiles %>% filter(!(tile %in% range_tiles$tile))
excluded_trmultimodal_tr = dbscan_tiles %>% anti_join(range_dist_tiles,by = "tile")
excluded_trmultimodal_tr_2 = dbscan_tiles %>% select(tile) %>% anti_join(a,by = "tile")
