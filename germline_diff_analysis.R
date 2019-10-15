library(dbscan)
### First analysis: different between Pair and Nvar analysis, what are the different tile? how much?

pair_cnv = read_rds("/mctp/users/noshadh/data/WGS/cnvex/C3L-00144/C3L-00144.rds")
nvar_cnv = read_rds("~/Data/C3L-00144.rds")
nvar_cnv = addSeg(nvar_cnv)
dummy_cnv = nvar_cnv

lr.diff = pair_cnv$tile$lr - nvar_cnv$tile$lr
ggplot()+geom_point(aes(x = 1:length(lr.diff),y = lr.diff, color = as.factor(nvar_cnv$tile$blacklist>0)), size = 0.5)

## find the outliers using dbscan
df = as_tibble(matrix(ncol = 2,c(lr.diff,1:length(lr.diff))))
colnames(df) = c("y","x")
df$y[is.na(df$y)] = 0
# db = dbscan(df,eps = 2.003,minPts = 5)

q1 = quantile(na.omit(lr.diff))[2]
q3 = quantile(na.omit(lr.diff))[4]

df = df %>% mutate(cluster = ifelse(y > q1*4.5 & y < q3*4,1,0))

# df = df %>% mutate(cluster = db$cluster)
df %>% ggplot()+geom_point(aes(x = x,y = y, color = as.factor(cluster)))+theme(legend.position = "none")

df = df %>% arrange(desc(abs(y))) %>% mutate(cluster = ifelse(row_number() > 2500,0,1))


##########RUN IT FOR EVERY SAMPLE
cnv.fns = read_delim("/mctp/users/noshadh/data/cohorts/kidney_cohort.txt",delim = "\n",col_names = F)
cnv.fns = paste0("/mctp/users/noshadh",cnv.fns$X1)
cnv.fns.new = str_replace(cnv.fns,"cnvex","kidney_rerun1")
res = foreach(i = 1:length(cnv.fns)) %dopar% {
  pair_cnv = read_rds(cnv.fns[i])
  pool_cnv = read_rds(cnv.fns.new[i])
  
  lr.diff = pair_cnv$tile$lr - pool_cnv$tile$lr
  
  df = as_tibble(matrix(ncol = 2,c(lr.diff,1:length(lr.diff))))
  colnames(df) = c("y","x")
  df$y[is.na(df$y)] = 0
  
  df = df %>% arrange(desc(abs(y))) %>% mutate(cluster = ifelse(row_number() > 2000,0,1)) %>% arrange(x) #find cluster (select top 2500 with most diff)
  df = df %>% mutate(y = ifelse(cluster == 1, y, NA_real_))
  
  return(df$y)
}
diffs = as.data.frame(do.call(cbind,res))
colnames(diffs) = paste0("s",1:ncol(diffs))
diffs_gather = diffs %>% gather()
diffs_gather = diffs_gather %>% mutate(x = rep(1:nrow(diffs),ncol(diffs)))
diffs_gather %>% ggplot()+geom_point(aes(x = x,y = value, color = key))+theme(legend.position = "none")+ggtitle("Pair-Pool Germlines")+theme_linedraw()

##only autosoms
na.omit(diffs_gather) %>% filter(x < 287510) %>% ggplot()+geom_histogram(aes(x = x,fill = ..count.. > 600),bins = 400)+theme_linedraw()
