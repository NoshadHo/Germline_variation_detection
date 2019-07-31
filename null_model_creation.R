# steps:

#choose the regions
chr = 22
start = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                     ungroup() %>% slice(chr) %>% select(tile))
end = as.numeric((as.data.frame(file$tile) %>% mutate(tile = row_number()) %>% group_by(seqnames) %>% slice(1))%>% 
                   ungroup() %>% slice(chr+1) %>% select(tile))

#variance raw is being calculated in blacklist_validation_plots.R
variance_raw %>% dplyr::slice(start:end) %>% ggplot()+
  geom_point(aes(x = start:end,y = variance,color = blacklist),size = 0.3)+
 # coord_cartesian(ylim = c(0,2))+
  scale_x_continuous(breaks=seq(start,end,1000))

#load the regions
variance_regions = c(2001:10001, 28900:32900, 50000:56200, 70000:72400, 90000:92000, 107000:111300, 124000:128500, 140600:142700, 154000:157100, 167550:171000, 181200:185800, 194400:197400, 209710:216000, 222100:229100, 233000:239800, 240150:241300, 249200:250600, 257800:258800, 266450:267500, 271400:273800, 278000:282000, 285400:287400,
                     16001:22001, 36900:43900, 62100:65120, 74500:87500, 96500:105500, 112500:120100, 130000:133200, 144100:147200, 160700:167100, 172500:174500, 187000:188800, 198400:200300, 216300:218700,245100:248500, 254000:257100, 259800:261900, 274700:277400,
                     44900:48000, 121600:122500, 133700:137300, 148000:153140, 176000:180500, 189800:193800, 200900:204380,262400:265000,
                     137700:138700, 204700:207300)

#calculate the variance of those regions(tiles)
variance_raw = variance_raw %>% mutate(tile = row_number())     #variance_row is calculated in blacklist_validation_plots
variance_region_variances = variance_raw %>% filter(tile %in% variance_regions)
colnames(variance_region_variances) = c('variance','tile')
#calculate the density of those variances
variance_region_variances %>% ggplot(aes(x = variance)) +
  geom_histogram(binwidth = 0.0001,aes(fill = (variance > 0.00012000 & variance < 0.018)))+
  coord_cartesian(xlim = c(0,0.1))+theme_linedraw()+
  ggtitle("variance null model: mean = 0.0137675, sd = 0.006562863, 1st Q=0.0127366, 2nd Q = 0.0132257 , 3rd Q = 0.0137781")
mean(variance_region_variances$variance)
sd(variance_region_variances$variance)
summary(variance_region_variances$variance)



#COVERAGE NULL MODEL -----------------------------------------------------------

coverage_raw = coverage_raw %>% mutate(tile = row_number())     #variance_row is calculated in blacklist_validation_plots
variance_region_coverage = coverage_raw %>% filter(tile %in% variance_regions)
colnames(variance_region_coverage) = c('coverage','tile')

#calculate the density of those coverage
variance_region_coverage %>% ggplot(aes(x = coverage)) +
  geom_histogram(binwidth = 0.001,aes(fill = (coverage > 0.3116 & coverage < 0.3171)))+
  coord_cartesian(xlim = c(0,0.6))+theme_linedraw()+
  ggtitle("coverage null model: mean = 0.3155561, sd = 0.1287049, 1st Q=0.3116, 2nd Q = 0.3145 , 3rd Q = 0.3171")
mean(variance_region_coverage$coverage)
sd(variance_region_coverage$coverage)
summary(variance_region_coverage$coverage)

#LR NULL MODEL --------------------------------------------------------------------
lr_raw = tile_lr %>% select(lr = !!25)
lr_raw = lr_raw %>% mutate(tile = row_number())     #variance_row is calculated in blacklist_validation_plots
colnames(variance_region_coverage) = c('lr','tile')
lr_null_regions = lr_raw %>% filter(tile %in% variance_regions)

