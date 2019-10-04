#adding segment to tile
seg_tiles <- queryHits(findOverlaps(cnv_old$tile,cnv_old$seg))
seg_ids <- subjectHits(findOverlaps(cnv_old$tile,cnv_old$seg))
cnv_old$tile$seg <- seg_ids

seg_tiles <- queryHits(findOverlaps(cnv_new$tile,cnv_new$seg))
seg_ids <- subjectHits(findOverlaps(cnv_new$tile,cnv_new$seg))
cnv_new$tile$seg <- seg_ids


#calculate the sd around each segment
old_sd = as.data.frame(cnv_old$tile) %>% group_by(seg) %>% summarise(sd = sd(lr,na.rm = T))
new_sd = as.data.frame(cnv_new$tile) %>% group_by(seg) %>% summarise(sd = sd(lr,na.rm = T))
mean(new_sd$sd,na.rm = T)
mean(old_sd$sd,na.rm = T)

#for all the samples:
agi.v4.files = read_delim(delim = "\n", file = "/data/cohorts/agi.v4.samples.txt", col_names = F)
olds = agi.v4.files$X1
news = paste0("/data/WXS/rerun1", substr(olds, start = 18, stop = nchar(olds)))

list = foreach(i = 1:length(olds)) %dopar% {
  cnv_old = read_rds(olds[i])
  cnv_new = read_rds(news[i])
  
  seg_tiles <- queryHits(findOverlaps(cnv_old$tile,cnv_old$seg))
  seg_ids <- subjectHits(findOverlaps(cnv_old$tile,cnv_old$seg))
  cnv_old$tile$seg <- seg_ids
  
  seg_tiles <- queryHits(findOverlaps(cnv_new$tile,cnv_new$seg))
  seg_ids <- subjectHits(findOverlaps(cnv_new$tile,cnv_new$seg))
  cnv_new$tile$seg <- seg_ids
  
  old_sd = as.data.frame(cnv_old$tile) %>% group_by(seg) %>% summarise(sd = sd(lr,na.rm = T))
  new_sd = as.data.frame(cnv_new$tile) %>% group_by(seg) %>% summarise(sd = sd(lr,na.rm = T))
  new = median(new_sd$sd,na.rm = T)
  old = median(old_sd$sd,na.rm = T)
  cat(as.character(i),file="/data/progress.txt",sep="\n",append=TRUE)
  return(c(old,new))
}
seg_medians = do.call(rbind,list)
colnames(seg_means) = c("old","new")
# these id's, sd has been increased
# [1]  49 118 127 130 133 134 136 159 160 179 185 191 193 200 207 213 216 217 218 219 220 221 222 223 224 225 226 227 228 230 231 239 240 241 244 246 249 250 252 253
# [41] 254 262 264 265 266 267 268 269 271 272 273 274 277 279 280 284 285 287 291 297 298 300 313 317 336 337 339 341 342 363

#Box plot
box_data = as.data.frame(seg_means) %>% gather()
box = box_data %>% ggplot() + geom_boxplot(aes(x = key, y = value), fill = c("#E69F00", "#56B4E9"))+scale_fill_brewer(palette="Dark2")
