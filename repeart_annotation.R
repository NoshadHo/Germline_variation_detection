repeat_file = read_tsv("~/Codes/Germline_variation_detection/genome_repeat_5.bed",col_names = TRUE)
blacklist_file = read_tsv("~/Codes/Germline_variation_detection/BLACKLIST.bed",col_names = FALSE)
colnames(blacklist_file) = c("seqnames","start", "end")
#generate a file for bed_to_tile function:
repeat_file_short = repeat_file %>% select(genoName,genoStart,genoEnd,repClass)
repeat_file_short = repeat_file_short %>% filter(nchar(genoName) == 4 | nchar(genoName) == 5)
colnames(repeat_file_short) = c('seqnames','start','end','repeat_class')  
repeat_file_short = repeat_file_short %>% arrange(seqnames,start)

#this code can be easily converted to the 'per region' annotation 
blacklist_length = as.numeric(blacklist_file %>% transmute(length = end-start) %>% summarise(sum(length)))
repeat_SINE = 0
repeat_LINE = 0
repeat_LTR = 0
repeat_SIMPLE = 0
repeat_proportion = 0

for (i in 1:dim(blacklist_file)[1]){
  chr = as.character(blacklist_file[i,1])
  start_bl = as.numeric(blacklist_file[i,2])
  end_bl = as.numeric(blacklist_file[i,3])
  repeats = repeat_file_short %>% filter(seqnames == chr & start > start_bl & end < end_bl)
  repeat_proportion = repeat_proportion + as.numeric(repeats %>% select(start,end) %>% summarise_all(sum) %>% transmute(end-start))
    #/(end_bl-start_bl)
  repeats = repeats %>% group_by(repeat_class) %>% summarise(count = n()) %>% 
    mutate(weighted_percent = (count/sum(count))*(end_bl-start_bl)/blacklist_length)
  if (!is.na(as.numeric(repeats %>%  filter(repeat_class == 'SINE') %>% select(weighted_percent)))){
  repeat_SINE = repeat_SINE + as.numeric(repeats %>%  filter(repeat_class == 'SINE') %>% select(weighted_percent))}
  
  if (!is.na(as.numeric(repeats %>%  filter(repeat_class == 'LINE') %>% select(weighted_percent)))){
  repeat_LINE = repeat_LINE + as.numeric(repeats %>%  filter(repeat_class == 'LINE') %>% select(weighted_percent))}
  
  if (!is.na(as.numeric(repeats %>%  filter(repeat_class == 'LTR') %>% select(weighted_percent)))){
  repeat_LTR = repeat_LTR + as.numeric(repeats %>%  filter(repeat_class == 'LTR') %>% select(weighted_percent))}
  
  if (!is.na(as.numeric(repeats %>%  filter(repeat_class == 'Simple_repeat') %>% select(weighted_percent)))){
  repeat_SIMPLE = repeat_SIMPLE + as.numeric(repeats %>%  filter(repeat_class == 'Simple_repeat') %>% select(weighted_percent))}
  print(i)
}
repeat_proportion = repeat_proportion/blacklist_length

# Pie Chart with Percentages
slices = c(1-repeat_proportion, 
            repeat_SINE*repeat_proportion, 
            repeat_LINE*repeat_proportion, 
            repeat_LTR*repeat_proportion,
            repeat_SIMPLE*repeat_proportion,
            (1 - (repeat_SINE+repeat_LINE+repeat_LTR+repeat_SIMPLE))*repeat_proportion)
lbls <- c("NO REPEAT", "SINE", "LINE", "LTR", "SIMPLE REPEATS", "OTHERS")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=terrain.colors(length(lbls)),
    main="Pie Chart of repeat annotations") 


#important ones: Simple_repeat,SINE,LINE,Satellite,LTR
#
# repeat_class  
# <chr>         
# 1 Simple_repeat 
# 2 Satellite     
# 3 LINE          
# 4 DNA           
# 5 SINE          
# 6 LTR           
# 7 Low_complexity
# 8 LTR?          
# 9 snRNA         
# 10 tRNA          
# 11 DNA?          
# 12 Retroposon    
# 13 srpRNA        
# 14 rRNA          
# 15 Unknown       
# 16 RC            
# 17 scRNA         
# 18 RNA           
# 19 RC?           
# 20 SINE? 
#