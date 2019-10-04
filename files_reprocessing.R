## Script to rerun cnvex samples

# GOAL IS TO CREAT A FILE CONTAINING THIS(BELOW) FOR WHOLE COHORT

# Rscript ./pipe/cords-cnvex/test_process.R \
# -i ~/Data/tmp/SI_12071.SI_12072.rds \
# -o ~/Data/tmp/agiv4_SI8989_8988_1.rds \
# -l ~/Data/WXS_prostate/results/agi.v4.pca.pool.rds \
# -s exome \
# -p 30 \
# -x male

## VARIABLES
SETTING = "genome"
CORES = 10
POOL = "~/Data/cnvex_pool/results/agi.v4.pca.V2.pool.rds"



## files addresses $$$$$$$$$$ this whole part needs modification for different cohorts
kidney_list = read_delim(delim = "\n", file = "/data/kidney_cancer_cohort_list.txt", col_names = F)
kidney_list$X1 = substr(kidney_list$X1, start = 1, stop = 9)

files = list.files("/data/WGS/cnvex",full.names = T)
kidney_files = files[files %in% paste0("/data/WGS/cnvex/",kidney_list$X1)]
kidney_files = list.files(kidney_files, pattern = paste0("\\d.rds$") ,full.names = T)
ilist = kidney_files
olist = paste0(substr(ilist,start = 1, stop = 10),"kidney_rerun1/", substr(ilist, start = 17, stop = nchar(ilist)))   


#agi.v4
agi.v4.files = read_delim(delim = "\n", file = "/data/cohorts/agi.v4.samples.txt", col_names = F)
ilist = agi.v4.files$X1
olist = paste0("/data/WXS/rerun1", substr(ilist, start = 18, stop = nchar(ilist)))   ###NOT WORKING

## FIRST LINE
line1 = sapply(X = 1:length(ilist), function(X) {
  return("Rscript ./pipe/cords-cnvex/test_process.R \\")
})
## SECOND LINE
line2 = sapply(X = (ilist), function(X) {
  string = paste0("-i ", X, " \\")
})
## THIRD LINE
line3 = sapply(X = (olist), function(X) {
  string = paste0("-o ", X, " \\")
})
## FORTH LINE
line4 = sapply(X = 1:length(ilist), function(X) {
  string = paste0("-s ", SETTING, " \\")
})
## FIFTH LINE
line5 = sapply(X = 1:length(ilist), function(X) {
  string = paste0("-p ", CORES, " \\")
})
# ## SIXTH LINE
# line6 = sapply(X = 1:length(ilist), function(X) {
#   string = paste0("-x ", "male", " \\")
# })
## SEVENTH LINE
line7 = sapply(X = 1:length(ilist), function(X) {
  string = paste0("-l ", POOL)
})
## EIGTH LINE
line8 = sapply(X = 1:length(ilist), function(X) {
  string = paste0("echo ", X, " >> /data/progress.txt"," \n")
})

## Making directories
system("for i in `ls /data/cords-cnvex/`; do mkdir /data/WXS/rerun1/$i; done")
# system2("mkdir",substr(X, start = 18, stop = length(X)))

pre_commands = paste("#!/usr/bin/bash","ls -l .*\n","rm /data/progress.txt","touch /data/progress.txt\n",sep = "\n")
commands = paste(line1,line2,line3,line4,line5,line7, sep = "")
commands = paste(commands, line8, sep = "\n")
write(pre_commands,"/data/commands/kidney_cohort_proocess.sh")
write(commands,"/data/commands/kidney_cohort_proocess.sh",append = T)

## in case you want it to look bettter, put sep = "\n"

## runing the file in parallel :
system("parallel -j 10 :::: /data/commands/agi.v4_process.sh")
