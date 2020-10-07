rm(list = ls())
#setwd('~/Project_Fangue_all_R1_files')#can't change working directory through script when 
#workin on ahi-do command line so are where you want to be, then run script
#getwd()
file.list <- list.files(pattern = "R1")
print(file.list)

for (i in file.list) {
  new.name <- paste(strsplit(i, c("-", "."))[[1]][2], ".fastq", sep = "")
  t <- scan(i, what = 'character',sep="\r",quote="")
  write(t, file = new.name, append=TRUE)
   unlink(i) #include this if want it to delete original files
}
