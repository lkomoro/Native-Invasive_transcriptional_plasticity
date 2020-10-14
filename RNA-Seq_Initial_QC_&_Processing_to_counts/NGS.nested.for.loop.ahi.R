
#setwd('~/lkomoroske/LK_for.loop.test/for.loop.cat.test')


parent <- '~/lkomoroske/LK_for.loop.test/for.loop.cat.test'
subdirs <- list.dirs(parent,full.names=TRUE,recursive=FALSE)#need to change recursive to false else it includes the parent directory as well in the list and messes up the length, etc.
outfile1 <- "combined_R1"
outfile2 <- "combined_R2"

for (i in 1:length(subdirs)){ 
  dir1 <- subdirs[i]
  sequences <- list.files(dir1,full.names=TRUE)#needs to be true so it has the full paths of the directories to go into
  R1s <- grep("R1", sequences)
  R2s <- grep("R2", sequences)
  
 	 for (k in 1:length(R1s)){
 	t <- scan(file = sequences[R1s[k]], what = 'character',sep="\r",quote="")
 	write(t, file= paste0(strsplit(dir1, '/')[[1]][6], "_",outfile1), append=TRUE)}
 
 	for (k2 in 1:length(R2s)){
	  t2 <- scan(file = sequences[R2s[k2]], what = 'character',sep="\r",quote="")
	write(t2, file= paste0(strsplit(dir1, '/')[[1]][6], "_",outfile2), append=TRUE)   
	} 
  
}