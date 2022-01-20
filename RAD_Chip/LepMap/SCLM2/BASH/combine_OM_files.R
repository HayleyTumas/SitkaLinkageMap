output <- #Path to Output
datafiles <- #Path to DataFiles

require(readxl)
require(tidyverse)

PATH=psate(output , "Top_OM", sep="")
#create a list of the files from your target directory
file_list <- list.files(path=PATH)

rt <- function(x){
  read.table(paste(PATH,"/", x, sep=""), header=F)
}

myfiles <- lapply(file_list, rt)
names(myfiles) <- sub(".txt", "", file_list)
rows <- lapply(myfiles, nrow)
for(i in 1:12){
  myfiles[[i]]$LG <- c(rep(paste("LG",sub("_", "",regmatches(file_list[[i]], regexpr("[0-9].", file_list[[i]]))), sep=""), rows[[i]]))
}
full <- do.call("rbind", myfiles)

subset <- full[,c(1:2, 13)]
colnames(subset) <- c("SNP","Position","LG")
write.csv(subset, file=paste(output, "RADSNP_SCLM2_Map.csv", sep=""))

#add in SNP IDs
ID <- read.table(file=paste(datafiles, "RADSNP_SNPIDs.txt", sep=""))
colnames(ID) <- c("CHR","POS","SNPID")
CallID <- read.table(file=paste(datafiles, "RADSNP_callfile.txt", sep=""), -1)
CallID$Order <- c(1:nrow(CallID))
withID <- merge(subset, CallID, by.x="SNP", by.y="Order")
withID$link <- paste(withID$CHR, withID$POS, sep=".")
ID$link <- paste(ID$CHR, ID$POS, sep=".")
withSNPID <- merge(withID, ID, by="link")
write.csv(withSNPID, file=paste(output, "RADSNP_SCLM2_Map_SNPID.csv", sep=""))
