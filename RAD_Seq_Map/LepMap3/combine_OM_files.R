require(readxl)
require(tidyverse)

output <- #path to Output/
  datafiles=#path to DataFiles/

PATH=#Output/Top_OM
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
write.csv(subset, file=paste(output, "RADQC_2FAMSS_Map.csv", sep=""))

#add in SNP IDs
ID <- read.table(file=paste(output, "2FAMSS_SNPIDs.txt", sep=""))
colnames(ID) <- c("CHR","POS","SNPID")
CallID <- read.table(file=paste(datafiles, "2FAMSS_callfile.txt", sep=""), -1)
CallID$Order <- c(1:nrow(CallID))
withID <- merge(subset, CallID, by.x="SNP", by.y="Order")
withID$link <- paste(withID$CHR, withID$POS, sep=".")
ID$link <- paste(ID$CHR, ID$POS, sep=".")
withSNPID <- merge(withID, ID, by="link")
write.csv(withSNPID, file=paste(output, "RADQC_2FAMSSSCLM2_Map_SNPID.csv", sep=""))
