#install.packages("flightcallr", repos="http://R-Forge.R-project.org")
#install.packages('seewave')
#install.packages('tuneR')
#install.packages('randomForest')
#install.packages('ROCR')



library('flightcallr')
library('seewave')


library("dplyr") 

Mooring <- "BS15_AU_02a"

SelectTabAll <- list.files(paste('E:/Selection tables/',Mooring,sep=""))
dir.create(paste('E:/Selection tables/',Mooring,"Sum/",sep=""))

SelectTabComb <- NULL

for(i in SelectTabAll){
  var1 <- read.delim(paste('E:/Selection tables/',Mooring,'/',i,sep=""))
  SelectTabComb <- rbind(SelectTabComb,var1)
}
OG_names <- colnames(SelectTabComb)
colnames(SelectTabComb)<-c("Selection","View","Channel","Begin Time (s)","End Time (s)","Low Freq (Hz)","High Freq (Hz)")
write.table(SelectTabComb,paste('E:/Selection tables/',Mooring,"Sum/",Mooring,"_All.txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
  
SelectTabComb <- SelectTabComb[SelectTabComb$View=="Spectrogram 1",]
colnames(SelectTabComb) <- OG_names

STC_meanH<- mean(SelectTabComb$High.Freq..Hz.)
STC_varH<-sd(SelectTabComb$High.Freq..Hz.)
STC_meanL<- mean(SelectTabComb$Low.Freq..Hz.)
STC_varL<-sd(SelectTabComb$Low.Freq..Hz.)
SelectTabComb$Range<- (SelectTabComb$High.Freq..Hz.-SelectTabComb$Low.Freq..Hz.)
STC_meanR<- mean(SelectTabComb$Range)
STC_varR<- sd(SelectTabComb$Range)
SelectTabComb$Duration<- (SelectTabComb$End.Time..s.-SelectTabComb$Begin.Time..s.)
STC_meanD<- mean(SelectTabComb$Duration)
STC_varD<- sd(SelectTabComb$Duration)

hist(SelectTabComb$Range)
hist(SelectTabComb$Duration)
hist(SelectTabComb$High.Freq..Hz.)
hist(SelectTabComb$Low.Freq..Hz.)




Summary <- rbind(c(Mooring,"HighFreq","LowFreq","RangeFreq","DurationTime"),c("Mean",STC_meanH,STC_meanL,STC_meanR,STC_meanD),c("SD",STC_varH,STC_varL,STC_varR,STC_varD))
write.table(Summary,paste('E:/Selection tables/',Mooring,"Sum/",Mooring,"_Summary.txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE,col.names=FALSE)

                 
#278 boxes from 325 files BS16_AU_02a. Why? 1. harder to see RW upcalls on Raven than on Soundchecker 2. Error in data clipping script including maybes 3. Amplified observed analyst error due to crypic nature of calls (seems like would be caught in review)
#^ It's because I was doing the paging wrong. Only seeing 1st 30 seconds of each sound file. This run was a subsampled review. 

list.files(system.file(package = 'flightcallr'), recursive = T, full.names = T) 
system.file()
?flightcallr
data(danby)
