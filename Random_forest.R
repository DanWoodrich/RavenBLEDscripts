install.packages("randomForest")
install.packages("seewave")
install.packages("tuneR")

library(randomForest)
library(seewave)
library(tuneR)



#Which data would you like to evaluate?
#species
spec<-"RW"

runname<-



##########
detfiles<-list.files(paste("E:/DetectorRunOutput/",runname,sep=""),pattern = "TPFPFN")  

mooringpat=NULL
for(n in 1:length(detfiles)){
  mooringpat<-c(mooringpat,substr(detfiles[1],1,12))
}

#Can automate this by digging out parameters, but not worth it right now
#was data whitened?(y or no)
whiten<-"y"
if(whiten=="y"){
  FO<-100 #filter order
  LMS<-.10 #LMS step size
  soundfile<-list.files(paste("E:/Combined_sound_files/",spec,sep=""),pattern =paste(LMS*100,"x_FO",FO,sep=""))
  soundfiles<-NULL
  for(n in 1:length(dir(paste("E:/Combined_sound_files/",spec,"/",soundfile,sep=""))))
  dir(paste("E:/Combined_sound_files/",spec,"/",soundfile,sep=""))
}else
  soundfiles<-dir("E:/Combined_sound_files/RW/No_whiten")
}
  
  
data=NULL
for(n in 1:length(detfiles)){
  data2<-read.csv(paste("E:/DetectorRunOutput/",runname,detfiles[n],sep=""))
  data2<-cbind(soundfiles[n])
  data<-rbind(data,data2)
}

#translate response to binary 
data[which(data$detectionType=="TP"),8]<-1
data[which(data$detectionType=="FP"),8]<-0
data[which(data$detectionType=="FN"),8]<-2

#remove FN from data
data<-data[which(data$detectionType==0|data$detectionType==1),]



readWave()
spec()