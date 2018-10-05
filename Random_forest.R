install.packages("randomForest")
install.packages("seewave")
install.packages("tuneR")

library(randomForest)
library(seewave)
library(tuneR)



#Which data would you like to evaluate?
#species
samplingRate<-16384 #hz
yminn<-0 #for spec plotting. Should be the same as detector window preset
ymaxx<-1000 #" "

spec<-"RW"

runname<-



##########
detfiles<-list.files(paste("E:/DetectorRunOutput/",runname,sep=""),pattern = "TPFPFN")  

#extract mooring names from moorings used in run
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
}else{
  soundfile<-dir("E:/Combined_sound_files/RW/No_whiten")
}

#only choose soundfiles that match those used in run
for(n in 1:length(dir(paste("E:/Combined_sound_files/",spec,"/",soundfile,sep="")))){
  if(dir(paste("E:/Combined_sound_files/",spec,"/",soundfile,sep=""))[n] %in% mooringpat){
    soundfiles<-c(soundfiles,dir(paste("E:/Combined_sound_files/",spec,"/",soundfile,sep=""))[n])
  }
}

#make sure they are in same order
soundfiles<-sort(soundfiles)
detfiles<-sort(detfiles)
  
  
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

for(n in 1:nrow(data)){
  clipwav <- readWave(data[n,1],data[n,4],data[n,5],units="seconds")
  clipspec = spec(clipwav, plot=F, wl=samplingRate, PSD=T,ylim=c(yminn,ymaxx))
  
}

readWave()
spec()