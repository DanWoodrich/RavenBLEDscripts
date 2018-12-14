#copy full moorings so don't have to do it by hand.

GTmoorings<-c("BS15_AU_02a","BS14_AU_04b","AW12_AU_BS3","BS13_AU_04b","BS16_AU_02a","BS15_AU_02b","AW14_AU_BS3")
  
drive1<-"G"
drive2<-"F"

####################

dir.create(paste(drive2,":/Full_datasets/",sep=""))

drive1<-paste(drive1,":/Waves/",sep="")
drive2<-paste(drive2,":/Full_datasets/",sep="")

copyvec<-dir(drive1)[which(dir(drive1)%in%GTmoorings)]

for(i in copyvec){
  dir.create(paste(drive2,i,sep=""))
  for(n in dir(paste(drive1,i,sep="/"))){
    file.copy(from=paste(drive1,i,"/",n,"/",dir(paste(drive1,i,"/",n,sep="/")),sep=""),to=paste(drive2,i,sep=""))
  }

}
