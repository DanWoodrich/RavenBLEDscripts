library(imager)
biocLite("EBImage")
library(EBImage)
#install_keras(conda = "C:/Users/daniel.woodrich/AppData/Local/Continuum/anaconda3/Scripts/conda.exe")
species<-"GS"

drivepath<-"E:/"
libpath<-paste(drivepath,"/DetectorRunFiles/Image_library/No_whiten_decimate_by_8/GS",sep="")

labelledY<-vector("list",length(dir(paste(libpath,"Yes",sep="/"))))
labelledN<-list()
#load in images to training and test dataset
for(i in 1:length(dir(paste(libpath,"Yes",sep="/")))){
  labelledY[[i]]<-readImage(paste(libpath,"Yes",dir(paste(libpath,"Yes",sep="/")),sep="/")[i])
}

for(i in 1:length(dir(paste(libpath,"No",sep="/")))){
  labelledN[[i]]<-readImage(paste(libpath,"No",dir(paste(libpath,"No",sep="/")),sep="/")[i])
}


testData <- extract_feature("test1/", width, height, labelsExist = F)