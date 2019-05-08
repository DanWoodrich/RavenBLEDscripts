#script for subsetting data that is not labelled. Protocol is run an initial detector, take any PNG that features a detection. 

Mooring<-"BS13_AU_02a"
Month<-"July"
HGdrivePath<-"E:/"
FDdrivePath<-"G:/"
SelecPath<-"//nmfs/akc-nmml/CAEP/Acoustics/Projects/Dans Detectors/LF moan project/SelectionTabelsInitialRun"

#identify files in mooring
sound_files<-c()
sound_filesfullpath<-c()
allMonths<-dir(paste(FDdrivePath,"/",Mooring,sep=""))
for(i in allMonths){
  sound_files<-c(sound_files,dir(paste(FDdrivePath,"/",Mooring,"/",i,sep=""),pattern=".wav"))
  sound_filesfullpath<-c(sound_filesfullpath,paste(paste(FDdrivePath,"/",Mooring,sep=""),i,dir(paste(FDdrivePath,"/",Mooring,"/",i,sep=""),pattern=".wav"),sep="/"))
}
sound_files <- sound_files[order(sound_files)] #need to look at how mooring is structured but should work fine for sox with a list of full path files. 
sound_filesfullpath<-sound_filesfullpath[order(sound_filesfullpath)]

#identify files in selection table 

Selections<-read.delim(paste(SelecPath,"/",Mooring,"/",Month,".txt",sep=""))

sound_filesSelections<-unique(Selections$End.File)

sound_files_same<-which(sound_files %in% sound_filesSelections)

destFile<-paste("E:/HG_datasets/",Mooring,"/LM_yesUnion",sep="")

dir.create(destFile)

file.copy(sound_filesfullpath[sound_files_same], destFile)


#some data analysis on the first detector results: 

TPFPFN<-read.delim("E:/DetectorRunOutput/test LM_20190507102454/LM_BS13_AU_02a_files_38-122_TPFPFN_Tab_Ravenformat.txt")

Labeled<-read.csv("E:/DetectorRunFiles/LMProcessed_GT_data/test LM_20190507102454_processedGT.csv")

Labeled$duration<-Labeled$End.Time..s.-Labeled$Begin.Time..s.
TPFPFN$duration<-TPFPFN$End.Time..s.-TPFPFN$Begin.Time..s.


ggplot(Labeled, aes(duration, fill = as.factor(detectionType))) + geom_density(alpha = 0.2)

ggplot(TPFPFN, aes(duration, fill = as.factor(detectionType))) + geom_density(alpha = 0.2)
