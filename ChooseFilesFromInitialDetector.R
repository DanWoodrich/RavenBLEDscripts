#script for subsetting data that is not labelled. Protocol is run an initial detector, take any PNG that features a detection. 

Mooring<-"AL16_AU_BS3"
Month<-"January"
HGdrivePath<-"E:/"
FDdrivePath<-"E:/"
SelecPath<-"//nmfs/akc-nmml/CAEP/Acoustics/Projects/Dans Detectors/LF moan project/SelectionTabelsInitialRun"

#identify files in mooring
sound_files<-c()
sound_filesfullpath<-c()
allMonths<-dir(paste(FDdrivePath,"/Full_datasets/",Mooring,sep=""))
for(i in allMonths){
  sound_files<-c(sound_files,dir(paste(FDdrivePath,"/",Mooring,"/",i,sep=""),pattern=".wav"))
  sound_filesfullpath<-c(sound_filesfullpath,paste(paste(FDdrivePath,"/",Mooring,sep=""),i,dir(paste(FDdrivePath,"/",Mooring,"/",i,sep=""),pattern=".wav"),sep="/"))
}
sound_files <- sound_files[order(sound_files)] #need to look at how mooring is structured but should work fine for sox with a list of full path files. 
sound_filesfullpath<-sound_filesfullpath[order(sound_filesfullpath)]

#identify files in selection table 

Selections<-read.delim(paste(SelecPath,"/",Mooring,"/",Month,".txt",sep=""))

#Selections<-read.delim("E:/DetectorRunOutput/LM AL16_AU_BS3_20190513174739/01AL16_AU_BS3_files_All_AllSpecies_Model_Applied_probs.txt")

sound_filesSelections<-unique(Selections$End.File)

#sound_filesSelections<-unique(Selections$File)

sound_files<-dir("E:/Full_datasets/AL16_AU_BS3/AL16_AU_BS3_files_All_decimate_by_128")
sound_filesfullpath<-paste("E:/Full_datasets/AL16_AU_BS3/AL16_AU_BS3_files_All_decimate_by_128",sound_files,sep="/")

sound_files_same<-which(sound_files %in% sound_filesSelections)

destFile<-paste("E:/HG_datasets/",Mooring,"/LM_yesUnion",sep="")

destFile<-"E:/HG_datasets/AL16_AU_BS3/LM_yesUnion/AL16_AU_BS3_files_All_decimate_by_128"

dir.create(destFile)

file.copy(sound_filesfullpath[sound_files_same], destFile)

file.copy(sound_filesSelections)

#load mooring with probabilities and use it to 'high grade'
par(mgp=c(2.5,1,0)) 
par(cex.main=2,cex.lab=2,cex.axis=1.5)
Selections<-read.delim("E:/DetectorRunOutput/test LM_20190509093533/02BS12_AU_02a_files_All_AllSpecies_Model_Applied_probs.txt",stringsAsFactors = FALSE)
Selections<-read.delim("//nmfs/akc-nmml/CAEP/Acoustics/Projects/Dans Detectors/LF moan project/Results folders/LM AL16_AU)CL1_20190517170936/02AL16_AU_CL1_files_All_AllSpecies_Model_Applied_probs.txt",stringsAsFactors = FALSE)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

RT<-substrRight(as.character(Selections$File),17)
RT<-substr(RT,1,13)
if(grepl("-",substrRight(as.character(Selections$File),17)[1])){
  RT<-strptime(RT,format='%y%m%d-%H%M%S',tz="UTC") #some files have different naming convention 
}else{
  RT<-strptime(RT,format='%y%m%d_%H%M%S',tz="UTC")
}

Selections$RTfile<-RT
Selections$RTFb<-RT
Selections$RTFe<-RT+(Selections$End.Time..s.-Selections$Begin.Time..s.)

pointsDate<-as.POSIXlt((as.numeric(Selections$RTFb)+Selections$FileOffsetBegin), origin="1970-01-01",tz="UTC")

plot(Selections$RTFb,Selections$LM.prob)

#####

Selections<-Selections[which(Selections$LM.prob>0.70),]

sound_filesSelections<-unique(Selections$File)

sound_files<-dir("E:/Full_datasets/BS14_AU_04")
sound_filesfullpath<-paste("E:/Full_datasets/BS14_AU_04/",dir("E:/Full_datasets/BS14_AU_04"),sep="")

destFile<-paste("E:/HG_datasets/BS14_AU_04/BN_yesUnion",sep="")

dir.create(destFile)

sound_files_same<-which(sound_files %in% sound_filesSelections)

file.copy(sound_filesfullpath[sound_files_same], destFile)









#some data analysis on the first detector results: 

TPFPFN<-read.delim("E:/DetectorRunOutput/test LM_20190507102454/LM_BS13_AU_02a_files_38-122_TPFPFN_Tab_Ravenformat.txt")

Labeled<-read.csv("E:/DetectorRunFiles/LMProcessed_GT_data/test LM_20190507102454_processedGT.csv")

Labeled$duration<-Labeled$End.Time..s.-Labeled$Begin.Time..s.
TPFPFN$duration<-TPFPFN$End.Time..s.-TPFPFN$Begin.Time..s.


ggplot(Labeled, aes(duration, fill = as.factor(detectionType))) + geom_density(alpha = 0.2)

ggplot(TPFPFN, aes(duration, fill = as.factor(detectionType))) + geom_density(alpha = 0.2)

table<-read.delim("//nmfs/akc-nmml/CAEP/Acoustics/Projects/Dans Detectors/RavenBLEDscripts/Data/Selection tables/BN/AW12_AU_CL1Sum/Aw12_AU_CL1_files_1-384.txt")

hist(table$Low.Freq..Hz.)
hist(table$High.Freq..Hz.)
hist(table$End.Time..s.-table$Begin.Time..s.)



table<-read.delim("//nmfs/akc-nmml/CAEP/Acoustics/Projects/Dans Detectors/RavenBLEDscripts/Data/Selection tables/BN/AW12_AU_CL1Sum/Aw12_AU_CL1_files_1-384.txt")








############
#create R package with functions. Using hilary parker guide
#install.packages("devtools")
library("devtools")
devtools::install_github("klutometis/roxygen")
library(roxygen2)


#analysis for measurement of boings PRR

Boings<-read.csv("T:/Detector/DetectorRunFiles/BNProcessed_GT_data/BN GT test_20190522132348_processedGT.csv")
Boings<-Boings[,-84]

BoingsHG<-Boings

BoingsHG<-Boings[which(Boings$detectionType==1),]

BoingsHGnowarn1<-Boings[which(Boings$V85==0),]
BoingsHGnowarn1and2<-BoingsHGnowarn1[which(BoingsHGnowarn1$V86==0),]
BoingsHGnowarn2<-Boings[which(Boings$V86==0),]

#BoingsHG2<-BoingsHG[which(BoingsHG$V83>=1),]

#BoingsHG3<-BoingsHG[which(BoingsHG$V81<350),]

#super high high grade 
#BoingsHG4<-BoingsHG[which(BoingsHG$V81<100000000),]
#BoingsHG4<-BoingsHG[which(BoingsHG$V82<5),]

#less high grade 
BoingsHG4<-BoingsHG[which(BoingsHG$V75<100000000|BoingsHG$V77<100000000|BoingsHG$V79<100000000),]
BoingsHG4<-BoingsHG4[which(BoingsHG4$V76<5|BoingsHG4$V78<5|BoingsHG4$V80<5),]
BoingsHG4$CrudeAverage<-NA
BoingsHG4$CrudeSde<-NA

for(n in 1:nrow(BoingsHG4)){
  val1<-BoingsHG4$V75[n]
  val2<-BoingsHG4$V76[n]
  val3<-BoingsHG4$V77[n]
  val4<-BoingsHG4$V78[n]
  val5<-BoingsHG4$V79[n]
  val6<-BoingsHG4$V80[n]
  
  IntVec<-c(val1,val3,val5)
  sdeVec<-c(val2,val4,val6)
  
  IntVec<-IntVec[which(IntVec<100000000&sdeVec<5)]
  
  BoingsHG4$CrudeAverage[n]<-mean(IntVec)
  BoingsHG4$CrudeSde[n]<-mean(sdeVec)
}

BoingsHG4<-BoingsHG4[which(!is.na(BoingsHG4$CrudeAverage)),]

#compare relative counts from each mooring
#plot(prop.table(table(BoingsHG4$MooringID))
#plot(prop.table(table(Boings$MooringID))


hist(BoingsHG4$V75[which(BoingsHG4$V75<200)],breaks=250)
hist(BoingsHG4$V77[which(BoingsHG4$V77<200)],breaks=80)
hist(BoingsHG4$V79[which(BoingsHG4$V79<200)],breaks=80)

hist(BoingsHG4$V81[which(BoingsHG4$V81<200&BoingsHG4$MooringName=="AW12_AU_CL1")],col='green',breaks=80,ylab='Counts',xlab="Frequency (Hz)",main="Average interval value for high quality Boings")
hist(BoingsHG4$V81[which(BoingsHG4$V81<200&BoingsHG4$MooringName=="AW13_AU_CL1")],col='blue',breaks=80,add=TRUE)
hist(BoingsHG4$V81[which(BoingsHG4$V81<200&BoingsHG4$MooringName=="AW12_AU_KZ1")],col='red',breaks=80,add=TRUE)
hist(BoingsHG4$V81[which(BoingsHG4$V81<200&BoingsHG4$MooringName=="AW14_AU_PH1")],col='purple',breaks=80,add=TRUE)
hist(BoingsHG4$V81[which(BoingsHG4$V81<200&BoingsHG4$MooringName=="AL17_AU_CC2")],col='yellow',breaks=80,add=TRUE)

legend("topleft", c("dist_x1", "dist_x2"), fill=c("green", "red","blue","orange","yellow"))
#or so it is stacked instead of overriding 
library(ggplot2)
qplot(BoingsHG4$CrudeAverage, binwidth = 1, fill = factor(BoingsHG4$MooringName),ylab='Counts',xlab="Frequency (Hz)",main="Average 1st interval value for high quality Boings")


hist(BoingsHG4$V76[which(BoingsHG4$V76<200)],breaks=80)
hist(BoingsHG4$V78[which(BoingsHG4$V78<200)],breaks=80)
hist(BoingsHG4$V80[which(BoingsHG4$V80<200)],breaks=80)

hist(BoingsHG4$V82[which(BoingsHG4$V82<200)],breaks=80)



means<-c()
means<-c(means,mean(BoingsHG4$V75[which(BoingsHG4$V75<100000000)]))
means<-c(means,mean(BoingsHG4$V77[which(BoingsHG4$V77<100000000)]))
means<-c(means,mean(BoingsHG4$V79[which(BoingsHG4$V79<100000000)]))
means<-c(means,mean(BoingsHG4$V81[which(BoingsHG4$V81<100000000)]))


sd_errs<-c()
sd_errs<-c(sd_errs,mean(BoingsHG4$V76[which(BoingsHG4$V76<100000000)]))
sd_errs<-c(sd_errs,mean(BoingsHG4$V78[which(BoingsHG4$V78<100000000)]))
sd_errs<-c(sd_errs,mean(BoingsHG4$V80[which(BoingsHG4$V80<100000000)]))
sd_errs<-c(sd_errs,mean(BoingsHG4$V82[which(BoingsHG4$V82<100000000)]))

n<-c()
n<-c(n,length(BoingsHG4$V75[which(BoingsHG4$V75<100000000)]))
n<-c(n,length(BoingsHG4$V77[which(BoingsHG4$V77<100000000)]))
n<-c(n,length(BoingsHG4$V79[which(BoingsHG4$V79<100000000)]))
n<-c(n,length(BoingsHG4$V81[which(BoingsHG4$V81<100000000)]))

nsd<-c()
nsd<-c(nsd,length(BoingsHG4$V76[which(BoingsHG4$V76<100000000)]))
nsd<-c(nsd,length(BoingsHG4$V78[which(BoingsHG4$V78<100000000)]))
nsd<-c(nsd,length(BoingsHG4$V80[which(BoingsHG4$V80<100000000)]))
nsd<-c(nsd,length(BoingsHG4$V82[which(BoingsHG4$V82<100000000)]))

#look at a processed dataset to see if any variables are proxy for frequency (not freq normalized)

Labeled<-read.csv("T:/Detector/DetectorRunFiles/GSProcessed_GT_data/test new plot and stats GS_20190528175129_processedGT.csv")
min(Labeled$Low.Freq..Hz.)
hist(Labeled$Low.Freq..Hz.)
max(Labeled$High.Freq..Hz.)
hist(Labeled$High.Freq..Hz.)

for(f in 22:ncol(Labeled)){
  hist(Labeled[,f],main=colnames(Labeled)[f])
}
