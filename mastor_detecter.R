#Evaluate and score custom BLED detector(s) in R. Sox must be installed and pathed to accordingly

#catherine vandalism
###############################################


#RUN RAVEN DETECTORS


################################################

#install.packages("e1071") install.packages("Rtools",repos = "http://cran.r-project.org")install.packages("randomForest")install.packages("seewave")install.packages("tuneR")install.packages("plotrix")install.packages("aod")install.packages("ggplot2", dep = TRUE)install.packages("usdm")install.packages("ROCR")install.packages("e1071") install.packages("caret")install.packages("ModelMetrics")install.packages("stringi")install.packages("signal")install.packages("beepr")install.packages("Rraven")install.packages("flightcallr", repos="http://R-Forge.R-project.org")install.packages("plotrix") install.packages("oce") install.packages("imager") install.packages("obliqueRF") install.packages("fpc") install.packages("devtools") devtools::install_github("Azure/rAzureBatch") devtools::install_github("Azure/doAzureParallel")




library(e1071)   
library(foreach)
library(doParallel)
#library(caret)  
library(randomForest)
library(seewave)
library(tuneR)
library(plotrix)
library(aod)
library(ggplot2)
library(usdm)
library(flightcallr)
library(ROCR)
library(ModelMetrics)
library(Rraven)
library(seewave)
library(beepr)
library(stringr)
library(stringi)
library(signal)
library(oce)
library(imager)
library(Cairo)
library(obliqueRF)
library(fpc)
library(RSAGA)
#library(rAzureBatch)
#library(doAzureParallel)

###############################
#run this stuff so I can know whether to start azure parallel. Run it again later to populate cluster w variable. 
#dumb conditional so I don't have to change path from machine to machine
if(dir.exists("C:/Users/ACS-3")){
  user<-"ACS-3"
  drivepath<-"F:/"
  gitPath<-paste(drivepath,"RavenBLEDscripts/",sep="")
  
}else if(dir.exists("C:/Users/danby456")){
  user<-"danby456"
  drivepath<-"E:/"
  gitPath<-paste(drivepath,"RavenBLEDscripts/",sep="")
  
}else if(dir.exists("C:/Users/Daniel.Woodrich")){
  user<-"Daniel.Woodrich"
  drivepath<-"E:/"
  gitPath<-"//nmfs/akc-nmml/CAEP/Acoustics/Projects/Dans Detectors/RavenBLEDscripts/"
}

startcombpath<-paste(drivepath,"Combined_sound_files/",sep="")
BLEDpath<-paste("C:/Users/",user,"/Raven Pro 1.5/Presets/Detector/Band Limited Energy Detector/",sep="")
ravenpath<-paste("C:/Users/",user,"/Raven Pro 1.5",sep="")
outputpath<-paste(drivepath,"DetectorRunOutput/",sep="")
outputpathfiles<-paste(drivepath,"DetectorRunFiles/",sep="")

###############input parameters
#detector control
ControlTab<-read.csv(paste(gitPath,"Data/CallParams/Detector_control.csv",sep=""))
ControlTab[,3]<-as.character(ControlTab[,3])
#####################################

parallelType<- ControlTab[which(ControlTab[,2]=="parallelType"),3]

#if(parallelType=="azure"){
#  
#setwd(getwd()) #how the hell this line does anything I don't know. But it fixed a pathing error.
#setCredentials("credentials.json")
#
#cluz <- doAzureParallel::makeCluster("cluster.json")
#registerDoAzureParallel(cluz) 
#
#setAutoDeleteJob(FALSE)
#setHttpTraffic(TRUE)
#
#}

startLocalPar<-function(...){
  if(user=="ACS-3"){
    num_cores <- detectCores()
  }else{
    num_cores <- detectCores()-1
  }
  cluz <<- parallel::makeCluster(num_cores)
  registerDoParallel(cluz)
  
  clusterExport(cluz, c(...))
  
}
#generateCredentialsConfig("credentials.json")
#generateClusterConfig("cluster.json")
#setCredentials("credentials.json")

# Create your cluster if it does not exist; this takes a few minutes
#cluster <- doAzureParallel::makeCluster("cluster.json") 

#  parallel::stopCluster(cluster)
# Register your parallel backend 
#registerDoAzureParallel(cluster) 

# Check that the nodes are running 
#doAzureParallel::getDoParWorkers() 

formatModelData<-function(dataPostModel){
  #pickup here or leave it for now. 
  newSpec<-NULL
  for(s in 1:length(spec)){
    newSpec<-c(newSpec,(nospec+3*s-2))
  }
  
  dataPostModel<-dataPostModel[,c(1:nospec,newSpec)]  
  
  #put columns excluded from model back in ?
  
  #change data columns from factor to numeric that should be numeric:
  for(b in c(1,2,3,4,5,10,11)){
    dataPostModel[,b]<-as.numeric(as.character(dataPostModel[,b])) 
    
  }
  
  for(b in c(13,14)){
    dataPostModel[,b]<-as.integer(as.numeric(as.POSIXlt(dataPostModel[,b])))
  }
  
  Tab<-NULL
  Tab2<-NULL
  
  dataPostModellabs<-factorLevels(dataPostModel)
  
  for(s in 1:length(spec)){
    
    loadSpecVars(spec[s])
    
    dataPostModelspec<-dataPostModel[grep(spec[s],dataPostModel$Species),]
    if(s!=1){
      if(any(grepl(spec[s],dataPostModel$Species))){
      Tab<-Tab[-grep(spec[s],Tab$Species),]
      }
    }
    
    #adaptively combine detections based on probability
    dataPostModelspecLabs<-factorLevels(dataPostModelspec)
    
    dataPostModelspec$Unq_Id<-as.numeric(as.factor(paste(dataPostModelspec$Species,dataPostModelspec$MooringID)))
    
    for(i in 1:ncol(dataPostModelspec)){
      if(typeof(dataPostModelspec[,i])=="character"){
        dataPostModelspec[,i]<-as.factor(dataPostModelspec[,i])
      }
      
    }
    
    dataPostModelspecMat<- data.matrix(dataPostModelspec)
    
    #take out temporarily
    #dataPostModelspecMat<-adaptive_compare(dataPostModelspecMat) 
    
    #simulate context over time using probability scores.Make new tab with selection ID and all supporting variables- don't pin on data3mat. Can use in place of prob if wanted. 
    CS_output<-context_sim(dataPostModelspecMat)
    
    #remove missing rows:
    dataPostModelspecMat<-dataPostModelspecMat[order(dataPostModelspecMat[,1]),]
    dataPostModelspec<-data.frame(dataPostModelspecMat)
    dataPostModelspec$Unq_Id<-NULL
    dataPostModelspec<-applyLevels(dataPostModelspec,dataPostModelspecLabs)
    
    Tab<-rbind(Tab,dataPostModelspec)
    Tab2<-rbind(Tab2,CS_output)
  }
  
  dataPostModel<-Tab
  CS_output<-Tab2
  
  dataPostModel<-applyLevels(dataPostModel,dataPostModellabs)

  return(list(dataPostModel,CS_output))
}

factorLevels<-function(dataToMat){
  FactorTab<-lapply(dataToMat,as.numeric)
  datType<-lapply(dataToMat,typeof)
  datLab<-lapply(dataToMat,as.character)
  return(list(FactorTab,datType,datLab))
}

applyLevels<-function(matToData,datLevels){

  keep<-which(datLevels[[1]][[1]] %in% matToData[,1])

  for(i in 1:length(datLevels[[1]])){
    if(datLevels[[2]][[i]]=="integer"){
      matToData[,i]<-datLevels[[3]][[i]][keep]
    }
    
  }
  return(matToData)
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


varImpPlot_AVG <- function(x, sort=TRUE,
                       n.var=min(30, nrow(x)),
                       type=NULL, class=NULL, scale=TRUE, 
                       main=deparse(substitute(x)), ...) {
  imp <- x
  ## If there are more than two columns, just use the last two columns.
    ord <- if (sort) rev(order(imp[,1],
                               decreasing=TRUE)[1:n.var]) else 1:n.var
    xmin <- if (colnames(imp)[1] %in%
                c("IncNodePurity", "MeanDecreaseGini")) 0 else min(imp[ord, 1])
    dotchart(imp[ord,1], xlab=colnames(imp)[1], labels=rownames(giniAv)[ord],
             main=main,xlim=c(xmin, max(imp[,1])), ...)

    invisible(imp)
}

rotate <- function(x) t(apply(x, 2, rev))

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

lastFeature<-function(a,b){
  
  ## get frequency windows length for smoothing
  step <- a/128/1000
  fsmooth = 0.1
  fsmooth <- fsmooth/step
  
  ## number of samples
  n <- nrow(b)
  
  ## smoothing parameter
  FWL <- fsmooth - 1
  
  ## smooth 
  zx <- apply(as.matrix(1:(n - FWL)), 1, function(y) sum(b[y:(y + FWL), 2]))
  zf <- seq(min(b[,1]), max(b[,1]), length.out = length(zx))
  
  ## make minimum amplitude 0
  zx <- zx - min(zx)
  zx[zx < 0] <- 0
  
  ## normalize amplitude from 0 to 1
  zx <- zx/max(zx)
  
  # return low and high freq
  return(zf[which.max(zx)] + (step / 2))
}



#decimate dataset. 
decimateData<-function(m,oldPath,sfpath){
  print(paste("decimate folder",oldPath))
  
  fstart<-MoorInfo[m,2]
  fend<-MoorInfo[m,3]
    for(z in dir(oldPath,pattern=".wav")[fstart:fend]){
      print(paste(z,"to",sfpath))
      wav<-readWave(paste(oldPath,"/",z,sep=""),unit="sample")
      wav.samp.rate<-wav@samp.rate
      if(wav@samp.rate<16384){
        dfact<-16384/wav@samp.rate
        decdo<-prime.factor(decimationFactor/dfact)
        }else{
          decdo<-decimationSteps
        }
      
        wav<-wav@left
        
        for(h in decdo){
          wav<-signal::decimate(wav,h)
        }
        
        wav <- Wave(wav, right = numeric(0), samp.rate = 16384/decimationFactor)

        writeWave.nowarn(wav, filename=paste(sfpath,"/",z,sep=""),extensible = FALSE)
      }
      write.table(paste("This data has been decimated by factor of",decdo),paste(sfpath,"/decimationStatus.txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE,col.names=FALSE)
      
}


inputGT<-function(){
  if(useMasterGT=="y"){
    if(addToMaster=="n"){
      GTset<<-read.csv(paste(gitPath,"Data/GroundTruth.csv",sep=""))
      TPtottab<<-read.csv(paste(gitPath,"Data/TotalTP_GT.csv",sep="")) #produces most recently modifed file 
      MoorInfo<<-read.csv(paste(gitPath,"Data/MoorInfo.csv",sep=""))
    }else if(addToMaster=="y"){
      GTset<<-rbind(GTset,read.csv(paste(gitPath,"Data/GroundTruth.csv",sep="")))
      TPtottab<<-rbind(TPtottab,read.csv(paste(gitPath,"Data/TotalTP_GT.csv",sep="")))
      MoorInfo<<-rbind(MoorInfo,read.csv(paste(gitPath,"Data/MoorInfo.csv",sep="")))
    }
    
  }else if(useMasterGT=="n"){
    
    
    if(length(spec)==1){
      recentTab<<-file.info(list.files(paste(outputpathfiles,spec,"Processed_GT_data/",sep=""), full.names = T))
      recentPath<<-rownames(recentTab)[which.max(recentTab$mtime)]
      GTset<<-read.csv(recentPath) #produces most recently modifed file 
      
      recentTab<<-file.info(list.files(paste(outputpathfiles,spec,"TPtottab/",sep=""), full.names = T))
      recentPath<<-rownames(recentTab)[which.max(recentTab$mtime)]
      TPtottab<<-read.csv(recentPath) #produces most recently modifed file 
      
    }else if(length(spec)>1){
      recentTab<<-file.info(list.files(paste(outputpathfiles,"Processed_GT_data/",sep=""), full.names = T))
      recentPath<<-rownames(recentTab)[which.max(recentTab$mtime)]
      GTset<<-read.csv(recentPath) #produces most recently modifed file 
      
      recentTab<<-file.info(list.files(paste(outputpathfiles,"TPtottab/",sep=""), full.names = T))
      recentPath<<-rownames(recentTab)[which.max(recentTab$mtime)]
      TPtottab<<-read.csv(recentPath) #produces most recently modifed file 
    }
    
  }
  
  GTset<<-GTset[,c(1,4:length(GTset))]
  GTset$detectionType<<-as.factor(GTset$detectionType)
  
}


loadSpecVars<-function(whichSpec){
  #save species variables in global environment
  
  #species specific
  ParamsTab<-read.csv(paste(gitPath,"Data/CallParams/",whichSpec,".csv",sep=""))
  ParamsTab[,3]<-as.character(ParamsTab[,3])
  
  #Species specific 
  #Signal manipulation
  if(whichRun=="GT"){
  Decimate<<-ParamsTab[which(ParamsTab[,2]=="Decimate"),3] 
  decimationFactor<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="decimationFactor"),3] )
  whiten<<-ParamsTab[which(ParamsTab[,2]=="whiten"),3]
  FO<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="FO"),3] )
  LMS<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="LMS"),3])
  Filtype<<-ParamsTab[which(ParamsTab[,2]=="Filtype"),3] 
  
  }else if(whichRun=="NEW"){
    Decimate<<-ControlTab[which(ControlTab[,2]=="Decimate"),3] 
    decimationFactor<<-as.numeric(ControlTab[which(ControlTab[,2]=="decimationFactor"),3] )
    whiten<<-ControlTab[which(ControlTab[,2]=="whiten"),3]
    FO<<-as.numeric(ControlTab[which(ControlTab[,2]=="FO"),3] )
    LMS<<-as.numeric(ControlTab[which(ControlTab[,2]=="LMS"),3])
    Filtype<<-ControlTab[which(ControlTab[,2]=="Filtype"),3] 

  }
  #Raven Detectors
  spStart<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="spStart"),3])
  spEnd<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="spEnd"),3])
  
  #Context sim
  greatcallThresh<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="greatcallThresh"),3])
  maxBonus<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="maxBonus"),3])
  goodcallBonus<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="goodcallBonus"),3]) 
  maxPenalty<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="maxPenalty"),3] )
  badcallPenalty<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="badcallPenalty"),3] )
  
  #Detection Processing Spread (algo)
  grpsize<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="grpsize"),3] )
  allowedZeros<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="allowedZeros"),3]) 
  detskip<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="detskip"),3] )
  groupInt<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="groupInt"),3] )
  Maxdur<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="Maxdur"),3])
  Mindur<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="Mindur"),3])
  
  #RW algo
  if(whichSpec=="RW"){
    downsweepCompMod<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="downsweepCompMod"),3])
    downsweepCompAdjust<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="downsweepCompAdjust"),3])
    timesepGS<<-NULL #placeholder so that I can export variable for parallel 
    timesepBP<<-NULL
  }
  
  #GS algo
  if(whichSpec=="GS"|whichSpec=="BP"){
    if(whichSpec=="GS"){
    timesepGS<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="timesepGS"),3] )
    timesepBP<<-NULL
    }else if(whichSpec=="BP"){
    timesepBP<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="timesepBP"),3] )
    timesepGS<<-NULL #placeholder so that I can export variable for parallel 
    }
    downsweepCompMod<<-NULL #placeholder so that I can export variable for parallel 
    downsweepCompAdjust<<-NULL #placeholder so that I can export variable for parallel 
    
  }
  
  #Image analysis
  ImgThresh<<-paste(ParamsTab[which(ParamsTab[,2]=="ImgThresh"),3],"%",sep="")
  
  #adaptive Compare
  timediffself<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="timediffself"),3])
  probdist<<-as.numeric(ParamsTab[which(ParamsTab[,2]=="probdist"),3])
  
  #GT
  if(runGT=="y"){
    GTmoorings<<- str_split(ParamsTab[which(ParamsTab[,2]=="GTmoorings"),3],",",simplify=TRUE)
    GTsf<<-str_split(ParamsTab[which(ParamsTab[,2]=="GTsf"),3],",",simplify=TRUE)
    GTpath<<-str_split(ParamsTab[which(ParamsTab[,2]=="GTpath"),3],",",simplify=TRUE)
    GTpath2<<-str_split(ParamsTab[which(ParamsTab[,2]=="GTpath2"),3],",",simplify=TRUE)
    GTsourceFormat<<-str_split(ParamsTab[which(ParamsTab[,2]=="GTsourceFormat"),3],",",simplify=TRUE)
  }
  
  
  write.csv(ParamsTab,paste(outputpath,runname,"/",whichSpec,".csv",sep=""),row.names = FALSE)
  
  ####################some var related calculations
  detectorssprshort<<- list.files(BLEDpath)[spStart:spEnd] 
  
  detlist<<-NULL
  detlist2<<-substr(detectorssprshort[1],1,3)
  
  detnum<<-1
  
  if(whiten=="n"){
    FO<<-NA
    LMS<<-NA
  }
  
  detectorsspr<<-NULL
  for(j in 1:length(detectorssprshort)){
    detectorsspr[j]<<-paste(BLEDpath,detectorssprshort[j],sep="")
  }
  
  if(whichSpec=="RW"){
    ravenView<<-"RW_Upcalls"
  }else if(whichSpec=="GS"|whichSpec=="BP"){
    ravenView<<-"RW_GS"
  }
  
  if(Decimate=="y"){
    ravenView<<-paste(ravenView,"_",decimationFactor,"Decimate",sep="")
  }
  
  decimationSteps<<-prime.factor(decimationFactor)
  
  
  ###############################
}

makeMoorInfo<-function(moorings,sf,path,pathspec,sourceFormat,curSpec){

  status<-sf=="full"
  type<-rep("partial",length(sf))
  type[which(status)]<-"all"
  
  if(any(sf=="full")){
    for(n in which(sf=="full")){
      if(path[n]=="HG_datasets"){
        lookup<-paste(drivepath,path[n],moorings[n],paste(pathspec[n],"_yesUnion",sep=""),sep="/")
      }else if(path[n]=="Full_datasets"){
        lookup<-paste(drivepath,path[n],moorings[n],sep="/")
      }
      if(sourceFormat[n]=="month"){
        #save this section for when I am looking at a folder with months in it. 
      }else if(sourceFormat[n]=="open"){
        fileEnd<-length(list.files(lookup,pattern=".wav"))
        sf[n]<-paste(1,fileEnd,sep=":")
      }
    }
  }
  
  #test for file name formatting and convert to numeric
  if(any(!grepl("[^A-Za-z]", sf))){
    for(n in which(!grepl("[^A-Za-z]", sf))){
      sf1<-c(as.character(str_split(sf[n],":",simplify=TRUE)[1]))
      sf2<-c(as.character(str_split(sf[n],":",simplify=TRUE)[2]))
      if(path[n]=="HG_datasets"){
        lookup<-paste(drivepath,path[n],moorings[n],paste(pathspec[n],"_yesUnion",sep=""),sep="/")
      }else if(path[n]=="Full_datasets"){
        lookup<-paste(drivepath,path[n],moorings[n],sep="/")
      }
      if(sourceFormat[n]=="month"){
        #save this section for when I am looking at a folder with months in it.  
      }else if(sourceFormat[n]=="open"){
        fileStart<-which(list.files(lookup,pattern=".wav")==sf1)
        fileEnd<-which(list.files(lookup,pattern=".wav")==sf2)
        sf[n]<-paste(fileStart,fileEnd,sep=":")
      }
    }
  }
  
  sf1<-list()
  sf2<-list()
  for(n in 1:length(sf)){
    sf1[[n]]<-c(as.numeric(str_split(sf[n],":",simplify=TRUE)[1]))
    sf2[[n]]<-c(as.numeric(str_split(sf[n],":",simplify=TRUE)[2]))
  }
  
  #add file name format
  sf3<-list()
  sf4<-list()
  for(n in 1:length(sf)){
    if(path[n]=="HG_datasets"){
      lookup<-paste(drivepath,path[n],moorings[n],paste(pathspec[n],"_yesUnion",sep=""),sep="/")
    }else if(path[n]=="Full_datasets"){
      lookup<-paste(drivepath,path[n],moorings[n],sep="/")
    }
    if(sourceFormat[n]=="month"){
      #save this section for when I am looking at a folder with months in it. 
    }else if(sourceFormat[n]=="open"){
      sf3[[n]]<-list.files(lookup,pattern=".wav")[sf1[[n]]]
      sf4[[n]]<-list.files(lookup,pattern=".wav")[sf2[[n]]]
    }
  }
  
  MoorsUniqueIDS<-NULL
  for(n in 1:length(sf)){
    if(type[n]!="all"){
      MoorsUniqueIDS[n]<-paste(moorings[n],"_files_",sf1[n],"-",sf2[n],sep="")
    }else if(type[n]=="all")
      MoorsUniqueIDS[n]<-paste(moorings[n],"_files_All",sep="")
  }
  
  MoorsInfo<-cbind(t(moorings),as.numeric(sf1),as.numeric(sf2),as.character(sf3),as.character(sf4),t(t(as.numeric(sf2)-as.numeric(sf1))),t(path),t(sourceFormat),t(t(rep(curSpec,length(sf)))),t(t(MoorsUniqueIDS)),t(pathspec))
  
  return(MoorsInfo)
}


parAlgo<-function(dataaa){
  
  if(parallelType=="local"){
  startLocalPar("Matdata","detskip","downsweepCompMod","downsweepCompAdjust","allowedZeros","grpsize","RW_algo","GS_algo","BP_algo","timesepGS","timesepBP","s")
  }
  print("start algo")
  if(s=="RW"){
    wantedSelections<-foreach(grouppp=unique(dataaa[,2])) %dopar% {
      group<-dataaa[which(dataaa[,2]==grouppp),1:3]
      RW_algo(resltsTabmat=group)
    }
    wantedSelections<-as.integer(do.call('c', wantedSelections))
  }else if(s=="GS"){
    wantedSelections<-foreach(grouppp=unique(dataaa[,2])) %dopar% {
      group<-dataaa[which(dataaa[,2]==grouppp),c(1,2,4,5)]
      GS_algo(groupdat=group)
    }
    wantedSelections<-do.call('cbind', wantedSelections)
  }else if(s=="BP"){
    wantedSelections<-foreach(grouppp=unique(dataaa[,2])) %dopar% {
      group<-dataaa[which(dataaa[,2]==grouppp),c(1,2,4,5)]
      BP_algo(groupdat=group)
    }
    wantedSelections<-do.call('cbind', wantedSelections)
  }
  print("end algo")
  if(parallelType=="local"){
    parallel::stopCluster(cluz)
  }
  return(wantedSelections)
  
}

dataArrangeModel<-function(dataForModel){
  dataForModel$year<-format(as.Date(dataForModel$RTfile),"%y")
  dataForModel$month<-format(as.Date(dataForModel$RTfile),"%m")
  
  dataForModel$detectionType<-as.numeric(as.character(dataForModel$detectionType),sep="")
  
  dataForModel$detectionType<-as.factor(paste(dataForModel$Species,dataForModel$detectionType))
  
  modelDat<-dataForModel[,c(1,18:19,21:(ncol(dataForModel)-2))]
  modelDatFactors<-dataForModel[,c(8,17,ncol(dataForModel)-1,ncol(dataForModel))]
  modelDatFactors<-apply(modelDatFactors,2,function(x) as.factor(x))
  
  otherDat<-dataForModel[,c(1:7,9:16)]
  
  modelDat<-apply(modelDat,2,function(x) as.numeric(as.character(x)))
  removeMat<-apply(modelDat,2,function(x) !is.finite(x))
  removeVec<-which(apply(removeMat,1,function(x) any(x)))
  
  if(length(removeVec)>0){
  modelDat<-modelDat[-removeVec,]
  modelDatFactors<-modelDatFactors[-removeVec,]
  otherDat<-otherDat[-removeVec,]
  }
  
  modelDat<-apply(modelDat,2,function(x) na.roughfix(x))
  
  modelDat<-cbind(data.frame(modelDat),data.frame(modelDatFactors))
  
  return(list(modelDat,otherDat))
}

RW_algo<-function(resltsTabmat){
  wantedSelections<-c()
  groupdat<- resltsTabmat
  grpvec<-groupdat[,1]
  colClasses = c("numeric","numeric","numeric","numeric","numeric")
  runsum<- read.csv(text="start, ones, zeros, length,skip", colClasses = colClasses)
  groupdat<-cbind(groupdat,matrix(99,nrow(groupdat),nrow(groupdat)-(grpsize-1)))
  
  for(g in 1:(nrow(groupdat)-(grpsize-1))){
    RM<-groupdat[g,1]
    groupdat[g,3+g]<-2
    skipvec<-0
    for(h in g:(nrow(groupdat)-1)){
      rsltvec0s<-rle(groupdat[,3+g][! groupdat[,3+g] %in% 98])
      if(any(rsltvec0s$lengths[rsltvec0s$values==0]>allowedZeros)){
        for(q in 1:allowedZeros){
          groupdat[h-(q-1),4+g]<-99
        }
        break
      }  
      if(RM<grpvec[h+1]&RM+(detskip+1)>grpvec[h+1]&groupdat[h,3]!=groupdat[h+1,3]){
        groupdat[h+1,3+g]<-1
        skipvec<-c(skipvec,(grpvec[h+1]-RM))
        RM<-grpvec[h+1]
      }else if(groupdat[h,3]!=groupdat[h+1,3]&RM+(detskip+1)>grpvec[h+1]&RM-(detskip+1)<=grpvec[h+1]){
        groupdat[h+1,3+g]<-0
      }else if(groupdat[h,3]!=groupdat[h+1,3]&(RM+(detskip+1)<=grpvec[h+1]|RM-(detskip+1)>grpvec[h+1])){
        groupdat[h+1,3+g]<-98
      }
      if(groupdat[h,3]==groupdat[h+1,3]&groupdat[h,3+g]==0&RM<grpvec[h+1]&RM+(detskip+1)>grpvec[h+1]){
        groupdat[h+1,3+g]<-1
        groupdat[h,3+g]<-98
        skipvec<-c(skipvec,(grpvec[h+1]-RM))
        RM<-grpvec[h+1]
      }else if(groupdat[h,3]==groupdat[h+1,3]){
        groupdat[h+1,3+g]<-98
      }
    }
    #replace 1s that are more than 0.5s apart from last in sequence to 0s
      if(any(groupdat[,4]==1)){
        grouptemp<-groupdat[which(groupdat[,4]==2|groupdat[,4]==1),]
        for(c in 1:(nrow(grouptemp)-1)){
          if(grouptemp[c,3]+0.5<grouptemp[c+1,3]){
            grouptemp[(c+1):nrow(grouptemp),4]<-99
            groupdat[which(rownames(groupdat) %in% rownames(grouptemp)),3+g]<-grouptemp[,4]
            break
          }
        }
      grouptemp<-NULL
        }
      
    runsum[g,1]<-g
    runsum[g,2]<-sum(groupdat[,3+g]==1)
    runsum[g,3]<-sum(groupdat[,3+g]==0)
    runsum[g,4]<-(sum(groupdat[,3+g]==1)+sum(groupdat[,3+g]==0)+1)
    runsum[g,5]<-max(skipvec)
  }
  runsum<-runsum[which(runsum[,2]==max(runsum[,2])),] #choose w most ones
  runsum<-runsum[which(runsum[,5]==min(runsum[,5])),] #choose w smallest maximum skip (most gradual)
  runsum<-runsum[which(runsum[,3]==min(runsum[,3])),] #choose w least 0s
  runsum<-runsum[which(runsum[,4]==min(runsum[,4])),] #choose w least length
  runsum<-runsum[1,] #choose first one
  
  #find other good runs g
  
  #if run is less than 33% of boxes, build a downsweep. If the downsweep has equal or more ones disqualify it.
  kill="n"
  if((runsum[,2]+1)<grpsize){
    kill="y"
  }
  HowBoutNot<-TRUE
  if(((runsum[,2]+1)*downsweepCompMod)<nrow(groupdat)&kill=="n"&HowBoutNot!=TRUE){
    groupdat2<- resltsTabmat
    groupdat2<-groupdat2[order(groupdat2[,3],-groupdat2[,1]),]#reverse the order it counts stacks detections
    grpvec2<-groupdat2[,1]
    colClasses = c("numeric","numeric","numeric","numeric","numeric")
    runsum2<- read.csv(text="start, ones, zeros, length,skip", colClasses = colClasses)
    groupdat2<-cbind(groupdat,matrix(99,nrow(groupdat2),nrow(groupdat2)-(grpsize-1)))
    
    for(g in 1:(nrow(groupdat2)-(grpsize-1))){
      RM2<-groupdat2[g,1]
      groupdat2[g,3+g]<-2
      skipvec2<-0
      for(h in g:(nrow(groupdat2)-1)){
        rsltvec0s2<-rle(groupdat2[,3+g][! groupdat2[,3+g] %in% 98])
        if(any(rsltvec0s2$lengths[rsltvec0s2$values==0]>allowedZeros)){
          for(q in 1:allowedZeros){
            groupdat2[h-(q-1),4+g]<-99
          }
          break
        }  
        if(RM2>grpvec2[h+1]&RM2-(detskip+1)<grpvec2[h+1]&groupdat2[h,3]!=groupdat2[h+1,3]){
          groupdat2[h+1,3+g]<-1
          skipvec2<-c(skipvec2,(grpvec2[h+1]-RM2))
          RM2<-grpvec2[h+1]
        }else if(groupdat2[h,3]!=groupdat2[h+1,3]&RM2-(detskip+1)<grpvec2[h+1]&RM2+(detskip+1)>=grpvec2[h+1]){
          groupdat2[h+1,3+g]<-0
        }else if(groupdat2[h,3]!=groupdat2[h+1,3]&(RM2-(detskip+1)>=grpvec2[h+1])|RM2+(detskip+1)<grpvec2[h+1]){
          groupdat2[h+1,3+g]<-98
        }
        if(groupdat2[h,3]==groupdat2[h+1,3]&(groupdat2[h,3+g]==0|groupdat2[h,3+g]==98)&RM2>grpvec2[h+1]&(RM2-(detskip+1))<grpvec2[h+1]){
          groupdat2[h+1,3+g]<-1
          groupdat2[h,3+g]<-98
          skipvec2<-c(skipvec2,(grpvec2[h+1]-RM2))
          RM2<-grpvec2[h+1]
        }else if(groupdat2[h,3]==groupdat2[h+1,3]){
          groupdat2[h+1,3+g]<-98
        }
      }
      runsum2[g,1]<-g
      runsum2[g,2]<-sum(groupdat2[,3+g]==1)
      runsum2[g,3]<-sum(groupdat2[,3+g]==0)
      runsum2[g,4]<-(sum(groupdat2[,3+g]==1)+sum(groupdat2[,3+g]==0)+1)
      runsum2[g,5]<-max(-skipvec2)
    }
    runsum2<-runsum2[which(runsum2[,2]==max(runsum2[,2])),] #choose w most ones
    runsum2<-runsum2[which(runsum2[,3]==min(runsum2[,3])),] #choose w least 0s
    runsum2<-runsum2[which(runsum2[,5]==min(runsum2[,5])),] #choose w smallest maximum skip (most gradual)
    runsum2<-runsum2[which(runsum2[,4]==min(runsum2[,4])),] #choose w least length
    runsum2<-runsum2[1,] #choose first one
    
    downsweep<-groupdat2[,c(1:3,3+runsum2[,1])]
    downsweep<-downsweep[which(downsweep[,4]==2|downsweep[,4]==1),]
    if(is.null(nrow(downsweep))){
      kill="n"
    }else{
      downsweep<-downsweep[,c(1:3)]
      upsweep<-groupdat[,c(1:3,3+runsum[,1])]
      upsweep<-upsweep[which(upsweep[,4]==2|upsweep[,4]==1),]
      upsweep<-upsweep[,c(1:3)]
      if(mean(downsweep[,3])<mean(upsweep[,3])&downsweep[1,1]<upsweep[nrow(upsweep),1]&downsweep[nrow(downsweep),3]==upsweep[1,3]){
        kill="s"
        wantedSelections<-c(rownames(downsweep[,]))
        wantedSelections<-c(rownames(upsweep[2:nrow(upsweep),]))
      }else if(runsum2[,2]>=runsum[,2]+downsweepCompAdjust){
        kill="y"
      }else{
        kill="n"
      }
    }
  }
  
  
  if(kill=="n"){
    groupdat<-groupdat[,c(1:3,3+runsum[,1])]
    wantedSelections<-c(rownames(groupdat[which(groupdat[,4]==2|groupdat[,4]==1),]))
  }
  return(wantedSelections)
}

GS_algo<-function(groupdat){
  #for(f in unique(resltsTabmat[,2])){
  groupdat<-groupdat[order(-groupdat[,1],groupdat[,3]),]#reverse the order it counts stacks detection
  groupdat<-cbind(groupdat,matrix(99,nrow(groupdat),nrow(groupdat)-(grpsize-1)))
  
  for(g in 1:(nrow(groupdat)-(grpsize-1))){
    groupdat[g,4+g]<-2
    RT<-groupdat[g,3]
    RM<-groupdat[g,1]
    skipvec<-0
    if(any(groupdat[g,5:ncol(groupdat)]==3)&g>1){
      #do not compute run
    }else if(g>=1){
    for(h in g:(nrow(groupdat)-1)){
      if(RM>groupdat[h+1,1]&(((RT-groupdat[h+1,3]-(timesepGS/((groupdat[h+1,1]+1.25)^1.1)))<0&(RT-groupdat[h+1,3]+(timesepGS/((groupdat[h+1,1]+1.25)^1.1)))>0)|((RT-groupdat[h+1,4]-(timesepGS/((groupdat[h+1,1]+1.25)^1.1)))<0&(RT-groupdat[h+1,4]+(timesepGS/((groupdat[h+1,1]+1.25)^1.1)))>0))&(RM-groupdat[h+1,1])<(detskip+1)){
        if(((RT-groupdat[h+1,3]-(timesepGS/((groupdat[h+1,1]+1.25)^1.1)))<0&(groupdat[h+1,3]-RT+(timesepGS/((groupdat[h+1,1]+1.25)^1.1)))>0)&((timesepGS/((groupdat[h+1,1]+1.25)^1.1))<0.2|groupdat[h+1,3]>groupdat[h,3])){
          boxPos<-3
        }else{
          boxPos<-4
        }
        groupdat[h+1,4+g]<-boxPos
        RT<-groupdat[h+1,boxPos]
        skipvec<-c(skipvec,(RM-groupdat[h+1,1]))
        RM<-groupdat[h+1,1]
      }else if(RM==groupdat[h+1,1]){
        groupdat[h+1,4+g]<-98
      }else{
        groupdat[h+1,4+g]<-98
      }
      if(g>1){
      if(any(groupdat[h+1,c(5:(3+g))]==groupdat[h+1,4+g])&(groupdat[h+1,4+g]!=98&groupdat[h+1,4+g]!=99&groupdat[h+1,4+g]!=0)){
        groupdat[h+1,4+g]<-0
        break
      }
      }
    }
    }
  }
  p=1
  rownombres=NULL
  rowID=NULL
  for(k in 5:ncol(groupdat)){
    if(length(which(groupdat[,k]==3|groupdat[,k]==4))>=grpsize-1){
      rownombres<-c(rownombres,as.integer(names(groupdat[which(groupdat[,k]==3|groupdat[,k]==4|groupdat[,k]==2),k])),999999999)
      rowID<-c(rowID,rep(p,each=length(names(groupdat[which(groupdat[,k]==3|groupdat[,k]==4|groupdat[,k]==2),k]))),999999999)
      p=p+1
    }
  }
  if(!is.null(rowID)){
 return(as.matrix(rbind(rownombres,rowID)))
  }
# }
}

BP_algo<-function(groupdat){
  
  #for(f in unique(dataaa[,2])){
  #groupdat<-dataaa[which(dataaa[,2]==f),c(1,2,4,5)]
    
  groupdat<-groupdat[order(-groupdat[,1],groupdat[,3]),]#reverse the order it counts stacks detection
  groupdat<-cbind(groupdat,matrix(99,nrow(groupdat),nrow(groupdat)-(grpsize-1)))
  
  for(g in 1:(nrow(groupdat)-(grpsize-1))){
    groupdat[g,4+g]<-2
    RT<-groupdat[g,3]
    RM<-groupdat[g,1]
    skipvec<-0
    if(any(groupdat[g,5:ncol(groupdat)]==3)&g>1){
      #do not compute run
    }else if(g>=1){
      for(h in g:(nrow(groupdat)-1)){
        if(RM>groupdat[h+1,1]&((RT-groupdat[h+1,3]-timesepBP)<0&(RT-groupdat[h+1,3]+timesepBP)>0)&(RM-groupdat[h+1,1])<(detskip+1)){
          groupdat[h+1,4+g]<-3
          RT<-groupdat[h+1,3]
          skipvec<-c(skipvec,(RM-groupdat[h+1,1]))
          RM<-groupdat[h+1,1]
        }else{
          groupdat[h+1,4+g]<-98
        }
        if(g>1){
          if(any(groupdat[h+1,c(5:(3+g))]==groupdat[h+1,4+g])&(groupdat[h+1,4+g]!=98&groupdat[h+1,4+g]!=99&groupdat[h+1,4+g]!=0)){
            groupdat[h+1,4+g]<-0
            break
          }
        }
      }
    }
  }
  p=1
  rownombres=NULL
  rowID=NULL
  for(k in 5:ncol(groupdat)){
    if(length(which(groupdat[,k]==3|groupdat[,k]==4))>=grpsize-1){
      rownombres<-c(rownombres,as.integer(names(groupdat[which(groupdat[,k]==3|groupdat[,k]==2),k])),999999999)
      rowID<-c(rowID,rep(p,each=length(names(groupdat[which(groupdat[,k]==3|groupdat[,k]==2),k]))),999999999)
      p=p+1
    }
  }
  if(!is.null(rowID)){
    return(as.matrix(rbind(rownombres,rowID)))
  }
   #}
}

sox.write<-function(numPass,m,b,pathh,sound_filesfullpathB,combSound,durTab,pad,bigFile_breaks,sound_filesB,durTab2){
dir.create(paste(pathh,sep=""))
print(paste("Creating file ",MoorInfo[m,10],"_",bigFile_breaks[b],sep=""))
sox_alt(paste(noquote(paste(paste(sound_filesfullpathB,collapse=" ")," ",combSound,sep=""))),exename="sox.exe",path2exe=paste(drivepath,"Accessory/sox-14-4-2",sep=""))
durList<-duration_store(numPass,m,sound_filesfullpathB,durTab,pad,sound_filesB,durTab2)
if(numPass==1){
  return(durList[[1]])
}else if(numPass==2){
  return(durList[[2]])
}

}

write_4byte_unsigned_int <- function(x, con){
  if((!is.numeric(x) || is.na(x)) || (x < 0 || x > (2^32-1))) 
    stop("the field length must be a 4 byte unsigned integer in [0, 2^32-1], i.e. file size < 4 GB")
  big <- x %/% 256^2
  small <- x %% 256^2
  writeBin(as.integer(c(small, big)), con, size = 2, endian = "little")
}

write_longvector <- function(x, ..., bytes){
  if(bytes > (2^31-1)){
    index <- 1:floor(length(x)/2)
    writeBin(x[index], ...)
    writeBin(x[-index], ...)        
  }
  else writeBin(x, ...)
}



writeWave.nowarn <- 
  function(object, filename, extensible = TRUE) {
    if(!is(object, "WaveGeneral")) 
      stop("'object' needs to be of class 'Wave' or 'WaveMC'")
    validObject(object)
    if(is(object, "Wave")){
      object <- as(object, "WaveMC")
      colnames(object) <- c("FL", if(ncol(object) > 1) "FR")
    }
    if(ncol(object) > 2 && !extensible)
      stop("Objects with more than two columns (multi channel) can only be written to a Wave extensible format file.")
    
    cn <- colnames(object)
    if((length(cn) != ncol(object) || !all(cn %in% MCnames[["name"]])) || any(duplicated(cn)))
      stop("colnames(object) must be specified and must uniquely identify the channel ordering for WaveMC objects, see ?MCnames for possible channels")
    cnamesnum <- as.numeric(factor(colnames(object), levels=MCnames[["name"]]))
    if(is.unsorted(cnamesnum))
      object <- object[,order(cnamesnum)]
    dwChannelMask <- sum(2^(cnamesnum - 1))  ##
    
    l <- as.numeric(length(object)) # can be an int > 2^31
    sample.data <- t(object@.Data)
    dim(sample.data) <- NULL
    
    ## PCM or IEEE FLOAT
    pcm <- object@pcm                                 
    
    # Open connection
    con <- file(filename, "wb")
    on.exit(close(con)) # be careful ...
    
    # Some calculations:
    byte <- as.integer(object@bit / 8)
    channels <- ncol(object)
    block.align <- channels * byte
    bytes <- l * byte * channels
    
    if((!is.numeric(bytes) || is.na(bytes)) || (bytes < 0 || (round(bytes + if(extensible) 72 else 36) + 4) > (2^32-1)))
      stop(paste("File size in bytes is", round(bytes + if(extensible) 72 else 36) + 4, "but must be a 4 byte unsigned integer in [0, 2^32-1], i.e. file size < 4 GB"))
    
    ## Writing the header:
    # RIFF
    writeChar("RIFF", con, 4, eos = NULL) 
    write_4byte_unsigned_int(round(bytes + if(extensible) 72 else 36), con) # cksize RIFF
    
    
    # WAVE
    writeChar("WAVE", con, 4, eos = NULL)
    # fmt chunk
    writeChar("fmt ", con, 4, eos = NULL)
    if(extensible) { # cksize format chunk
      writeBin(as.integer(40), con, size = 4, endian = "little") 
    } else {
      writeBin(as.integer(16), con, size = 4, endian = "little")
    }    
    if(!extensible) { # wFormatTag
      writeBin(as.integer(if(pcm) 1 else 3), con, size = 2, endian = "little")
    } else {
      writeBin(as.integer(65534), con, size = 2, endian = "little") # wFormatTag: extensible   
    }
    writeBin(as.integer(channels), con, size = 2, endian = "little") # nChannels
    writeBin(as.integer(object@samp.rate), con, size = 4, endian = "little") # nSamplesPerSec
    writeBin(as.integer(object@samp.rate * block.align), con, size = 4, endian = "little") # nAvgBytesPerSec
    writeBin(as.integer(block.align), con, size = 2, endian = "little") # nBlockAlign
    writeBin(as.integer(object@bit), con, size = 2, endian = "little") # wBitsPerSample
    # extensible
    if(extensible) {
      writeBin(as.integer(22), con, size = 2, endian = "little") # cbsize extensible
      writeBin(as.integer(object@bit), con, size = 2, endian = "little") # ValidBitsPerSample
      writeBin(as.integer(dwChannelMask), con, size = 4, endian = "little") #  dbChannelMask
      writeBin(as.integer(if(pcm) 1 else 3), con, size = 2, endian = "little") # SubFormat 1-2
      writeBin(as.raw(c(0,   0,   0,  0,  16,   0, 128,   0 ,  0, 170,   0,  56, 155, 113)), con) # SubFormat 3-16
      # fact
      writeChar("fact", con, 4, eos = NULL)
      writeBin(as.integer(4), con, size = 4, endian = "little") # cksize fact chunk
      writeBin(as.integer(l), con, size = 4, endian = "little") # dwSampleLength
    }
    # data
    writeChar("data", con, 4, eos = NULL)
    write_4byte_unsigned_int(round(bytes), con)
    
    # Write data:
    # PCM format
    if(pcm) { 
      if(byte == 3){
        sample.data <- sample.data + 2^24 * (sample.data < 0)
        temp <- sample.data %% (256^2)
        sample.data <- sample.data %/% 256^2
        a2 <- temp %/% 256
        temp <- temp %%  256
        write_longvector(as.integer(rbind(temp, a2, sample.data)), con, size = 1, endian = "little", bytes=bytes)
      } else {
        write_longvector(as.integer(sample.data), con, size = byte, endian = "little", bytes=bytes)
      }
    } else {
      write_longvector(as.numeric(sample.data), con, size = byte, endian = "little", bytes=bytes)
    }
    
    invisible(NULL)
  }

prime.factor <- function(x){
  n=c()
  i=2
  r=x
  while(prod(n)!=x){
    if(!r%%i) {n=c(n,i);r=r/i;i=1}
    i=i+1
  }
  return(n)
}

freqstat.normalize<- function(freqstat,lowFreq,highFreq){
  newstat<-((freqstat-lowFreq)/(highFreq-lowFreq))
  return(newstat)
}

raven_batch_detec <- function(raven.path = NULL, sound.files, path = NULL, detector = "Amplitude detector", relabel_colms = TRUE, pb = TRUE, dpreset="Default",vpreset="Default")
{
  
  #check path to working directory
  if (is.null(path)) path <- getwd() else if (!dir.exists(path)) stop("'path' provided does not exist") 
  
  # reset working directory 
  wd <- getwd()
  on.exit(setwd(wd), add = TRUE)
  on.exit(file.remove(file.path(raven.path, "temp.bcv.txt")), add = TRUE)
  
  if (is.null(raven.path))
    stop("Path to 'Raven' folder must be provided")  else
      if (!dir.exists(raven.path)) stop("'raven.path' provided does not exist")
  
  setwd(raven.path)
  
  sf <- sound.files <- as.character(sound.files)
  
  #return warning if not all sound files were found
  recs.wd <- list.files(path = path, pattern = "\\.wav$|\\.aif$|\\.flac$|\\.mp3$", ignore.case = TRUE)
  
  #count number of sound files in working directory and if 0 stop
  sound.files <- sound.files[sound.files %in% recs.wd]
  if (length(sound.files) == 0)
    stop("The sound files are not in the working directory")
  
  # remove sound files not found
  if (length(sound.files) != length(sf)) 
    cat(paste(length(sf) - length(sound.files), ".wav file(s) not found"))
  
  
  # check if sound file names contains directory and fix
  if (basename(sound.files[1]) == sound.files[1])
    sound.files <- file.path(path, sound.files)
  
  if (pb) lply <- pbapply::pblapply else lply <- lapply
  
  out <- lply(sound.files, function(x) {
    
    if (Sys.info()[1] == "Windows")
    {  
      comnd <- paste(shQuote(file.path(raven.path, "Raven.exe"), type = "cmd"), paste0("-detPreset:",dpreset), paste0("-viewPreset:",vpreset), paste0("-detType:", detector), shQuote(x), "-detTable:temp.bcv.txt -x")
    } else
    {
      if (Sys.info()[1] == "Linux")
        comnd <- paste(file.path(raven.path, "Raven"), paste0("-detPreset:",dpreset), paste0("-viewPreset:",vpreset), paste0("-detType:", detector), x, "-detTable:temp.bcv.txt -x") else
          comnd <- paste("open Raven.app --args", x, paste0("-detPreset:",dpreset), paste0("-viewPreset:",vpreset), paste0("-detType:", detector), "-detTable:temp.bcv.txt -x") # OSX
    }
    
    # run raven
    system(command = comnd, ignore.stderr = TRUE, intern = TRUE)
    
    output <- utils::read.table("temp.bcv.txt", sep = "\t",  header = TRUE)
    
    if (nrow(output) > 1)
      output$sound.files <- basename(x) else output <- vector(length = 0)
    
    return(output)
  })
  
  out <- out[sapply(out, is.data.frame)]
  
  if (length(out) > 0)
  {
    output <- do.call(rbind, out)
    
    if (relabel_colms)
      output <- relabel_colms(output, waveform = !any(grepl("Spectrogram", output$View)))
    
    output <- output[, c(ncol(output), 2:(ncol(output) - 1))]  
    
    return(output)} else return(NULL)
  
  
}

plot_jpeg = function(path, add=FALSE)
{
  require('jpeg')
  jpg = readJPEG(path, native=T) # read the file
  res = dim(jpg)[2:1] # get the resolution, [x, y]
  if (!add) # initialize an empty plot area if add==FALSE
    plot(1,1,xlim=c(1,res[1]),ylim=c(1,res[2]),asp=1,type='n',xaxs='i',yaxs='i',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  rasterImage(jpg,1,1,res[1],res[2])
}

sox_alt <- function (command, exename = NULL, path2exe = NULL, argus = NULL, shQuote_type = NULL)
{
  
  if (is.null(exename)) {
    exename <- "sox"
  } else {
    exename = exename
  }
  
  if (is.null(path2exe)) {
    exe <- exename
  } else {
    path2exe <- normalizePath(path = path2exe, winslash = "/", mustWork = TRUE)
    exe <- paste(path2exe, exename, sep = "/")
  }
  
  if (.Platform$OS.type == "windows") {
    
    if(is.null(shQuote_type)) {
      shQuote_type = "cmd"
    } else {
      shQuote_type = "cmd2"
    }
    
  } else {  # .Platform$OS.type == "unix" # + Apple OS X + other
    
    if(is.null(shQuote_type)) {
      shQuote_type = "sh"
    } else {
      shQuote_type = "csh"
    }
    
  }
  
  exe <- shQuote(exe, type = shQuote_type)
  system(paste(exe, command, argus, sep = " "), ignore.stderr = TRUE)
  
}

process_model<-function(stuff2,Moddata2){

  if(modelType=="rf"){
  giniTab<-data.frame(stuff2[[1]][[2]])
  for(i in 2:CV){
    giniTab<-cbind(giniTab,stuff2[[i]][[2]])
  }
    
  #for some reason have to save giniAv as global variable to reassign rownames...
  giniAv<<-data.frame(apply(giniTab,1,mean))
  giniRows<<-c(colnames(Moddata2[,2:ncol(Moddata2)]))[which(!c(colnames(Moddata2[,2:ncol(Moddata2)])) %in% c('detectionType','MooringCode','year','month'))]
  rownames(giniAv)<<-giniRows
  colnames(giniAv)<<-"MeanDecreaseGini"
    
  varImpPlot_AVG(giniAv,  
                 sort = TRUE,
                 main="Variable Importance random forest")
  giniAv<-NULL
  }
  
  
  #probDeets<-NULL
  CUTmean<-list()
  
  for(s in 1:length(mSpec)){

  CUTvec<-NULL
      
  for(i in 1:CV){
    CUTvec<-c(CUTvec,stuff2[[i]][[3]][[s]])
  }
  
  CUTmean[[s]]<-mean(unlist(CUTvec))
  
  #build probability data frame 
  probstab<-NULL
  for(i in 1:CV){
    probstab<-cbind(probstab,stuff2[[i]][[1]][[s]][,2])
  }

  probmean<-NULL
  probstderr<-NULL
  n<-NULL
  for(x in 1:nrow(probstab)){
    probmean[x]<-mean(as.numeric(probstab[x,2:ncol(probstab)]),na.rm=TRUE)
    probstderr[x]<-std.error(as.numeric(probstab[x,2:ncol(probstab)]),na.rm=TRUE)
    n[x]<-length(which(!is.na(probstab[x,2:ncol(probstab)])))
  }
  
  ##assuming $detection type is already in this data NOTE not 
  tempNames<-colnames(Moddata2)
  Moddata2<-cbind(Moddata2,probmean,probstderr,n)
  colnames(Moddata2)<-c(tempNames,paste(mSpec[s],"prob"),paste(mSpec[s],"stderr"),paste(mSpec[s],"n"))
  
  #tempNames<-colnames(probDeets)
  #probDeets<-cbind(probDeets,probmean)
  #probDeets<-cbind(probDeets,probstderr)
  #probDeets<-cbind(probDeets,n)
    
  #colnames(probDeets)<-c(tempNames,paste(spec[s],"prob"),paste(spec[s],"stderr"),paste(spec[s],"n"))
  }
  
  return(list(Moddata2,CUTmean))
  probDeets<<-NULL
}

runObliqueRandomForest<-function(Moddata,method){
  
  Moddata<<-Moddata
  
  if(length(unique(Moddata$detectionType))==2){
    
    #appears oblique RF only supports binary classification and doesn't allow formula style input. Need to take row I don't want to be predictors out before fitting model. 
    
    print(paste("creating oblique random forest models with CV",CV))
    
    if(parallelType=="local"){
    startLocalPar("Moddata","CV","splitdf","TPRthresh")
    }
    
    stuff<-foreach(p=1:CV,.packages=c("obliqueRF","ROCR","stats")) %dopar% {
    train<-splitdf(Moddata,weight = 2/3)
    trainModdataPred<-as.matrix(train[[1]][,c(8,9,11:ncol(train[[1]]))])
    testModdataPred<-as.matrix(train[[2]][,c(8,9,11:ncol(train[[1]]))])
    trainModdataResponse<-as.numeric(train[[1]][,7])
    
    
    data.orf<-obliqueRF(x=trainModdataPred,y=trainModdataResponse,mtry=11,training_method=method) 
    
    #stats::predict?
    pred<-predict(data.orf,testModdataPred,type="prob")
    pred<-cbind(pred,train[[2]]$Selection)

    ROCRpred<-prediction(pred[,2],train[[2]]$detectionType)
    prob.perf = performance(ROCRpred, "tpr","fpr")
    
    TPR<-NULL
    TPR <- data.frame(cut=prob.perf@alpha.values[[1]], tpr=prob.perf@y.values[[1]])
    CUT <- max(TPR[which(TPR$tpr>=TPRthresh),1])
    
    probstab<-data.frame(Moddata$Selection)
    probstab[,2]<-NA
    #giniTab<-as.numeric(data.orf$importance)  spit out nonsense for ridge. Assuming this will not work for most methods so ignoring for now. 
    for(n in 1:nrow(pred)){
      probstab[which(probstab$Moddata.Selection==pred[n,3]),2]<-pred[n,2]
    }
    placeHolder<-TRUE
    return(list(probstab[,2],placeHolder,CUT))
    
    }
    
  if(parallelType=="local"){
  parallel::stopCluster(cluz)
  }

}else{
    
   #do other one (need up update full mooring on this)
}
  
  stuff<-process_model(stuff,Moddata)
  return(stuff)
  Moddata<<-NULL
}

runRavenDetector<-function(m,filePath,combname,resltsTab){
  for(r in detectorssprshort){
    print(paste("Running detector for",MoorInfo[m,9],MoorInfo[m,10]))
    resltVar <- raven_batch_detec(raven.path = ravenpath, sound.files = combname, path = filePath ,detector = "Band Limited Energy Detector",dpreset=r,vpreset=ravenView)
    resltVar$MooringID<-MoorInfo[m,10]
    resltVar$MooringName<-MoorInfo[m,1]
    resltVar$MooringCode<-sub("[:alpha:]$", "",unlist(strsplit(sub("(_)(?=[^_]+$)", " ", MoorInfo[m,1], perl=T), " "))[c(FALSE, TRUE)])
    resltVar$detector<-r
    resltVar$Species<-MoorInfo[m,9]
    if(is.null(nrow(resltVar))==FALSE){
      resltsTab<- rbind(resltsTab,resltVar)
    }
  resltVar<-NULL 
  }
  return(resltsTab)
}

runRandomForest<-function(Moddata,GTdataset){
  
  Moddata<<-Moddata
  
  GTdataset<<-GTdataset
  
  print(paste("creating random forest models with CV",CV))

  if(parallelType=="local"){
  startLocalPar("Moddata","CV","splitdf","TPRthresh","mSpec","whichRun")
  }
  
  stuff<-foreach(p=1:CV,.packages=c("randomForest","ROCR","stats")) %dopar% {

  train<-splitdf(GTdataset,weight = 2/3)
  data.rf<-randomForest(formula=detectionType ~ . -Selection -MooringCode -year -month,data=train[[1]],mtry=11) #-meanfreq,-freqrange,   na.action = na.roughfix

  #if(length(unique(Moddata$detectionType))>1){
  #train[[2]]<-splitdf()
  
  pred<-stats::predict(data.rf,train[[2]],type="prob")
  pred<-cbind(pred,train[[2]]$Selection)
  
  if(whichRun=="NEW"){
  predreal<-stats::predict(data.rf,Moddata,type="prob")
  predreal<-cbind(predreal,Moddata$Selection)
  }
  
  p=1
  probstab<-list()
  CUT<-list()
  giniTab<-as.numeric(data.rf$importance)
    
  for(s in (2*(1:length(mSpec))-2)){
    
  oneSpecDetType<-as.numeric(train[[2]]$detectionType)-(1+s)
  
  oneSpecDetType[which(oneSpecDetType>1|oneSpecDetType<1)]<-0
                 
  ROCRpred<-prediction(pred[,s+2],oneSpecDetType)
  prob.perf = performance(ROCRpred, "tpr","fpr")
  
  TPR<-NULL
  TPR <- data.frame(cut=prob.perf@alpha.values[[1]], tpr=prob.perf@y.values[[1]])
  CUT <- c(CUT,max(TPR[which(TPR$tpr>=TPRthresh[which(s %in% (2*(1:length(mSpec))-2))]),1]))
  
  probstab[[p]]<-data.frame(Moddata$Selection)
  probstab[[p]][,2]<-NA
  
  if(whichRun=="GT"){
  for(n in 1:nrow(pred)){
    probstab[[p]][which(probstab[[p]]$Moddata.Selection==pred[n,(length(mSpec)*2+1)]),2]<-pred[n,2+s]
  }
  }else if(whichRun=="NEW"){
    for(n in 1:nrow(predreal)){
      probstab[[p]][which(probstab[[p]]$Moddata.Selection==predreal[n,(length(mSpec)*2+1)]),2]<-predreal[n,2+s]
    }
  }
  
  
  p=p+1
  }
  
  return(list(probstab,giniTab,CUT))
}
  if(parallelType=="local"){
  parallel::stopCluster(cluz)
  }

  stuff<-process_model(stuff,Moddata)
  return(stuff)
  Moddata<<-NULL
}

context_sim <-function(sdata){
  CUTmeanspec<-CUTmean[[s]]
  datTab<-matrix(,ncol=8,nrow=0)
  colnames(datTab)<-paste(spec[s],c("Selection","mRT","Oprob","Fmod","Fprob","Bmod","Bprob","Tprob"))
  #context simulator- add or subtract % points based on how good neighboring calls were. Only useful for full mooring dataset. 
  for(w in 1:length(unique(sdata[,6]))){
    datVar<-sdata[which(sdata[,6]==unique(sdata[,6])[w]),]
    mRT<-datVar[,13]+((datVar[,10]+datVar[,11])/2)
    datVar<-cbind(datVar[,1],mRT,datVar[,nospec+s])
    datVar<-cbind(datVar,matrix(0,nrow=nrow(datVar),ncol=5))
    datVar[,5]<-datVar[,3]
    for(n in 1:(nrow(datVar)-1)){
      if(datVar[n,5]>=greatcallThresh){
        datVar[n+1,4]<-maxBonus
        if(datVar[n+1,4]>maxBonus){
          datVar[n+1,4]<-maxBonus
        }
      }else if(datVar[n,5]>=(-maxPenalty+CUTmeanspec)&datVar[n,5]<greatcallThresh){
        datVar[n+1,4]<-datVar[n,4]+(datVar[n,5]*goodcallBonus)
        if(datVar[n+1,4]>maxBonus){
          datVar[n+1,4]<-maxBonus
        }
      }else{
        datVar[n+1,4]<-datVar[n,4]+((mRT[n+1]-mRT[n])*badcallPenalty) #*datVar[n,3])
        if(datVar[n+1,4]<maxPenalty){
          datVar[n+1,4]<-maxPenalty
        }
      }
      if(datVar[n,5]+datVar[n,4]>0 & datVar[n,5]+datVar[n,4]<1){
        datVar[n,5]<-datVar[n,5]+datVar[n,4]
      }else if(datVar[n,5]+datVar[n,4]<0){
        datVar[n,5]<-0
      }else if(datVar[n,5]+datVar[n,4]>1){
        datVar[n,5]<-1
      }
    }
  #same but backwards through data 
  datVar[,7]<-datVar[,3]
  for(n in (nrow(datVar)-1):2){
    if(datVar[n,7]>=greatcallThresh){
      datVar[n-1,6]<-maxBonus
      if(datVar[n-1,6]>maxBonus){
        datVar[n-1,6]<-maxBonus
      }
    }else if(datVar[n,7]>(-maxPenalty)&datVar[n,7]<greatcallThresh){
      datVar[n-1,6]<-datVar[n,6]+(datVar[n,7]*goodcallBonus)
      if(datVar[n-1,6]>maxBonus){
        datVar[n-1,6]<-maxBonus
      }
    }else{
      datVar[n-1,6]<-datVar[n,6]+((mRT[n]-mRT[n-1])*badcallPenalty)
      if(datVar[n-1,6]<maxPenalty){
        datVar[n-1,6]<-maxPenalty
      }
    }
    if(datVar[n,7]+datVar[n,6]>0 & datVar[n,7]+datVar[n,6]<1){
      datVar[n,7]<-datVar[n,7]+datVar[n,6]
    }else if(datVar[n,7]+datVar[n,6]<0){
      datVar[n,7]<-0
    }else if(datVar[n,7]+datVar[n,6]>1){
      datVar[n,7]<-1
    }
  }
  
  datVar[,8]<-pmax(datVar[,5],datVar[,7]) #max
  #datVar$probmean<-pmax(datVar[n,4],datVar[n,pos+3]) #max
  #datVar$probmean<-pmin(datVar[n,4],datVar[n,pos+3]) #min
  
  datTab<-rbind(datTab,datVar)
  }
  return(datTab)
}

#should do table for all species/mooring, and for each species/mooring
after_model_write <-function(mdata){
  #all species/mooring 

  #in theory, could be a conflict here if two seperate species HG data share the same name by coincedence. 
  for(v in 1:length(unique(MoorInfo[,10]))){
    print(paste("first loop",v))
    MoorVar1<-mdata[which(mdata$MooringID == unique(MoorInfo[,10])[v]),]
    #only cut out detections where all species are below the threshold. 
    
    MoorVar1$remove<-0
    for(s in 1:length(spec)){
      for(n in 1:nrow(MoorVar1)){
        if(MoorVar1[n,nospec+s]<CUTmean[[s]])
        MoorVar1$remove[n]<-MoorVar1$remove[n]+1
      }
    }
    
    MoorVar1<-MoorVar1[which(MoorVar1$remove!=length(mSpec)),]
    MoorVar1$remove<-NULL
    
    MoorVar1<-MoorVar1[order(MoorVar1$Begin.Time..s.),]
    #dont think we need to report stats for these (total # of sucessful classifications). Pretty meaningless when you pool all of them together. Just print the table
    RavenExport<-data.frame(seq(1,nrow(MoorVar1),by=1))
    RavenExport[,2]<-"Spectrogram 1"
    RavenExport[,3]<-1
    RavenExport[,4]<-MoorVar1[,2]
    RavenExport[,5]<-MoorVar1[,3]
    RavenExport[,6]<-MoorVar1[,4]
    RavenExport[,7]<-MoorVar1[,5]
    RavenExport[,8]<-MoorVar1$detectionType
    RavenExport[,9]<-MoorVar1$File
    RavenExport[,10]<-MoorVar1$FileOffsetBegin
    for(s in 1:length(mSpec)){
      RavenExport<-cbind(RavenExport,MoorVar1[,nospec+s])
    }
    
    colnames(RavenExport)<-c("Selection","View","Channel","Begin Time (s)","End Time (s)","Low Freq (Hz)","High Freq (Hz)","detType","File","FileOffsetBegin",paste(mSpec,"prob"))
    
    write.table(RavenExport,paste(outputpath,runname,"/",unique(MoorInfo[,10])[v],"_AllSpecies_Model_Applied_probs",".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
  }

  for(v in 1:length(unique(MoorInfo[,10]))){
    print(paste("2nd loop",v))
    
  #now do it for each species individually and each unique mooring ID. Compare with GTtot if available. 
    for(s in 1:length(mSpec)){
      MoorVar1<-mdata[which(mdata$MooringID == unique(MoorInfo[,10])[v]),]
      
      MoorVar1<-MoorVar1[which(MoorVar1[,nospec+s]>=CUTmean[[s]]),]
      
      if(nrow(MoorVar1)>0){
      
      MoorVar1<-MoorVar1[order(MoorVar1$Begin.Time..s.),]
      #save data frame 
      RavenExport<-data.frame(seq(1,nrow(MoorVar1),by=1))
      RavenExport[,2]<-"Spectrogram 1"
      RavenExport[,3]<-1
      RavenExport[,4]<-MoorVar1[,2]
      RavenExport[,5]<-MoorVar1[,3]
      RavenExport[,6]<-MoorVar1[,4]
      RavenExport[,7]<-MoorVar1[,5]
      RavenExport[,8]<-MoorVar1$detectionType
      RavenExport[,9]<-MoorVar1$File
      RavenExport[,10]<-MoorVar1$FileOffsetBegin
      RavenExport[,11]<-MoorVar1[,nospec+s]
      
      
      colnames(RavenExport)<-c("Selection","View","Channel","Begin Time (s)","End Time (s)","Low Freq (Hz)","High Freq (Hz)","detType","File","FileOffsetBegin",paste(mSpec[s],"prob"))
      
      write.table(RavenExport,paste(outputpath,runname,"/",unique(MoorInfo[,10])[v],"_",mSpec[s],"_Model_Applied_probs",".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
      
      #get together variables to write to detectorRunLog
      if(modelType=="orf"){
        name<-paste(modelType,modelMethod,sep=":")
      }else if(modelType=="rf"){
        name<-modelType
      }
      
      #calculate stats: 
      detTotal<-nrow(MoorVar1)
      
      #check to see if TPtottab entry exists (if so, compare it and calculate performance. If not, treat as new data and do not calculate performance) 
      if(any(paste(mSpec[s],unique(MoorInfo[,10])[v])==TPtottab$MoorCor)){
        numTP<-sum(MoorVar1$detectionType==paste(mSpec[s],"1"))
        numTPtruth<-TPtottab[which(paste(mSpec[s],unique(MoorInfo[,10])[v])==TPtottab$MoorCor),2]
        numFP<-detTotal-numTP
        numFN<-numTPtruth-numTP
        
        TPR <- numTP/(numTP+numFN)
        FPR <- numFP/(numTP+numFP)
        TPdivFP<- numTP/numFP
        
        detecEval<-read.csv(paste(outputpath,"DetectorRunLog.csv",sep=""))
        detecEval<-detecEval[0,]
        detecEval<-data.frame(lapply(detecEval,as.character),stringsAsFactors = FALSE)

        detecEval[1,]<-c(mSpec[[s]],paste(name,MoorVar1[1,"MooringName"]),paste(detnum,paste(detlist2,collapse="+"),sep=";"),"spread",runname,detTotal,numTP,numFP,numFN,TPR,FPR,TPdivFP,AUCadj,paste(CV,TPRthresh,sep=","),paste(greatcallThresh,-maxPenalty,sep=","),paste(maxBonus,goodcallBonus,badcallPenalty,sep=","), paste(allowedZeros,collapse=","),paste(grpsize,collapse=","),paste(downsweepCompMod,downsweepCompAdjust,sep=","),paste(detskip,collapse=","),paste(groupInt,collapse=","),NA,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")

        
        detecEval2<-read.csv(paste(outputpath,"DetectorRunLog.csv",sep=""))
        detecEvalFinal<-rbind(detecEval2,detecEval)
        write.csv(detecEvalFinal,paste(outputpath,"DetectorRunLog.csv",sep=""),row.names=FALSE)
      }else{
        numTP<-"uk"
        numTPtruth<-"uk"
        numFP<-"uk"
        numFN<-"uk"
        
        TPR <- "uk"
        FPR <- "uk"
        TPdivFP<-"uk"
        
        detecEval<-read.csv(paste(outputpath,"DetectorRunLog.csv",sep=""))
        detecEval<-detecEval[0,]
        detecEval<-data.frame(lapply(detecEval,as.character),stringsAsFactors = FALSE)
        
        detecEval[1,]<-c(mSpec[[s]],paste(name,MoorVar1[1,"MooringName"]),paste(detnum,paste(detlist2,collapse="+"),sep=";"),"spread",runname,detTotal,numTP,numFP,numFN,TPR,FPR,TPdivFP,AUCadj,paste(CV,TPRthresh,sep=","),paste(greatcallThresh,-maxPenalty,sep=","),paste(maxBonus,goodcallBonus,badcallPenalty,sep=","), paste(allowedZeros,collapse=","),paste(grpsize,collapse=","),paste(downsweepCompMod,downsweepCompAdjust,sep=","),paste(detskip,collapse=","),paste(groupInt,collapse=","),NA,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")
        
        detecEval2<-read.csv(paste(outputpath,"DetectorRunLog.csv",sep=""))
        detecEvalFinal<-rbind(detecEval2,detecEval)
        write.csv(detecEvalFinal,paste(outputpath,"DetectorRunLog.csv",sep=""),row.names=FALSE)
      }
    }
    }
  }
  for(s in 1:length(mSpec)){
  print(paste("3rd loop",s))
    
    
  MoorVar1<-mdata[which(mdata[,nospec+s]>=CUTmean[[s]]),]
  
  detTotal<-nrow(MoorVar1)
  
  MoorVar1 <-MoorVar1[which(paste(mSpec[[s]],MoorVar1$MooringID) %in% TPtottab$MoorCor),]
    
    numTP<-sum(MoorVar1$detectionType==paste(mSpec[[s]],"1"))
    numTPtruth<-sum(TPtottab[which(TPtottab$MoorCor %in% paste(mSpec[[s]],unique(MoorVar1$MooringID))),2])
    numFP<-detTotal-numTP
    numFN<-numTPtruth-numTP
    
    TPR <- numTP/(numTP+numFP)
    FPR <- numFP/(numTP+numFP)
    TPdivFP<- numTP/numFP
    
    detecEval<-read.csv(paste(outputpath,"DetectorRunLog.csv",sep=""))
    detecEval<-detecEval[0,]
    detecEval<-data.frame(lapply(detecEval,as.character),stringsAsFactors = FALSE)
    
    detecEval[1,]<-c(mSpec[[s]],paste(name,"all"),paste(detnum,paste(detlist2,collapse="+"),sep=";"),"spread",runname,detTotal,numTP,numFP,numFN,TPR,FPR,TPdivFP,AUCadj,paste(CV,TPRthresh,sep=","),paste(greatcallThresh,-maxPenalty,sep=","),paste(maxBonus,goodcallBonus,badcallPenalty,sep=","), paste(allowedZeros,collapse=","),paste(grpsize,collapse=","),paste(downsweepCompMod,downsweepCompAdjust,sep=","),paste(detskip,collapse=","),paste(groupInt,collapse=","),NA,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")
    
    
    detecEval2<-read.csv(paste(outputpath,"DetectorRunLog.csv",sep=""))
    detecEvalFinal<-rbind(detecEval2,detecEval)
    write.csv(detecEvalFinal,paste(outputpath,"DetectorRunLog.csv",sep=""),row.names=FALSE)
  
  }
}

adaptive_compare<-function(Compdata){
  for(a in 1:10){#go through twice in case there are mulitple boxes close to one another. 
  for(o in unique(Compdata[,6])){
    print(paste("for mooring",s,o))
    CompVar<-Compdata[which(Compdata[,6]==o),]
    CompVar<-CompVar[order(CompVar[,2]),]
    n=0
    newrow<-matrix(0,ncol=nospec+length(spec)+1)
    IDvec<-NULL
    for(q in 1:(nrow(CompVar)-1)){
      if(CompVar[q+1,2]<=(CompVar[q,2]+timediffself)){
        if(CompVar[q+1,nospec+s]+probdist<CompVar[q,nospec+s]|CompVar[q+1,nospec+s]-probdist>CompVar[q,nospec+s]){#take only the best one
          newdat<-CompVar[c(q,q+1),]
          newdat<-newdat[order(newdat[,nospec+s]),]
          IDvec<-c(IDvec,as.numeric(newdat[,1]))
          newrow<-rbind(newrow,newdat[2,])
        }else{
          newdat<-CompVar[c(q,q+1),]
          IDvec<-c(IDvec,as.numeric(newdat[,1]))
          sel1<-newdat[1,1]
          strt<-as.numeric(min(newdat[,2]))
          e<-as.numeric(max(newdat[,3]))
          l<-as.numeric(min(newdat[,4]))
          h<-as.numeric(max(newdat[,5]))
          MID<-newdat[1,6]
          Mname<-newdat[1,7]
          Sf<-newdat[1,8]
          File<-newdat[1,9]
          FOB<-newdat[1,10]
          FOE<-newdat[1,11]
          RTf<-newdat[1,12]
          RTb<-newdat[1,13]
          RTe<-newdat[1,14]
          species<-newdat[1,15]
          sel2<-newdat[1,16]
          mf<-(h+l)/2
          fr<-(h-l)
          row<-c(sel1,strt,e,l,h,MID,Mname,Sf,File,FOB,FOE,RTf,RTb,RTe,species,sel2,mf,fr)
          
          Mcode<-newdat[1,nospec-3]
          dt<-max(newdat[,nospec-2])
          yr<-newdat[1,nospec-1]
          month<-newdat[1,nospec]
          
          afterrow<-c(Mcode,dt,yr,month)
          
          newSpec<-NULL
          for(i in 1:length(spec)){
          newSpec<-c(newSpec,mean(newdat[,nospec+i]))
          }
          
          print(paste("adaptive compare index:",newdat[1,nospec+length(spec)+1]))
          print(paste("adaptive compare lookup:",which(newdat[1,nospec+length(spec)+1]==as.numeric(factor(paste(MoorInfo[,9],MoorInfo[,10]))))))
          print(paste("adaptive compare MoorID:",MoorInfo[which(newdat[1,nospec+length(spec)+1]==as.numeric(factor(paste(MoorInfo[,9],MoorInfo[,10])))),10]))
            
          #print(c(newdat[1,nospec+length(spec)+1],row[c(2,3,4,5)],(as.numeric(RTb)+as.numeric(FOB)),dt))
          row2<-unlist(spectral_features(c(newdat[1,nospec+length(spec)+1],row[c(2,3,4,5)],(as.numeric(RTb)+as.numeric(FOB)),dt)))
          
          row<-c(row,row2[c(8:length(row2))],afterrow,newSpec,newdat[1,nospec+length(spec)+1])
          
          newrow<-rbind(newrow,row)
          if(newrow[1,1]==0){
            newrow<-matrix(newrow[2,],ncol=length(newrow[2,]))
          }
        }
          
        }

      }
      if(newrow[1,1]!=0){
        CompVar<-subset(CompVar,!(CompVar[,1] %in% IDvec))
        CompVar<-rbind(CompVar,newrow)
        CompVar<-CompVar[order(CompVar[,2]),]
        n=n+1
      }
    if(n>0){
      Compdata<-Compdata[-which(Compdata[,6]==o),]
      Compdata<-rbind(Compdata,CompVar)      
    }

  }
  }
 return(Compdata) 
}

  
duration_store<- function(numPass,m,sound_filesfullpathB,durTab,pad,sound_filesB,durTab2){
  if(numPass==1){
  durVar<-NULL
  durVar<-data.frame(SFfp=character(),
                     SFsh=character(),
                     Duration=numeric(), 
                     CumDur=numeric(),
                     CombSF=character(), 
                     MooringID=character(),
                     MooringName=character(),
                     stringsAsFactors=FALSE) 
for(f in 1:length(sound_filesfullpathB)){
  durVar[f,1]<-sound_filesfullpathB[f]
  durVar[f,2]<-sound_filesB[f]
  audio<-readWave(sound_filesfullpathB[f], header=TRUE)
  durVar[f,3]<-round(audio$samples / audio$sample.rate, 2)
}
  durVar[,4]<-cumsum(durVar$Duration)
  durVar[,5]<-paste(pad,MoorInfo[m,10],".wav",sep="")
  durVar[,6]<-MoorInfo[m,10]
  durVar[,7]<-MoorInfo[m,1]

durTab<-rbind(durTab,durVar)
durTab2<-1 #workaround to R not allowing dataframes to be "null". 
  }else if(numPass==2){
    durVar<-durTab[which(durTab$CombSF %in% sound_filesB),]
   # durVar<-durVar[-which(duplicated(durVar$SFsh)),]
    durVar$CombSF<-paste(pad,MoorInfo[m,10],".wav",sep="")
    durVar$CumDur<-cumsum(durVar$Duration)
    durTab2<-rbind(durTab2,durVar)
  }
durList<-list(durTab,durTab2)  
  return(durList)
}

spectral_features<- function(specdata){
  
  specdata<<-specdata

  specpath<<-NULL
  specTab<<-NULL
  
  if(is.null(nrow(specdata))){
  index<-which(specdata[1]==as.numeric(factor(paste(MoorInfo[,9],MoorInfo[,10]))))
  moors<-paste(MoorInfo[index,9],MoorInfo[index,10])
  }else{
    moors<-paste(MoorInfo[,9],MoorInfo[,10])
  }
  
for(m in moors){
  
  loadSpecVars(MoorInfo[which(paste(MoorInfo[,9],MoorInfo[,10])==m),9])
  print(paste("spectral features index:",which(paste(MoorInfo[,9],MoorInfo[,10])==m)))
  print(paste("spectral features MoorID:",MoorInfo[which(paste(MoorInfo[,9],MoorInfo[,10])==m),10]))
  
  specVar<<-NULL
  whiten2<<-NULL
  specpath<<-NULL

  if(MoorInfo[which(paste(MoorInfo[,9],MoorInfo[,10])==m),7]=="HG_datasets"){
    if(whiten=="y"){
      specpath<<-paste(startcombpath,"/",MoorInfo[which(paste(MoorInfo[,9],MoorInfo[,10])==m),11],"/",Filtype,"p",LMS*100,"x_FO",FO,sep="")
    }else{
      specpath<<-paste(startcombpath,"/",MoorInfo[which(paste(MoorInfo[,9],MoorInfo[,10])==m),11],"/No_whiten",sep="")
    }
    
    if(Decimate=="y"){
      specpath<<-paste(specpath,"_decimate_by_",decimationFactor,sep="")
    }
    
    if(is.null(nrow(specdata))){
      rowcount<<-1
      specVar<<-c(specdata,matrix(1,rowcount,numFeatures))
      specVar<<-rbind(specVar,matrix(1,rowcount,numFeatures+7)) #make 
      noPar<-TRUE
      
    }else{
      specVar<<-specdata[which(specdata[,1]==as.numeric(factor(moors))[which(m==moors)]),]
      rowcount<<-nrow(specVar)
      specVar<<-cbind(specVar,matrix(1,rowcount,numFeatures))
      noPar<-FALSE
    
    }
    
  }else{
    
    if(whiten=="y"){
      whiten2<<-(paste("",Filtype,"p",LMS*100,"x_FO",FO,sep=""))
    }else{
      whiten2<<-"No_whiten"
    }
    if(Decimate=="y"){
      whiten2<<-paste(whiten2,"_decimate_by_",decimationFactor,sep="")
    }
    
    specpath<<-paste(startcombpath,whiten2,sep="")
      
    if(is.null(nrow(specdata))){
      rowcount<<-1
      specVar<<-c(specdata,matrix(1,rowcount,numFeatures))
      specVar<<-rbind(specVar,matrix(1,rowcount,numFeatures+7)) #make
      noPar<-TRUE
      
    }else{
      specVar<<-specdata[which(specdata[,1]==as.numeric(factor(moors))[which(m==moors)]),]
      rowcount<<-nrow(specVar)
      specVar<<-cbind(specVar,matrix(1,rowcount,numFeatures))
      noPar<-FALSE
      
      
    }
    
  }
  
  if(noPar==FALSE){
    print(paste("      for mooring",m))   
    
    if(parallelType=="local"){
    startLocalPar("ImgThresh","MoorInfo","specVar","specpath","rowcount","readWave","freqstat.normalize","lastFeature","std.error","specDo","specgram","imagep","outputpathfiles","jpeg","whichRun")
    }
    
#print("extracting spectral parameters")
specVar2<<-foreach(z=1:rowcount, .packages=c("seewave","tuneR","imager","fpc","cluster")) %dopar% {
  #for(z in 1:rowcount){
  specRow<-unlist(specVar[z,])
  unlist(specDo(z,specRow,specpath))
  #}
}

  if(parallelType=="local"){
  parallel::stopCluster(cluz)
  }

prevdir<-getwd()
setwd(paste(outputpathfiles,"/Image_temp/",sep=""))
unlink('*.jpg')
setwd(prevdir)

specVar<<-do.call('rbind', specVar2)

  }else if(noPar==TRUE){
    z=1
    specRow<-unlist(specVar[1,])
    print(paste("specdo lookup",specRow[1]))
    print(paste("specdo index:",which(specRow[1]==as.numeric(factor(paste(MoorInfo[,9],MoorInfo[,10]))))))
    print(paste("specdo MoorID:",MoorInfo[which(specRow[1]==as.numeric(factor(paste(MoorInfo[,9],MoorInfo[,10])))),10]))
    
    specVar<<-unlist(specDo(z,specRow,specpath))
    

  }
  
specTab<<-rbind(specTab,specVar)
}

  
  return(specTab)
}

specDo<-function(z,featList,specpathh){
  
  #store reused calculations to avoid indexing 
  Start<-featList[2]
  End<-  featList[3]
  if(End-Start<0.1){
    End<-End+(0.9-(End-Start))
    Start<-Start-0.1
    
  }
  Low<-featList[4]
  High<-featList[5]
  
  foo <-readWave(paste(specpathh,"/",MoorInfo[which(featList[1]==as.numeric(factor(paste(MoorInfo[,9],MoorInfo[,10])))),10],".wav",sep=""),Start,End,units="seconds")
  
  fs<-foo@samp.rate
  #foo<-ffilter(foo,from=Low,to=High,output="Wave",wl=512)
  samples<-length(foo@left)
  # demean to remove DC offset
  snd = foo@left - mean(foo@left)
  
  # number of points to use for the fft
  nfft=2024
  
  window<-132
  
  overlap=128
  
  foo.spec <- spec(foo,plot=F, PSD=T,wl=128)
  #foo.spec <- foo.spec[which(foo.spec[,1]<(High/1000)&foo.spec[,1]>(Low/1000)),]#,ylim=c(specVar$Low.Freq..Hz.[z],specVar$High.Freq..Hz.[z])
  foo.specprop <- specprop(foo.spec) #
  foo.meanspec = meanspec(foo, plot=F,ovlp=50,wl=128)#not sure what ovlp parameter does but initially set to 90 #
  #foo.meanspec.db = meanspec(foo, plot=F,ovlp=50,dB="max0",wl=128)#not sure what ovlp parameter does but initially set to 90 #,flim=c(specVar$Low.Freq..Hz.[z]/1000,specVar$High.Freq..Hz.[z]/1000)
  foo.autoc = autoc(foo, plot=F,ovlp=50,wl=64) #
  foo.dfreq = dfreq(foo, plot=F, ovlp=50,wl=64) #tried bandpass argument, limited dfreq to only 2 different values for some reason. Seemed wrong. 
  Startdom<-foo.dfreq[,2][1]
  Enddom<-foo.dfreq[,2][length(foo.dfreq[,2])]
  Mindom <- min(foo.dfreq, na.rm = TRUE)
  Maxdom <- max(foo.dfreq, na.rm = TRUE)
  Dfrange <- Maxdom - Mindom
  featList[8] = rugo(foo@left / max(foo@left)) #rugosity
  featList[9] = crest(foo,wl=128)$C #crest factor
  foo.env = seewave:::env(foo, plot=F) 
  featList[10] = th(foo.env) #temporal entropy
  featList[11] = sh(foo.spec) #shannon entropy
  featList[12] = roughness(foo.meanspec[,2]) #spectrum roughness
  featList[13] = freqstat.normalize(mean(foo.autoc[,2], na.rm=T),Low,High) #autoc mean 
  featList[14] = freqstat.normalize(median(foo.autoc[,2], na.rm=T),Low,High) #autoc.median
  featList[15] = std.error(foo.autoc[,2], na.rm=T) #autoc se
  featList[16] = freqstat.normalize(mean(foo.dfreq[,2], na.rm=T),Low,High) #dfreq mean
  featList[17] = std.error(foo.dfreq[,2], na.rm=T) #dfreq se
  featList[18] = freqstat.normalize(foo.specprop$mean[1],Low,High) #specprop mean
  featList[19] = foo.specprop$sd[1] #specprop sd
  featList[20] = foo.specprop$sem[1] #specprop sem
  featList[21] = freqstat.normalize(foo.specprop$median[1],Low,High) #specprop median
  featList[22] = freqstat.normalize(foo.specprop$mode[1],Low,High) #specprop mode
  featList[23] = foo.specprop$Q25[1] # specprop q25
  featList[24] = foo.specprop$Q75[1] #specprop q75
  featList[25] = foo.specprop$IQR[1] #specprop IQR
  featList[26] = foo.specprop$cent[1] #specrop cent
  featList[27] = foo.specprop$skewness[1] #specprop skewness
  featList[28] = foo.specprop$kurtosis[1] #specprop kurtosis
  featList[29] = foo.specprop$sfm[1] #specprop sfm
  featList[30] = foo.specprop$sh[1] #specprop sh
  featList[31] = foo.specprop$prec[1] #specprop prec
  featList[32] = M(foo,wl=128) #amp env median
  featList[33] = H(foo,wl=128) #total entropy
  #featList[34]<-Q(foo.meanspec.db,plot=F,wl=128)$Q #0s introduced
  #warbler params
  featList[35]<- (sum(sapply(2:length(foo.dfreq[,2]), function(j) abs(foo.dfreq[,2][j] - foo.dfreq[,2][j - 1])))/(Dfrange)) #modinx
  featList[36]<-freqstat.normalize(Startdom,Low,High) #startdom
  featList[37]<-freqstat.normalize(Enddom,Low,High) #enddom 
  featList[38]<-freqstat.normalize(Mindom,Low,High) #mindom
  featList[39]<-freqstat.normalize(Maxdom,Low,High) #maxdom
  featList[40]<-Dfrange #dfrange
  featList[41]<-((Enddom-Startdom)/(End-Start)) #dfslope
  featList[42]  <- freqstat.normalize(lastFeature(fs,foo.meanspec),Low,High)
  
  
  #FEATURES FROM IMAGE 
  
  dir.create(paste(outputpathfiles,"Image_library/",MoorInfo[which(featList[1]==as.numeric(factor(paste(MoorInfo[,9],MoorInfo[,10])))),9],sep=""))
  dir.create(paste(outputpathfiles,"Image_library/",MoorInfo[which(featList[1]==as.numeric(factor(paste(MoorInfo[,9],MoorInfo[,10])))),9],"/",c("No","Yes")[featList[7]+1],sep=""))
  
  im_file_name<-paste(sub("[:alpha:]$", "",unlist(strsplit(sub("(_)(?=[^_]+$)", " ", MoorInfo[which(featList[1]==as.numeric(factor(paste(MoorInfo[,9],MoorInfo[,10])))),1], perl=T), " "))[c(FALSE, TRUE)]),
                      "_B",format(as.POSIXlt(featList[6],origin="1970-01-01",tz="UTC"),"%Y%m%d-%H%M%OS3"),"E",format(as.POSIXlt(featList[6]+End-Start,origin="1970-01-01",tz="UTC"),"%Y%m%d-%H%M%OS3"),"L",as.character(Low),"H",as.character(High),".jpg",sep="")
  if(!file.exists(paste(outputpathfiles,"Image_library/",MoorInfo[which(featList[1]==as.numeric(factor(paste(MoorInfo[,9],MoorInfo[,10])))),9],"/",c("No","Yes")[featList[7]+1],"/",im_file_name,sep=""))){
  
  # create spectrogram
  spec.gram = specgram(x = snd,
                  Fs = fs,
                  window=window,
                  overlap=overlap
  )
  
  # discard phase information
  P = abs(spec.gram$S)
  
  # normalize
  P = P/max(P)
  
  # convert to dB
  P = 15*log10(P)
  
  # config time axis
  t = spec.gram$t
  
  #can use curSpec here to direct to species/detectionID folders. 
  # plot spectrogram
  
  if(whichRun=="GT"){
  
    jpeg(paste(outputpathfiles,"Image_library/",MoorInfo[which(featList[1]==as.numeric(factor(paste(MoorInfo[,9],MoorInfo[,10])))),9],"/",c("No","Yes")[featList[7]+1],"/",im_file_name,sep=""),quality=100)
    imagep(x = t,y = spec.gram$f,z = t(P),col = gray(0:255/255),axes=FALSE,decimate = F,ylim=c(Low,High), drawPalette = FALSE,mar=c(0,0,0,0))
    dev.off()
    
  }else if(whichRun=="NEW"){
    jpeg(paste(outputpathfiles,"Image_temp/",z,".jpg",sep=""),quality=100)
    imagep(x = t,y = spec.gram$f,z = t(P),col = gray(0:255/255),axes=FALSE,decimate = F,ylim=c(Low,High), drawPalette = FALSE,mar=c(0,0,0,0))
    dev.off()
  }
  
  }else{
   #do nothing 
  }
  if(whichRun=="GT"){
  image1<-load.image(paste(outputpathfiles,"Image_library/",MoorInfo[which(featList[1]==as.numeric(factor(paste(MoorInfo[,9],MoorInfo[,10])))),9],"/",c("No","Yes")[featList[7]+1],"/",im_file_name,sep=""))
  }else if(whichRun=="NEW"){
  image1<-load.image(paste(outputpathfiles,"Image_temp/",z,".jpg",sep=""))
  }
  image1<-grayscale(image1, method = "Luma", drop = TRUE)
  f <- ecdf(image1)
  
  image1<-threshold(image1,ImgThresh) 
  image1<-clean(image1,5) %>% imager::fill(1) 
  
  #TEST SECTION
  #jpeg(paste(outputpathfiles,"/Image_temp/Spectrogram",z,".jpg",sep=""),quality=100)
  
   #calculate area chunks x and y 
  chunks<-5
  areaX<- vector("list", length = chunks)
  num<-seq(1,5,1)
  areaX<-lapply(num,function(x) sum(as.matrix(image1)[(((x*480/chunks)-95):(x*480/chunks)),1:480]))
  areaX<-unlist(areaX)
    
  areaY<- vector("list", length = chunks)
  areaY<-lapply(num,function(x) sum(as.matrix(image1)[1:480,(((x*480/chunks)-95):(x*480/chunks))]))
  areaY<-unlist(areaY)
  
  #distinguish islands and calculate area
  labels<-label(image1)
  labels<-as.matrix(labels[1:480,1:480])
  labelsW<-labels #make seperate object for just white regions
  threshBool<-as.matrix(image1[1:480,1:480])
  for(i in 1:length(labelsW)){
    if(threshBool[i]){
      labelsW[i]<-labelsW[i]+1000000
    }
  }
  labelsW[which(labelsW<1000000)]<-0
  
  areaW<-data.frame(table(labelsW))
  areaW<-areaW[which(areaW$labelsW!=0),]
  areaW<-as.vector(areaW[,2])
  
  #calculate convexity of each shape:
  if( length(areaW)<5){
    worthyones<-length(areaW)
  }else{
    worthyones<-5
  }
  IslIndex<-as.numeric(as.character(data.frame(table(labelsW))$labelsW))[-which(as.numeric(as.character(data.frame(table(labelsW))$labelsW))==0)][1:worthyones]
  shapeSlopes<-vector(mode="numeric", length=worthyones)
  shapeCentDistance<-vector(mode="numeric", length=worthyones)
  for(i in 1:worthyones){
    Island<-labelsW
    ind<-IslIndex[i]
    Island[which(Island!=ind)]<-0
    Island<-t(Island) #doesn't plot right but correct orientation for calculating stuff
    #actually want to have upsidedown image so line works... I think 
    #Island<-t(apply(Island, 1, rev)) 
    #par(mar=c(0,0,0,0))
    #image(Island)
    Island_df <- as.matrix(RSAGA::grid.to.xyz(Island))
    Island_df<-Island_df[which(Island_df[,3]==ind),]
    #
    slope<-vector(mode="numeric", length=20)
    b<-vector(mode="numeric", length=20)
    for(a in 1:20){
    Island_df_minx<-NULL
    Island_df_maxx<-NULL
    Island_df_minx<-Island_df[which(Island_df[,1]==min(Island_df[,1])),][sample(1:nrow(Island_df[which(Island_df[,1]==min(Island_df[,1])),]))[1],]
    Island_df_maxx<-Island_df[which(Island_df[,1]==max(Island_df[,1])),][sample(1:nrow(Island_df[which(Island_df[,1]==max(Island_df[,1])),]))[1],]
    slope[a]<-(Island_df_maxx[2]-Island_df_minx[2])/(Island_df_maxx[1]-Island_df_minx[1])
    b[a]<-Island_df_minx[2]-(slope*Island_df_minx[1])
    }
    slope<-mean(slope)
    b<-mean(b)
    
    shapeSlopes[i]<-atan(slope)

    positionsX <- apply(Island[1:480,1:480], 1, function(x) which(x==ind))
    positionsY <- apply(Island[1:480,1:480], 2, function(x) which(x==ind))
    cX<-mean(unlist(positionsX,recursive = TRUE),na.rm=TRUE)#xavg
    cY<-mean(unlist(positionsY,recursive = TRUE),na.rm=TRUE)#yavg
    
    cent<-c(cX,cY)
    linep1<-c(b,0)
    linep2<-c(Island_df_minx[1],Island_df_minx[2])
    
    v1 <- linep1 - linep2
    v2 <- cent - linep1
    m <- cbind(v1,v2)
    d <- det(m)/sqrt(sum(v1*v1))
    
    shapeCentDistance[i]<-d
    
  }
  
  perConcave<-length(shapeCentDistance[which(shapeCentDistance>0)])/length(shapeCentDistance)
  
  #TEST SECTION
  #labelsW[which(labelsW!=1000001)]<-0
  #plot(raster(t(labelsTest)))
  #plot(raster(t(labelsW)))
  #plot(smooth(rasterToPolygons(raster(t(labelsW)))))

  #hough lines
  test9<-hough_line(image1,data.frame = TRUE)
  test9<- cbind(test9,(-(cos(test9$theta)/sin(test9$theta))),test9$rho/sin(test9$theta))
  test9[which(is.infinite(test9[,4])&test9[,4]<0),4]<-min(c(-2.76824e+18,min(test9[-which(is.infinite(test9[,4])),4])))
  test9[which(is.infinite(test9[,4])&test9[,4]>0),4]<-max(c(2.76824e+18,max(test9[-which(is.infinite(test9[,4])),4])))
  test9[which(is.infinite(test9[,5])&test9[,5]<0),5]<-min(c(-2.76824e+18,min(test9[-which(is.infinite(test9[,5])),5])))
  test9[which(is.infinite(test9[,5])&test9[,5]>0),5]<-max(c(2.76824e+18,max(test9[-which(is.infinite(test9[,5])),5])))
  
  Bestline<-test9[which.max(test9$score),]
  Bestlines<-test9[which(test9$score>=max(test9$score)*.7),]
  
  #calculate centroid of area (could also add centroid of largest area)
  positionsX <- apply(image1[1:480,1:480], 1, function(x) which(x==TRUE))
  positionsY <- apply(image1[1:480,1:480], 2, function(x) which(x==TRUE))
  
  #calculate stats on horizontal and vertical switches
  hpEven<-c()
  hpEven<-c(hpEven,sum(diff(image1[1:480,24]) == 1) + sum(diff(image1[1:480,24]) == -1))
  hpEven<-c(hpEven,sum(diff(image1[1:480,96]) == 1) + sum(diff(image1[1:480,96]) == -1))
  hpEven<-c(hpEven,sum(diff(image1[1:480,192]) == 1) + sum(diff(image1[1:480,192]) == -1))
  hpEven<-c(hpEven,sum(diff(image1[1:480,288]) == 1) + sum(diff(image1[1:480,288]) == -1))
  hpEven<-c(hpEven,sum(diff(image1[1:480,384]) == 1) + sum(diff(image1[1:480,384]) == -1))
  hpEven<-c(hpEven,sum(diff(image1[1:480,456]) == 1) + sum(diff(image1[1:480,456]) == -1))
  vpEven<-c()
  vpEven<-c(vpEven,sum(diff(image1[24,1:480]) == 1) + sum(diff(image1[24,1:480]) == -1))
  vpEven<-c(vpEven,sum(diff(image1[96,1:480]) == 1) + sum(diff(image1[96,1:480]) == -1))
  vpEven<-c(vpEven,sum(diff(image1[192,1:480]) == 1) + sum(diff(image1[192,1:480]) == -1))
  vpEven<-c(vpEven,sum(diff(image1[288,1:480]) == 1) + sum(diff(image1[288,1:480]) == -1))
  vpEven<-c(vpEven,sum(diff(image1[384,1:480]) == 1) + sum(diff(image1[384,1:480]) == -1))
  vpEven<-c(vpEven,sum(diff(image1[456,1:480]) == 1) + sum(diff(image1[456,1:480]) == -1))
  
  #add new variables
  featList[43]<-which.max(areaX) #areaXmaxP
  featList[44]<-max(areaX) #areaXmax
  featList[45]<-max(areaX)/sum(areaX) #areaXdom
  featList[46]<-std.error(areaX) #areaXstd
  
  featList[47]<-which.max(areaY) #areaYmaxP
  featList[48]<-max(areaY) #areaYmax
  featList[49]<-max(areaY)/sum(areaY)#areaYdom
  featList[50]<-std.error(areaY)#areaYstd
  
  #featList[51]<-std.error(areaW) #Areaspread
  featList[52]<-max(areaW)#AreaTop
  featList[53]<-max(areaW)/(sum(areaW))#AreaTopDom
  featList[54]<-if(length(areaW)>=3){sum(-sort(-areaW)[1:3])/sum(areaW)}else{1}#AreaTop3Dom
  featList[55]<-length(areaW)#NumShapes

  featList[56]<-Bestline[2]#bestSlopeHough radian
  featList[57]<-Bestline[3]#bestBHough radian
  featList[58]<-Bestline[4]#bestSlopeHough
  featList[59]<-Bestline[5]#bestBHough
  featList[61]<-median(Bestlines[,2])#medSlope radian
  featList[62]<-median(Bestlines[,3])#medB radian 
  featList[63]<-median(Bestlines[,4])#medSlope
  featList[64]<-median(Bestlines[,5])#medB
  featList[65]<-nrow(Bestlines)#numGoodlines
  

  featList[66]<-mean(unlist(positionsX,recursive = TRUE),na.rm=TRUE)#xavg
  featList[67]<-mean(unlist(positionsY,recursive = TRUE),na.rm=TRUE)#yavg
  
  featList[68]<-mean(hpEven)#switchesX
  featList[69]<-std.error(hpEven)#switchesXreg
  featList[70]<-max(hpEven)#switchesXmax
  featList[71]<-min(hpEven)#switchesXmin
  
  featList[72]<-mean(vpEven)#switchesY
  featList[73]<-std.error(vpEven)#switchesYreg
  featList[74]<-max(vpEven)#switchesYmax
  featList[75]<-min(vpEven)#switchesYmin
  
  #individual shape features:
  featList[76]<-mean(shapeSlopes) #avg slope. using arctan(slope) hoping it makes numeric comparisons better since slope is nonlinear especially at large values
  featList[77]<-std.error(shapeSlopes) #var slope
  featList[78]<-sum(shapeCentDistance) #sum of all centroid distances (hopefully weighted towards larger shapes that have better potential for concavity)
  featList[79]<-mean(shapeCentDistance) #other way of comparing centroid distances
  featList[80]<-perConcave #% positive distances.  
  
  return(featList)
}

dataConflicts<-function(cData){
  
#remove 0s that conflict with 1s to make sure model is not fed good calls as examples from other category. 
cData$RTb<-cData$RTFb+cData$FileOffsetBegin
cData$RTe<-cData$RTFb+cData$FileOffsetEnd
cData$remove<-0

#remove conflicts so don't feed model that positives of a sound are + for a negative category
if(any(cData$detectionType==1)){

for(h in 1:nrow(cData)){
  if(cData$detectionType[h]==1){
    #do nothing
  }else{
    gvec <- cData[which(cData$detectionType==1&cData$File==cData$File[h]&cData$Species!=cData$Species[h]),]
    if(nrow(gvec)>0){
      for(g in 1:nrow(gvec)){
        if((gvec$RTb[g]>=cData$RTb[h]&gvec$RTb[g]<cData$RTe[h])|(gvec$RTe[g]<=cData$RTe[h]&gvec$RTe[g]>cData$RTb[h])|(gvec$RTe[g]<=cData$RTe[h]&gvec$RTb[g]>=cData$RTb[h])|(gvec$RTe[g]>=cData$RTe[h]&gvec$RTb[g]<=cData$RTb[h])){
          cData$remove[h]<-1
        }
        break
      }
    }
  }
}

cData<-cData[which(cData$remove==0),]

cData$remove<-NULL

}

#combine boxes with overlap and true detection (0.06 corresponding to lowest timediffself parameter GSs- works at whatever value though)
cData$combine<-0
p=1

if(any(cData$detectionType==1)){
for(h in 1:length(unique(cData$MooringID))){
  
  cDataPos<-cData[which(cData$detectionType==1&cData$MooringID==unique(cData$MooringID)[h]),]
  cDataPos<-cDataPos[order(cDataPos$`Begin Time (s)`),]
  
  #merge all overlapping 1s. 
  RollEnd<-cDataPos$`End Time (s)`[1]
  for(i in 1:(nrow(cDataPos)-1)){
    if(cDataPos$`Begin Time (s)`[i+1]<RollEnd){
      cDataPos$combine[i]<-p
      cDataPos$combine[i+1]<-p
      if(RollEnd<cDataPos$`End Time (s)`[i+1]){
        RollEnd<-cDataPos$`End Time (s)`[i+1]
      }
    }else{
      p=p+1
      RollEnd<-cDataPos$`End Time (s)`[i+1]
    }
  }
  
  cData<-cData[-which(cData$detectionType==1&cData$MooringID==unique(cData$MooringID)[h]),]
  cData<-rbind(cData,cDataPos)
}
  #combine boxes with overlap and no detection
  p=max(cData$combine)
}

cData$RTb<-NULL
cData$RTe<-NULL
cData$remove<-NULL




for(h in 1:length(unique(cData$MooringID))){
  
  cDataPos<-cData[which(cData$detectionType==0&cData$MooringID==unique(cData$MooringID)[h]),]
  cDataPos<-cDataPos[order(cDataPos$`Begin Time (s)`),]
  
  #merge all overlapping 1s. 
  RollEnd<-cDataPos$`End Time (s)`[1]
  for(i in 1:(nrow(cDataPos)-1)){
    if(cDataPos$`Begin Time (s)`[i+1]<RollEnd){
      cDataPos$combine[i]<-p
      cDataPos$combine[i+1]<-p
      if(RollEnd<cDataPos$`End Time (s)`[i+1]){
        RollEnd<-cDataPos$`End Time (s)`[i+1]
      }
    }else{
      p=p+1
      RollEnd<-cDataPos$`End Time (s)`[i+1]
    }
  }
  
  cData<-cData[-which(cData$detectionType==0&cData$MooringID==unique(cData$MooringID)[h]),]
  cData<-rbind(cData,cDataPos)
  
}

#combine boxes 
cData$meantime<-NULL

if(any(cData$combine>0)){

for(g in unique(cData$combine)[which(unique(cData$combine)!=0)]){
  
  cDataGroup<-cData[which(cData$combine==g),]
  
  cDataRow<-data.frame(cDataGroup[1,])
  
  cDataRow[,4]<-min(cDataGroup$`Begin Time (s)`)
  cDataRow[,5]<-max(cDataGroup$`End Time (s)`)
  cDataRow[,6]<-min(cDataGroup$`Low Freq (Hz)`)
  cDataRow[,7]<-max(cDataGroup$`High Freq (Hz)`)
  cDataRow[,13]<-min(cDataGroup$FileOffsetBegin)
  cDataRow[,14]<-max(cDataGroup$FileOffsetEnd)
  cDataRow[,18]<-paste(unique(cDataGroup$Species),collapse=",")
  cDataRow[,19]<-cDataGroup$detectionType[1]#internally consistent
  cDataRow[,20]<-0

  cData<-cData[-which(cData$combine==g),]
  names(cDataRow)<-names(cData)
  cData<-rbind(cData,cDataRow)
}
  
}

cData$combine<-NULL
#add frequency stats to cData 
cData$meanfreq<- (cData$`Low Freq (Hz)`+cData$`High Freq (Hz)`)/2
cData$freqrange<- (cData$`High Freq (Hz)`-cData$`Low Freq (Hz)`)
cData$meantime<- (cData$`Begin Time (s)`+cData$`End Time (s)`)/2

cData$Selection<-seq(1,nrow(cData))

return(cData)
}

process_data<-function(){
  
  allRSTF<-NULL
#Combine and configure spread detectors. 

  #seperate by species
resltsTabTemp<-resltsTab[which(resltsTab$Species==s),]

if(useMasterGT=="y" & addToMaster=="y"){
  resltsTabTemp<-resltsTabTemp[-which(paste(resltsTabTemp$Species,resltsTabTemp$MooringID) %in% moorsCompleted),]
}

print(paste("Combining spread for"))
  
for(e in unique(resltsTabTemp$sound.files)){
  resltsVar<-resltsTabTemp[which(resltsTabTemp$sound.files==e),]
  print(paste("    ",e))
  for(f in 1:length(unique(resltsVar$bottom.freq))){
    resltsVar[resltsVar$bottom.freq==sort(unique(resltsVar$bottom.freq))[f],14]<-f
  }
  colnames(resltsVar)[14] <- "detectorRank"
  resltsVar$detectorRank<-as.numeric(resltsVar$detectorRank)
  resltsVar$group[1]<-1
      
  #make sure these values are what they say they are- encountered some behaviors to suggest there was some behind the scenes decimals earlier. This may not be necessary. 
  resltsVar$start<-as.numeric(as.character(resltsVar$start))
  resltsVar$end<-as.numeric(as.character(resltsVar$end))
  resltsVar$meantime<-(resltsVar$start+resltsVar$end)/2

    if(s=="RW"){
      #need to order chronologically. 
      resltsVar<-resltsVar[order(resltsVar$sound.files,resltsVar$meantime,resltsVar$bottom.freq),]
      
      #assign groups based on groupInt value
      f<-1
      
      #index columns to process faster:
      gTime<-resltsVar[,16]
      gGroup<-resltsVar[,15]
      
    print("assigning group values")
    for(z in 1:(nrow(resltsVar)-1)){
      if(gTime[z]+groupInt>=gTime[z+1]){
        gGroup[z+1]<-f
      }else{
        f<-f+1
        gGroup[z+1]<-f
      }
    }
    }else if(s=="GS"|s=="BP"){
      
      #need to order chronologically. 
      resltsVar<-resltsVar[order(resltsVar$sound.files,resltsVar$start,resltsVar$top.freq),]
      
      #assign groups based on groupInt value
      f<-1
      
      #index columns to process faster:
      gTimeS<-resltsVar$start        
      gGroup<-resltsVar[,15]
      nexstart<-gTimeS[1]
      print("assigning group values")
      for(z in 1:(nrow(resltsVar)-1)){
        if(gTimeS[z]+groupInt>=gTimeS[z+1]){
          gGroup[z+1]<-f
        }else{
          f<-f+1
          gGroup[z+1]<-f
        }
      }
    }
    
    resltsVar[,15]<-gGroup
    
    #remove groups based on grpsize value
    removegrp <- table(resltsVar$group)
    resltsVar <- subset(resltsVar, group %in% names(removegrp[removegrp > (grpsize-1)]))
    Matdata<<-data.matrix(resltsVar[,c(14,15,16,4,5)])
  
    #updated algorithm, optimized for performance. avoids r bind
    print(paste("calculating best runs for each group"))
    #print(system.time(parAlgo(Matdata)))
    wantedSelections<-parAlgo(Matdata)
    
    Matdata<<-NULL
  print("start process wanted selections")
    if(s=="RW"){
    resltsVar<-resltsVar[which(as.integer(rownames(resltsVar)) %in% as.integer(wantedSelections)),]
    }else if(s=="GS"|s=="BP"){
    wantedSelections<-t(wantedSelections)
    c=1
    
    for(u in 1:nrow(wantedSelections)-1){
      wantedSelections[u,2]<-c
      if(wantedSelections[u+1,2]==999999999){
      c=c+1
    }
    }
    wantedSelections<-wantedSelections[which(wantedSelections[,1]!=999999999),]
    wantedSelections<-wantedSelections[which(!duplicated(wantedSelections[,1])),] #throw out duplicates, later assign boxes just using first position of highest in group and last end time in group
    resltsVar<-resltsVar[which(as.integer(rownames(resltsVar)) %in% wantedSelections[,1]),]
    print("start create new groups values")
    #create new groups values
    
    resltsVar<<-resltsVar
    wantedSelections<<-wantedSelections

    if(parallelType=="local"){
      startLocalPar("wantedSelections","resltsVar")
    }
    
    vec<-foreach(u=1:nrow(resltsVar),.combine="c") %dopar% {
      as.numeric(wantedSelections[which(wantedSelections[,1]==as.integer(rownames(resltsVar[u,]))),2] )
    }

    resltsVar[,15]<-vec
    
    if(parallelType=="local"){
      parallel::stopCluster(cluz)
    }
    
    }
    
    if(nrow(resltsVar)==0){
      write.table("FINAL There were no detections",paste(outputpath,runname,"/",e,"/FINAL_Summary_spread_",substr(resltsVar$detector[1],1,3),"_",length(detectorsspr),".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE,col.names=FALSE)
    }else{
      print("start create new table")
      colClasses = c("numeric", "character","numeric","numeric", "numeric","numeric","numeric","character","character","character")
      resltsTabFinal <- read.csv(text="Selection,View,Channel,Begin Time (s),End Time (s),Low Freq (Hz),High Freq (Hz), MooringID, MooringName, MooringCode, sound.files", colClasses = colClasses)
      colnames(resltsTabFinal)<-c("Selection","View","Channel","Begin Time (s)","End Time (s)","Low Freq (Hz)","High Freq (Hz)","MooringID","MooringName","MooringCode","sound.files")
      
      if(s=="RW"){
      p=1
      for(j in unique(resltsVar$group)){
        grpminfreq <- min(resltsVar[resltsVar$group==j,6])
        grpmaxfreq <- max(resltsVar[resltsVar$group==j,7])
        grpstarttime <- min(resltsVar[resltsVar$group==j,4])
        grpendtime <- max(resltsVar[resltsVar$group==j,5])
        
        resltsTabFinal[p,1]<-p
        resltsTabFinal[p,4]<-grpstarttime
        resltsTabFinal[p,5]<-grpendtime
        resltsTabFinal[p,6]<-grpminfreq
        resltsTabFinal[p,7]<-grpmaxfreq
        resltsTabFinal[p,8]<-as.character(resltsVar$MooringID[1])
        resltsTabFinal[p,9]<-as.character(resltsVar$MooringName[1])
        resltsTabFinal[p,10]<-sub("[:alpha:]$", "",unlist(strsplit(sub("(_)(?=[^_]+$)", " ", resltsVar$MooringName[1], perl=T), " "))[c(FALSE, TRUE)])
        resltsTabFinal[p,11]<-as.character(resltsVar$sound.files[1])
        
        p<-p+1
      }
        
      }else if(s=="GS"|s=="BP"){
        p=1
        for(j in unique(resltsVar$group)){
          grp<-resltsVar[resltsVar$group==j,]
          grpminfreq <- min(grp[,6])
          grpmaxfreq <- max(grp[,7])
          if(any(duplicated(grp[,4]))){             #min(grp[which.max(grp[,13]),4]) #for start of 1st detector in run
          grpstarttime <- getmode(grp[,4])
          }else{
          grpstarttime <- median(grp[,4])
          }
          grpendtime <- max(grp[,5])
          
          resltsTabFinal[p,1]<-p
          resltsTabFinal[p,4]<-grpstarttime
          resltsTabFinal[p,5]<-grpendtime
          resltsTabFinal[p,6]<-grpminfreq
          resltsTabFinal[p,7]<-grpmaxfreq
          resltsTabFinal[p,8]<-as.character(resltsVar$MooringID[1])
          resltsTabFinal[p,9]<-as.character(resltsVar$MooringName[1])
          resltsTabFinal[p,10]<-sub("[:alpha:]$", "",unlist(strsplit(sub("(_)(?=[^_]+$)", " ", resltsVar$MooringName[1], perl=T), " "))[c(FALSE, TRUE)])
          resltsTabFinal[p,11]<-as.character(resltsVar$sound.files[1])
          
          p<-p+1
        }
        
        if(s=="GS"){
          #remove fragments (100hz or under) that are resonably high up 
        resltsTabFinal<-resltsTabFinal[which(!(((resltsTabFinal[,7]-resltsTabFinal[,6])<=100)&(resltsTabFinal[,6]>=225))),] #&resltsTabFinal[,6]>=225)
        
        #remove fragments (100hz or under) also on the low end
        resltsTabFinal<-resltsTabFinal[which(!(((resltsTabFinal[,7]-resltsTabFinal[,6])<=100)&(resltsTabFinal[,7]<150)&(resltsTabFinal[,5]-resltsTabFinal[,4]<=0.5))),] #&resltsTabFinal[,6]>=225)
        
        }
        
        #change end time of call to start of next call if they overlap. 
        resltsTabFinal<-resltsTabFinal[order(resltsTabFinal[,4]),]
        for(b in 1:(nrow(resltsTabFinal)-1)){
          if(resltsTabFinal[b,5]>resltsTabFinal[b+1,4])
            resltsTabFinal[b,5]<-resltsTabFinal[b+1,4]
        }
        
      }
        
        #remove ones with no data 
        resltsTabFinal<-resltsTabFinal[which((resltsTabFinal[,5]-resltsTabFinal[,4])!=0),]
        
        resltsTabFinal$View<-"Spectrogram 1"
        resltsTabFinal$Channel<-1
        
        resltsTabFinal<-resltsTabFinal[order(resltsTabFinal$`Begin Time (s)`),]
        
        #remove detections that do not fit min/max duration parameters  #
        resltsTabFinal$remove<-0
        for(f in 1:nrow(resltsTabFinal)){
          if((resltsTabFinal[f,5]-resltsTabFinal[f,4])>Maxdur|(resltsTabFinal[f,5]-resltsTabFinal[f,4])<Mindur){
            resltsTabFinal$remove[f]<-1
          }
        }
        
        resltsTabFinal<-resltsTabFinal[which(resltsTabFinal$remove==0),]
        
        resltsTabFinal$remove<-NULL
        
        #######
        
        
        ##############do the durtab stuff
        if(whiten=="y"){
          whiten2<-(paste("",Filtype,"p",LMS*100,"x_FO",FO,sep=""))
        }else{
          whiten2<-"No_whiten"
        }
        if(Decimate=="y"){
          whiten2<-paste(whiten2,"_decimate_by_",decimationFactor,sep="")
        }
        
        resltsTabFinal$Species<-s #needed to add this here to assess unique data 

        curSpec<-MoorInfo[which(paste(MoorInfo[,9],MoorInfo[,10])==paste(resltsTabFinal$Species[1],resltsTabFinal[1,8])),11]
        print(curSpec)
        print(paste(resltsTabFinal[1,12],resltsTabFinal[1,8]))
        print(paste(MoorInfo[,9],MoorInfo[,10]))
        if(MoorInfo[which(paste(MoorInfo[,9],MoorInfo[,10])==paste(resltsTabFinal[1,12],resltsTabFinal[1,8])),7]=="HG_datasets"){
          pathh<-paste(startcombpath,curSpec,sep="")
        }else if(MoorInfo[which(paste(MoorInfo[,9],MoorInfo[,10])==paste(resltsTabFinal[1,12],resltsTabFinal[1,8])),7]=="Full_datasets"){
          pathh<-startcombpath   
        }
        
        resltsTabFinal$Species<-NULL #add this back later to not mess up future indexes that rely on position. 
        
        filePath<-paste(pathh,whiten2,sep="/")
        
        durTab<-read.csv(paste(filePath,"/",resltsTabFinal[1,8],"_SFiles_and_durations.csv",sep=""))
        
        #add information on original sound files and calculate time since file start
        resltsTabFinal$File<-""
        FileStartSec<-c()
        resltsTabFinal$FileOffsetBegin<-0
        resltsTabFinal$FileOffsetEnd<-0
        
        print("calculate file ID and begin time and end time relative to file for sound.file")
          for(c in 1:nrow(resltsTabFinal)){
            index<-findInterval(resltsTabFinal[c,4], c(0,durTab$CumDur))
            resltsTabFinal$File[c]<-as.character(durTab[index,2])
            durations<-durTab[index,3]
            FileStartSec[c]<-durTab[index,4]-durTab[index,3]
          }
        resltsTabFinal$FileOffsetBegin<-resltsTabFinal$`Begin Time (s)`-FileStartSec
        resltsTabFinal$FileOffsetEnd<-resltsTabFinal$`End Time (s)`-FileStartSec  
        
        #add real times! 
        
        resltsTabFinal$Selection<-seq(1,nrow(resltsTabFinal))
        
        RT<-substrRight(as.character(resltsTabFinal$File),19)
        RT<-substr(RT,1,15)
        RT<-strptime(RT,format='%Y%m%d_%H%M%S',tz="UTC")
        
        resltsTabFinal$RTfile<-RT
        resltsTabFinal$RTFb<-RT
        resltsTabFinal$RTFe<-RT+durations
        
        ###
        resltsTabFinal$Species<-s

        #for HG,
          
        
        allRSTF<- rbind(allRSTF,resltsTabFinal)
        
        write.table(resltsTabFinal,paste(outputpath,runname,"/",resltsTabFinal$MooringID[1],"FINAL_Summary_spread_",substr(resltsTabFinal$detector[1],1,3),"_",length(detectorsspr),".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
        
      }
  }

############
return(allRSTF)
}

combineDecRaven<-function(){
  
  resltsTab<-NULL
  
  #decimate and whiten are not 
  if(whiten=="y"){
    whiten2<-(paste("",Filtype,"p",LMS*100,"x_FO",FO,sep=""))
  }else{
    whiten2<-"No_whiten"
  }
  
  if(Decimate=="y"){
    whiten2<-paste(whiten2,"_decimate_by_",decimationFactor,sep="")
  }
  
  #run sound files:
  decNeeded<-"n"
  for(m in 1:nrow(MoorInfo)){
    
    curSpec<-paste(MoorInfo[m,11])
    
    durTab<-NULL
    durTab2<-NULL
    
    if(whiten=="y"){
      if(MoorInfo[m,7]=="HG_datasets"){
        if(file.exists(paste(startcombpath,curSpec,"/",whiten2,"/",MoorInfo[m,10],"_SFiles_and_durations.csv",sep=""))){
        }else{
          stop("First run without whitening, then use whitening filter in Raven and path to folder accordingly")
        }
      }else if(MoorInfo[m,7]=="Full_datasets"){
        if(file.exists(paste(startcombpath,"/",whiten2,"/",MoorInfo[m,10],"_SFiles_and_durations.csv",sep=""))){
        }else{
          stop("First run without whitening, then use whitening filter in Raven and path to folder accordingly")
        }
      }
    }
    
    if(whiten=="n"){
      if(MoorInfo[m,7]=="HG_datasets"){
        if(Decimate=="y"&file.exists(paste(drivepath,MoorInfo[m,7],MoorInfo[m,1],paste(curSpec,"_yesUnion",sep=""),paste(MoorInfo[m,10],"_decimate_by_",decimationFactor,sep=""),sep="/"))){
          sfpath<-paste(drivepath,MoorInfo[m,7],MoorInfo[m,1],paste(curSpec,"_yesUnion",sep=""),paste(MoorInfo[m,10],"_decimate_by_",decimationFactor,sep=""),sep="/")
          decNeeded<-"n"
        }else if(Decimate=="y"&!file.exists(paste(drivepath,MoorInfo[m,7],MoorInfo[m,1],paste(curSpec,"_yesUnion",sep=""),paste(MoorInfo[m,10],"_decimate_by_",decimationFactor,sep=""),sep="/"))){
          sfpath<-paste(drivepath,MoorInfo[m,7],MoorInfo[m,1],paste(curSpec,"_yesUnion",sep=""),sep="/")
          decNeeded<-"y"
        }else if(Decimate=="n"){
          sfpath<-paste(drivepath,MoorInfo[m,7],MoorInfo[m,1],paste(curSpec,"_yesUnion",sep=""),sep="/")
          decNeeded<-"n"
        }
      }else if(MoorInfo[m,7]=="Full_datasets"){
        if(Decimate=="y"&file.exists(paste(drivepath,MoorInfo[m,7],MoorInfo[m,1],paste(MoorInfo[m,10],"_decimate_by_",decimationFactor,sep = ""),sep="/"))){
          sfpath<-paste(drivepath,MoorInfo[m,7],MoorInfo[m,1],paste(MoorInfo[m,10],"_decimate_by_",decimationFactor,sep = ""),sep="/")
          decNeeded<-"n"
        }else if(Decimate=="y"&!file.exists(paste(drivepath,MoorInfo[m,7],MoorInfo[m,1],paste(MoorInfo[m,10],"_decimate_by_",decimationFactor,sep = ""),sep="/"))){
          sfpath<-paste(drivepath,MoorInfo[m,7],"/",MoorInfo[m,1],sep = "")
          decNeeded<-"y"
        }else if(Decimate=="n"){
          sfpath<-paste(drivepath,MoorInfo[m,7],MoorInfo[m,1],sep = "/")
          decNeeded<-"n"
        }
      }
      
      #maybe put whiten check right here? 
      if(MoorInfo[m,8]=="open"){
        sound_files <- dir(sfpath,pattern=".wav")[which(dir(sfpath,pattern=".wav")==MoorInfo[m,4]):which(dir(sfpath,pattern=".wav")==MoorInfo[m,5])]
        sound_filesfullpath<-paste(sfpath,sound_files,sep = "")
      }else if(MoorInfo[m,8]=="month"&decNeeded=="y"){ #assuming month wil only be used when acthing on full moorings. If a different use case conditionals may not work. 
        sound_files<-c()
        sound_filesfullpath<-
        allMonths<-dir(sfpath)
        for(i in allMonths){
          sound_files<-c(soundfiles,dir(paste(sfpath,i,sep="/"),pattern=".wav"))
          sound_filesfullpath<-c(sound_filesfullpath,paste(sfpath,i,dir(paste(sfpath,i,sep="/"),pattern=".wav"),sep="/"))
        }
        sound_files <- NULL #need to look at how mooring is structured but should work fine for sox with a list of full path files. 
        sound_filesfullpath<-NULL
      }
      
      if(decNeeded=="y"){
        oldPath<-sfpath
        if(MoorInfo[m,7]=="HG_datasets"){
          sfpath<-paste(drivepath,MoorInfo[m,7],MoorInfo[m,1],paste(curSpec,"_yesUnion",sep=""),paste(MoorInfo[m,10],"_decimate_by_",decimationFactor,sep=""),sep="/")
          dir.create(sfpath)
        }else if(MoorInfo[m,7]=="Full_datasets"){
          sfpath<-paste(drivepath,MoorInfo[m,7],MoorInfo[m,1],paste(MoorInfo[m,10],"_decimate_by_",decimationFactor,sep = ""),sep="/")
          dir.create(sfpath)
        }
        decimateData(m,oldPath,sfpath)
      }
      
      if(length(sound_files)>=fileCombinesize){
        fileSizeInt<-fileCombinesize
        iterate_SF<-c(1,2)
        fileSizeInt2<-fileCombinesize2ndIt 
      }else{
        iterate_SF<-1
      }
      
      bigFile_breaks<-c(seq(1,length(sound_files),fileCombinesize),length(sound_files)) 
      
      for(b in 1:(length(bigFile_breaks)-1)){
        sound_filesB <- dir(sfpath)[bigFile_breaks[b]:(bigFile_breaks[b+1]-1)]
        if(b==length(bigFile_breaks)-1){
          sound_filesB <- dir(sfpath)[bigFile_breaks[b]:length(sound_files)]
        }
        sound_filesfullpathB <- paste(sfpath,"/",sound_filesB,sep = "")
        if(MoorInfo[m,7]=="HG_datasets"){
          pathh<-paste(startcombpath,curSpec,sep="")
        }else if(MoorInfo[m,7]=="Full_datasets"){
          pathh<-startcombpath   
        }
        if(length(iterate_SF)>1){
          filePathNoTemp<-pathh
          pathh<-paste(pathh,"/temp",sep="")
          pad<-sprintf("%02d",b)
        }else{
          filePathNoTemp<-pathh
          pad<-""
        }
        filePath<-paste(pathh,whiten2,sep="/")
        filePathNoTemp<-paste(filePathNoTemp,whiten2,sep="/")#look in final folder to see if SFiles is populated. 
        dir.create(filePath)
        combSound<-paste(filePath,"/",pad,MoorInfo[m,10],".wav",sep="")
        if(file.exists(paste(filePathNoTemp,"/",MoorInfo[m,10],"_SFiles_and_durations.csv",sep=""))){
          durTab <-read.csv(paste(filePathNoTemp,"/",MoorInfo[m,10],"_SFiles_and_durations.csv",sep=""))  
          stopRun<-TRUE
        }else if(!file.exists(paste(filePathNoTemp,"/",MoorInfo[m,10],"_SFiles_and_durations.csv",sep=""))){
          durTab<-sox.write(1,m,b,pathh,sound_filesfullpathB,combSound,durTab,pad,bigFile_breaks,sound_filesB,durTab2)
          stopRun<-FALSE
        }
        
      }
      
      did2<-FALSE
      
      if(length(iterate_SF)>1&stopRun==TRUE){
        sfpath<-filePathNoTemp
        
        sound_files <- dir(sfpath,pattern=".wav")[grep(MoorInfo[m,10],dir(sfpath,pattern=".wav"))]
        sound_filesfullpath<-paste(sfpath,sound_files,sep = "")
        bigFile_breaks<-c(seq(1,length(sound_files),fileSizeInt2),length(sound_files)) 
        did2<-TRUE
        
      }
      #go through again if more than certain # of files 
      if(length(iterate_SF)>1&stopRun!=TRUE){
        sfpath<-filePath
        
        sound_files <- dir(sfpath,pattern=".wav")
        sound_filesfullpath<-paste(sfpath,sound_files,sep = "")
        
        bigFile_breaks<-c(seq(1,length(sound_files),fileSizeInt2),length(sound_files)) 
        
        for(b in 1:(length(bigFile_breaks)-1)){
          sound_filesB <- dir(sfpath)[bigFile_breaks[b]:(bigFile_breaks[b+1]-1)]
          if(b==length(bigFile_breaks)-1){
            sound_filesB <- dir(sfpath)[bigFile_breaks[b]:length(sound_files)]
          }
          sound_filesfullpathB <- paste(sfpath,"/",sound_filesB,sep = "")
          if(MoorInfo[m,7]=="HG_datasets"){
            pathh<-paste(startcombpath,curSpec,sep="")
          }else if(MoorInfo[m,7]=="Full_datasets"){
            pathh<-startcombpath   
          }
          pathh<-paste(pathh,sep="")
          pad<-sprintf("%02d",b)
          
          filePath<-paste(pathh,whiten2,sep="/")
          dir.create(filePath)
          if(length(bigFile_breaks)<=2){
            combSound<-paste(filePath,"/",MoorInfo[m,10],".wav",sep="")
          }else if(length(bigFile_breaks)>2){
            combSound<-paste(filePath,"/",pad,MoorInfo[m,10],".wav",sep="")
          }
          durTab2<-sox.write(2,m,b,pathh,sound_filesfullpathB,combSound,durTab,pad,bigFile_breaks,sound_filesB,durTab2)
        }
        
        unlink(paste(pathh,"/temp/",whiten2,sep=""),recursive=TRUE)
        durTab<-durTab2
        
        did2<-TRUE
        
      }
      
      if(did2==FALSE&length(bigFile_breaks)>2){
        bigFile_breaks<-c(1,fileCombinesize) #if data is able to be combined in 1 file and ran again it will not update bF_b so this is a fix 
      }
      
      if(onlyPopulate=="n"){
        for(b in 1:(length(bigFile_breaks)-1)){
          if(length(bigFile_breaks)>2&did2==TRUE){
            combname<- paste(sprintf("%02d",b),MoorInfo[m,10],"_files",bigFile_breaks[b],".wav",sep="")
          }else{
            combname<- paste(MoorInfo[m,10],".wav",sep="")
          }
          #run detector(s)
          resltsTab<-runRavenDetector(m,filePathNoTemp,combname,resltsTab)
        }
      }
      #write durtab to file
      write.csv(durTab,paste(filePath,"/",MoorInfo[m,10],"_SFiles_and_durations.csv",sep=""),row.names = F)
      
    }else if(whiten=="y"){
      if(onlyPopulate=="y"){
        stop("you cannot populate and whiten. First populate, then whiten in raven, then use whiten argument to specify you have whitened data")
      }
      if(MoorInfo[m,7]=="HG_datasets"){
        filePath<-paste(startcombpath,curSpec,"/",whiten2,sep="")
      }else if(MoorInfo[m,7]=="Full_datasets"){
        filePath<-paste(startcombpath,"/",whiten2,sep="")   
      }  
      for(i in intersect(list.files(filePathNoTemp,pattern = paste(MoorInfo[m,10])), list.files(filePath,pattern = ".wav"))){
        combname<-i
        
        resltsTab<-runRavenDetector(m,filePathNoTemp,combname,resltsTab)
      }
      
    }
    
  }
  
  return(resltsTab)
}

#############################################################################


###################### ACTUALLY RUN THE SCRIPT ##############################


#############################################################################

#dumb conditional so I don't have to change path from machine to machine
if(dir.exists("C:/Users/ACS-3")){
  user<-"ACS-3"
  drivepath<-"F:/"
  gitPath<-paste(drivepath,"RavenBLEDscripts/",sep="")
  
}else if(dir.exists("C:/Users/danby456")){
  user<-"danby456"
  drivepath<-"E:/"
  gitPath<-paste(drivepath,"RavenBLEDscripts/",sep="")
  
}else if(dir.exists("C:/Users/Daniel.Woodrich")){
  user<-"Daniel.Woodrich"
  drivepath<-"E:/"
  gitPath<-"//nmfs/akc-nmml/CAEP/Acoustics/Projects/Dans Detectors/RavenBLEDscripts/"
}

startcombpath<-paste(drivepath,"Combined_sound_files/",sep="")
BLEDpath<-paste("C:/Users/",user,"/Raven Pro 1.5/Presets/Detector/Band Limited Energy Detector/",sep="")
ravenpath<-paste("C:/Users/",user,"/Raven Pro 1.5",sep="")
outputpath<-paste(drivepath,"DetectorRunOutput/",sep="")
outputpathfiles<-paste(drivepath,"DetectorRunFiles/",sep="")

###############input parameters
#detector control
ControlTab<-read.csv(paste(gitPath,"Data/CallParams/Detector_control.csv",sep=""))
ControlTab[,3]<-as.character(ControlTab[,3])

#General
runname<- ControlTab[which(ControlTab[,2]=="runname"),3]
spec<-str_split(ControlTab[which(ControlTab[,2]=="spec"),3],",",simplify=TRUE)
runGT<-ControlTab[which(ControlTab[,2]=="runGT"),3]
if(runGT=="y"){
runGTsections<-str_split(ControlTab[which(ControlTab[,2]=="runGTsections"),3],",",simplify=TRUE)
addToMaster<<-ControlTab[which(ControlTab[,2]=="addToMaster"),3] 
}else{
runGTsections<-c("n","n","n")
addToMaster<-"n"
}
useMasterGT<<-ControlTab[which(ControlTab[,2]=="useMasterGT"),3] 

runNEW<-ControlTab[which(ControlTab[,2]=="runNEW"),3]
#NEW
if(runNEW=="y"){
  NEWmoorings<- str_split(ControlTab[which(ControlTab[,2]=="NEWmoorings"),3],",",simplify=TRUE)
  NEWsf<-str_split(ControlTab[which(ControlTab[,2]=="NEWsf"),3],",",simplify=TRUE)  
  NEWpath<-str_split(ControlTab[which(ControlTab[,2]=="NEWpath"),3],",",simplify=TRUE)
  NEWpath2<-str_split(ControlTab[which(ControlTab[,2]=="NEWpath2"),3],",",simplify=TRUE)
  NEWsourceFormat<-str_split(ControlTab[which(ControlTab[,2]=="NEWsourceFormat"),3],",",simplify=TRUE)
}

fileCombinesize<-as.numeric(ControlTab[which(ControlTab[,2]=="fileCombinesize"),3] )
fileCombinesize2ndIt<-as.numeric(ControlTab[which(ControlTab[,2]=="fileCombinesize2ndIt"),3] )
onlyPopulate<-ControlTab[which(ControlTab[,2]=="onlyPopulate"),3] 
parallelType<- ControlTab[which(ControlTab[,2]=="parallelType"),3]

#model
CV<- ControlTab[which(ControlTab[,2]=="CV"),3]
TPRthresh<- str_split(ControlTab[which(ControlTab[,2]=="TPRthresh"),3],",",simplify=TRUE)  
modelType<- ControlTab[which(ControlTab[,2]=="modelType"),3]
modelMethod<- ControlTab[which(ControlTab[,2]=="modelMethod"),3]
numFeatures<- as.numeric(ControlTab[which(ControlTab[,2]=="numFeatures"),3])


#assign runname and make run folder 
runname<-paste(runname,gsub("\\D","",Sys.time()),sep="_")
dir.create(paste(outputpath,runname,sep=""))

###########################make txt file of params for run:
write.csv(ControlTab,paste(outputpath,runname,"/Detector_control.csv",sep=""),row.names = FALSE)

################Script function

if(runGT=="y"){
  runRavenGT<-runGTsections[1]
  runProcessGT<-runGTsections[2]
  runTestModel<-runGTsections[3]
}else{
  runRavenGT<-'n'
  runProcessGT<-'n'
  runTestModel<-'n'
}

whichRun<-"GT"

##################start script#################
if(runRavenGT=="y"){
  
resltsTabF <- NULL
MoorInfoMspec<-NULL
for(s in spec){
    
  #populate global env with species specific variables
  loadSpecVars(s)
  
  #############################
  
  MoorInfo<-makeMoorInfo(GTmoorings,GTsf,GTpath,GTpath2,GTsourceFormat,s)
  MoorInfoMspec<-rbind(MoorInfoMspec,MoorInfo)
  
  resltsTabS<-combineDecRaven()
  resltsTabF<-rbind(resltsTabF,resltsTabS)
  
  #Save raven output 
  write.csv(resltsTabS,paste(paste(outputpathfiles,s,"Unprocessed_GT_data/",sep=""),runname,"_UnprocessedGT.csv",sep=""),row.names = FALSE)
  write.csv(MoorInfo,paste(paste(outputpathfiles,s,"Unprocessed_GT_data/",sep=""),runname,"_UnprocessedGTMoorInfo.csv",sep=""),row.names = FALSE,col.names = FALSE)
  
}
  
if(spec>1){
write.csv(resltsTabF,paste(paste(outputpathfiles,"Unprocessed_GT_data/",sep=""),runname,"_UnprocessedGT.csv",sep=""),row.names = FALSE)
write.csv(MoorInfoMspec,paste(paste(outputpathfiles,"Unprocessed_GT_data/",sep=""),runname,"_UnprocessedGTMoorInfo.csv",sep=""),row.names = F,col.names = FALSE)
}

MoorInfo<-MoorInfoMspec
resltsTab<-resltsTabF

}else if(runGT=="y"&runProcessGT=="y"){
  
  if(length(spec)==1){
  recentTab<-file.info(list.files(paste(outputpathfiles,spec,"Unprocessed_GT_data/",sep=""), full.names = T,pattern = "GT.csv"))
  recentPath<-rownames(recentTab)[which.max(recentTab$mtime)]
  resltsTab<-read.csv(recentPath) #produces most recently modifed file 
  
  recentTab<-file.info(list.files(paste(outputpathfiles,spec,"Unprocessed_GT_data/",sep=""), full.names = T,pattern = "Info.csv"))
  recentPath<-rownames(recentTab)[which.max(recentTab$mtime)]
  MoorInfo<-read.csv(recentPath) #produces most recently modifed file 
  
  
  }else if(length(spec)>1){
  recentTab<-file.info(list.files(paste(outputpathfiles,"Unprocessed_GT_data/",sep=""), full.names = T,pattern = "GT.csv"))
  recentPath<-rownames(recentTab)[which.max(recentTab$mtime)]
  resltsTab<-read.csv(recentPath) #produces most recently modifed file 
  
  recentTab<-file.info(list.files(paste(outputpathfiles,"Unprocessed_GT_data/",sep=""), full.names = T,pattern = "Info.csv"))
  recentPath<-rownames(recentTab)[which.max(recentTab$mtime)]
  MoorInfo<-read.csv(recentPath) #produces most recently modifed file 
  
  }

}

if(useMasterGT=="y"){
  if(addToMaster=="y"){
    MoorInfoMaster<-read.csv(paste(gitPath,"Data/MoorInfo.csv",sep=""))
    moorsCompleted<-paste(MoorInfoMaster[,9],MoorInfoMaster[,10])
    if(any(!(paste(resltsTabTemp$Species,resltsTabTemp$MooringID) %in% moorsCompleted))){
      runProcessGT<-"n"
    }
  }else{
    runProcessGT<-"n"
  }
}

if(runProcessGT=="y"){

  
#Combine and configure spread detectors. 
  
DetecTab<-NULL
for(s in spec){
  
  loadSpecVars(s)
  
  DetecTab<-rbind(DetecTab,process_data())
  #resltsTab<-NULL

}

#depending on how RF handles multiple species classification- add step here that isolates boxes with overlap between species and allows them to be treated differently

DetecTab$detectionType<-0

#assess accuracy compared to GT 

GTset=NULL

GTtot=0
GTtot2=NULL
TPtot=NULL
MoorCor=NULL

for(s in spec){
GT<-list()
DetecTab2<-DetecTab[which(DetecTab$Species %in% s),]
p=1
for(f in 1:length(unique(DetecTab2$MooringName))){
  DetecVar<-DetecTab2[which(DetecTab2$MooringName==unique(DetecTab2$MooringName)[f]),]
  for(g in 1:length(unique(DetecVar$MooringID))){
    GT[[p]] <- read.delim(paste(gitPath,"Data/Selection tables/",s,"/",DetecVar$MooringName[1],"Sum/",unique(DetecVar$MooringID)[g],".txt",sep=""))
    GT[[p]] <- GT[[p]][GT[[p]]$View=="Spectrogram 1",]
    p=p+1
  }
}

#Define table for later excel file export.

colClasses = c("character","character","character","character","character","numeric","numeric","numeric","numeric", "numeric","numeric","numeric","numeric","character","character","character","character","character","character","character","character","character","character","character","character","numeric","numeric","character")
detecEvalFinal <- read.csv(text="Species, Moorings, Detectors, DetType, RunName, detTotal, numTP, numFP, numFN, TPR, FPR, TPdivFP,AUCav,CV_TPRthresh,Greatcall_goodcall,Max_modifier_penalty,ZerosAllowed,GroupSize,DownsweepThresh_DownsweepDiff,SkipAllowance,GroupInterval,TimeDiff,TimeDiffself,MinMaxDur,numDetectors,FO,LMS,Notes", colClasses = colClasses)

for(v in 1:length(unique(paste(DetecTab2$Species,DetecTab2$MooringID)))){
  print(paste("Comparing ground truth of",unique(paste(DetecTab2$Species,DetecTab2$MooringID))[v],"with final detector"))   
  MoorVar<-DetecTab2[which(paste(DetecTab2$Species,DetecTab2$MooringID)==unique(paste(DetecTab2$Species,DetecTab2$MooringID))[v]),]
  MoorVar$meantime<-(MoorVar[,4]+MoorVar[,5])/2
  MoorVar<-MoorVar[order(MoorVar$meantime),]
  MoorVar$Selection<-seq(1:nrow(MoorVar))
  
  #Define useful comlumns in both MoorVar and GT
  GT[[v]]$meantime<-(as.numeric(GT[[v]][,4])+as.numeric(GT[[v]][,5]))/2
  GT[[v]]$View<-as.character(GT[[v]]$View)
  GT[[v]]$detectionType<-0
  
  #sum rows of GT to do stats with later
  GTtot<-sum(GTtot,nrow(GT[[v]]))
  
  #change GT names to match what Raven accepts
  colnames(GT[[v]]) <- c("Selection","View","Channel","Begin Time (s)","End Time (s)","Low Freq (Hz)","High Freq (Hz)", "meantime", "detectionType")
  colClasses = c("numeric","character", "numeric","numeric","numeric", "numeric","numeric","numeric","numeric")
  
  #define tables for Tp/FN/FPs
  OutputCompare <- read.csv(text="Selection,View,Channel,Begin Time (s),End Time (s),Low Freq (Hz),High Freq (Hz), meantime, detection type", colClasses = colClasses)
  colnames(OutputCompare)<- c("Selection","View","Channel","Begin Time (s)","End Time (s)","Low Freq (Hz)","High Freq (Hz)", "meantime", "detectionType")
  OutputCompare2 <- read.csv(text="Selection,View,Channel,Begin Time (s),End Time (s),Low Freq (Hz),High Freq (Hz),meantime, detection type", colClasses = colClasses)
  colnames(OutputCompare2)<- c("Selection","View","Channel","Begin Time (s)","End Time (s)","Low Freq (Hz)","High Freq (Hz)", "meantime", "detectionType")
  
  #Identify TPs in data. Criteria is if meantime of detection is between that of GT start and end time
  p=1
  for(h in 1:nrow(MoorVar)){
    gvec <- which(GT[[v]]$meantime<(MoorVar$meantime[h]+Maxdur+1)&GT[[v]]$meantime>(MoorVar$meantime[h]-Maxdur-1))
    if(length(gvec)>0){
    for(g in min(gvec):max(gvec)){
      if((MoorVar$meantime[h]>GT[[v]][g,4]) & (MoorVar$meantime[h]<GT[[v]][g,5])|(GT[[v]]$meantime[g]>MoorVar[h,4] & GT[[v]]$meantime[g]<MoorVar[h,5])){
        OutputCompare[p,]<-MoorVar[h,c(1:7,19,20)]
        OutputCompare$detectionType<-"TP"
        p=p+1
      }
    }
  }
  }
  #remove duplicates:
  OutputCompare<-OutputCompare[which(!duplicated(OutputCompare$Selection)),]

  
  #Identify and add FPs. if selection in MoorVar row does not match that in Output compare, add it to Output compare under designation FP.  
  OutputCompare <- rbind(OutputCompare,MoorVar[-which(MoorVar$Selection %in% OutputCompare$Selection),c(1:7,19,20)])
  OutputCompare[which(OutputCompare$detectionType!="TP"),9]<-"FP"
  
  #Add rows where GT meantime was in between 
  p=1
  for(h in 1:nrow(GT[[v]])){
    gvec <- which(MoorVar$meantime<(GT[[v]]$meantime[h]+Maxdur+1)&MoorVar$meantime>(GT[[v]]$meantime[h]-Maxdur-1))
    if(length(gvec)>0){
      for(g in min(gvec):max(gvec)){
        if((GT[[v]][h,8]>MoorVar[g,4] & GT[[v]][h,8]<MoorVar[g,5])|((MoorVar$meantime[g]>GT[[v]][h,4]) & (MoorVar$meantime[g]<GT[[v]][h,5]))){
          OutputCompare2[p,]<-GT[[v]][h,]
          OutputCompare2[p,9]<-"TP truth"
          p=p+1
        }
      }
    }
  }
  
  OutputCompare2<-OutputCompare2[which(!duplicated(OutputCompare2$Selection)),]
  
  
  #Identify and add FNs. if selection in GT row does not match that in OutputCompare2, add it to Output compare under designation FN.  
  OutputCompare2 <- rbind(OutputCompare2,GT[[v]][-which(GT[[v]]$Selection %in% OutputCompare2$Selection),])
  OutputCompare2[which(OutputCompare2$detectionType!="TP truth"),9]<-"FN"
  
  #Combine tables and remove GT TPs from dataset. Save TPs in other vector 
  OutputCompareDet<-OutputCompare
  OutputCompareDet$meantime<-(OutputCompareDet[,4]+OutputCompareDet[,5])/2  
  OutputCompareDet<-OutputCompareDet[order(OutputCompareDet$meantime),]
  TPFPs<-OutputCompareDet$detectionType
  OutputCompare<-rbind(OutputCompare,OutputCompare2)
  OutputCompare$meantime<-as.numeric(OutputCompare$meantime)
  OutputCompare<-OutputCompare[order(OutputCompare$meantime),]
  OutputCompare[which(OutputCompare$detectionType=="TP truth"),9]<-"x"
  OutputCompare <- subset(OutputCompare,detectionType!="x")
  OutputCompare$Selection<-seq(1:nrow(OutputCompare))
  
  #Ready table for Raven and save. 
  colnames(OutputCompare)[8]<-'TP/FP/FN'
  OutputCompare<- OutputCompare[,-8]
  OutputCompare<- OutputCompare[,1:8]
  
  OutputCompare$Selection<-seq(1:nrow(OutputCompare))
  write.table(OutputCompare,paste(outputpath,runname,"/",s,"_",MoorVar$MooringID[1],"_TPFPFN_Tab_Ravenformat.txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
  
  #Make summary table of statistics for table comparison. 
  numTP <- nrow(OutputCompare[which(OutputCompare[,8]=="TP"),])
  numFN <- nrow(OutputCompare[which(OutputCompare[,8]=="FN"),])
  numFP <- nrow(OutputCompare[which(OutputCompare[,8]=="FP"),])
  numTPtruth<- nrow(GT[[v]])
  detTotal<-numTP+numFP

  #store this for comparison with full mooring later
  TPtot<-c(TPtot,numTP)
  MoorCor<-c(MoorCor,unique(paste(DetecTab2$Species,DetecTab2$MooringID))[v])
  GTtot2<-c(GTtot2,nrow(GT[[v]]))
  
  TPR <- numTP/(numTP+numFN)
  FPR <- numFP/(numTP+numFP)
  TPdivFP<- numTP/numFP

  #save stats and parameters to excel file
  detecEval<-detecEvalFinal[0,]
  detecEval[1,]<-c(s,unique(paste(DetecTab2$Species,DetecTab2$MooringID))[v],paste(detnum,paste(detlist2,collapse="+"),sep=";"),"spread",runname,detTotal,numTP,numFP,numFN,TPR,FPR,TPdivFP,NA,NA,NA,NA,paste(allowedZeros,collapse=","),paste(grpsize,collapse=","),"depreciated",paste(detskip,collapse=","),paste(groupInt,collapse=","),NA,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")
  detecEvalFinal <- rbind(detecEvalFinal,detecEval)
  
  MoorVar$detectionType<-TPFPs
  
  
  #combine all data:
  
  GTset<-rbind(GTset,MoorVar)
}

#Make summary table of whole run statistics for table comparison. 
numTP <- sum(as.numeric(detecEvalFinal[,7]))
numFN <- sum(as.numeric(detecEvalFinal[,9]))
numFP <- sum(as.numeric(detecEvalFinal[,8]))
numTPtruth<- GTtot
detTotal<-numTP+numFP


TPR <- numTP/(numTP+numFN)
FPR <- numFP/(numTP+numFP)
TPdivFP<- numTP/numFP

#save stats and parameters to excel file
detecEval<-detecEvalFinal[0,]
detecEval[1,]<-c(s,"all",paste(detnum,paste(detlist2,collapse="+"),sep=";"),"spread",runname,detTotal,numTP,numFP,numFN,TPR,FPR,TPdivFP,NA,NA,NA,NA,paste(allowedZeros,collapse=","),paste(grpsize,collapse=","),"depreciated",paste(detskip,collapse=","),paste(groupInt,collapse=","),NA,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")

detecEvalFinal <- rbind(detecEvalFinal,detecEval)
 
detecEval2<-read.csv(paste(outputpath,"DetectorRunLog.csv",sep=""))
detecEvalFinal<-rbind(detecEval2,detecEvalFinal)

write.csv(detecEvalFinal,paste(outputpath,"DetectorRunLog.csv",sep=""),row.names=FALSE)

###################

}

#define this table to compare counts after running model
TPtottab<-data.frame(TPtot,GTtot2,MoorCor)

#translate response to binary 
GTset$detectionType<- as.character(GTset$detectionType)
GTset[which(GTset$detectionType=="TP"),which(colnames(GTset)=="detectionType")]<-1
GTset[which(GTset$detectionType=="FP"),which(colnames(GTset)=="detectionType")]<-0
GTset[which(GTset$detectionType=="FN"),which(colnames(GTset)=="detectionType")]<-2

#######1 mooring test######
if(length(spec)>1){
#GTset<-GTset[which(GTset$`soundfiles[n]`=="BS15_AU_02a_files1-104.wav"),]

GTset<-dataConflicts(GTset)


#if multiclass: throw out double positives from the ground truth. Comment out line if choose to go with binary. 

GTset<-GTset[-grep(",",GTset$Species),]

}else{
  GTset$meanfreq<- (GTset$`Low Freq (Hz)`+GTset$`High Freq (Hz)`)/2
  GTset$freqrange<- (GTset$`High Freq (Hz)`-GTset$`Low Freq (Hz)`)
  GTset$meantime<- (GTset$`Begin Time (s)`+GTset$`End Time (s)`)/2
  
  GTset$Selection<-seq(1,nrow(GTset))
}

#"vectorize" GTset frame. 

GTset$combID<-as.factor(paste(GTset$Species,GTset$MooringID))
dataMat<- data.matrix(cbind(GTset[,c(23,4:7)],as.numeric(GTset$RTFb+GTset$FileOffsetBegin),as.numeric(GTset$detectionType)))
print("extracting features from FFT of each putative call")

stop()
dataMat<-spectral_features(dataMat)

GTset$combID<-NULL
GTset<-cbind(GTset,dataMat[,c(8:ncol(dataMat))])

GTset<-apply(GTset,2,function(x) unlist(x))

if(length(spec)==1){
dir.create(paste(outputpathfiles,spec,"Processed_GT_data/",sep=""))
dir.create(paste(outputpathfiles,spec,"TPtottab/",sep=""))

write.csv(GTset,paste(outputpathfiles,spec,"Processed_GT_data/",runname,"_processedGT.csv",sep=""),row.names = F)
write.csv(TPtottab,paste(outputpathfiles,spec,"TPtottab/",runname,"_processedGT.csv",sep=""),row.names = F)

}else if(length(spec)>1){
  dir.create(paste(outputpathfiles,"Processed_GT_data/",sep=""))
  dir.create(paste(outputpathfiles,"TPtottab/",sep=""))
  
  write.csv(GTset,paste(outputpathfiles,"Processed_GT_data/",runname,"_processedGT.csv",sep=""),row.names = F)
  write.csv(TPtottab,paste(outputpathfiles,"TPtottab/",runname,"_processedGT.csv",sep=""),row.names = F)
  
}

GTset<-GTset[,c(1,4:ncol(GTset))]
GTset<-data.frame(GTset)
GTset$detectionType<-as.factor(GTset$detectionType)

}else if(runGT=="y"){
  
inputGT()
  
}

###############################################


#RUN TEST RANDOM FOREST MODEL 

#################

if(runTestModel=="y"){
  
modData<-dataArrangeModel(GTset)

if(length(spec)>1){
  mSpec<-unique(substr(modData[[1]]$detectionType,1,2))
}else if(length(spec==1)){
  mSpec<-spec
}

if(modelType=="rf"){
  modelOutput<-runRandomForest(modData[[1]],modData[[1]])
}else if(modelType=='orf'){
  modelOutput<-runObliqueRandomForest(dataForModel,method=modelMethod)
f}


#end model function. export dataset, cutmean

#return dataset from random forest 

######################
data3<-cbind(modData[[2]],modelOutput[[1]])

CUTmean<-modelOutput[[2]]

nospec<-ncol(data3)-(length(spec)*3)

for(s in 1:length(spec)){
  
dataSPEC<-data3[which(data3$Species==spec[s]),]
CUTmeanspec<-CUTmean[[s]]


probIndex<-nospec+(3*s-2)
varIndex<-nospec+(3*s-1)
  
#get these out of the way so I can remove mostly useless std.err and n from dataframe
plot(dataSPEC[which(dataSPEC$detectionType==paste(spec[s],0)),probIndex],dataSPEC[which(dataSPEC$detectionType==paste(spec[s],0)),varIndex], col = "red",main=spec[s],cex=0.25)
abline(v=CUTmeanspec)

plot(dataSPEC[which(dataSPEC$detectionType==paste(spec[s],1)),probIndex],dataSPEC[which(dataSPEC$detectionType==paste(spec[s],1)),varIndex], col = "blue",main=spec[s],cex=0.25)
abline(v=CUTmeanspec)

plot(dataSPEC[,probIndex],dataSPEC[,varIndex], col = ifelse(dataSPEC$detectionType==paste(spec[s],1),'blue','red'),main=spec[s],cex=0.25)
abline(v=CUTmeanspec)
}

#pickup here or leave it for now. 
newSpec<-NULL
for(s in 1:length(spec)){
  newSpec<-c(newSpec,(nospec+3*s-2))
}

data3<-data3[,c(1:nospec,newSpec)]  

#put columns excluded from model back in ?

#change data columns from factor to numeric that should be numeric:
for(b in c(1,2,3,4,5,10,11)){
  data3[,b]<-as.numeric(as.character(data3[,b])) 
  
}

for(b in c(13,14)){
  data3[,b]<-as.integer(as.numeric(as.POSIXlt(data3[,b])))
}

#change to factor
data3[,15]<-as.factor(data3[,15])

#adaptively combine detections based on probability
data3Labs<-factorLevels(data3)

data3$Unq_Id<-as.numeric(as.factor(paste(data3$Species,data3$MooringID)))

data3Mat<- data.matrix(data3)

data3Mat<-adaptive_compare(data3Mat) 

#simulate context over time using probability scores.Make new tab with selection ID and all supporting variables- don't pin on data3mat. Can use in place of prob if wanted. 
CS_output<-context_sim(data3Mat)

#remove missing rows:
data3Mat<-data3Mat[order(data3Mat[,1]),]
data3<-data.frame(data3Mat)
data3$Unq_Id<-NULL
data3<-applyLevels(data3,data3Labs)

#graphs
AUCadj<-NULL
for(s in 1:length(spec)){
  
dataSPEC<-data3[which(data3$Species==spec[s]),]
CUTmeanspec<-CUTmean[[s]]


dataSPEC$detectionType<-as.numeric(as.factor(dataSPEC$detectionType))-1

pp2<-as.vector(dataSPEC[,nospec+s])
ll2<-dataSPEC$detectionType
predd2<-prediction(pp2,ll2)
perff2<-performance(predd2,"tpr","fpr")

#with permutations on probs
plot(perff2, avg = "threshold",  xaxs="i", yaxs="i", spread.scale=2,
     lwd = 2, main = paste("Threshold avg",spec[s]),colorize=T)
abline(a=0, b= 1)
auc.perf = performance(predd2, measure = "auc",plot=F)
print(auc.perf@y.values)

AUCadj<-c(AUCadj,auc.perf@y.values)

#plot of probabilities after context sim:
for(m in unique(dataSPEC[,6])){
data3Matmoors<-dataSPEC[which(dataSPEC[,6]==m),]
CS_outputDate<-CS_output[which(CS_output[,1] %in% as.numeric(data3Matmoors[,1])),]
#pointsDate<-data3Matmoors$Begin.Time..s.
pointsDate<-as.POSIXlt((as.numeric(data3Matmoors$RTFb)+data3Matmoors$FileOffsetBegin), origin="1970-01-01",tz="UTC")
plot(x=pointsDate,y=data3Matmoors[,nospec+1], col=as.factor(data3Matmoors$detectionType),main=paste(spec[s],m))
abline(h=CUTmean[s],col="red")
abline(h=0.5,lty=3)
#lines(lowess(data3Matmoors[,pos+5]))
#lines(lowess(data3Matmoors[,pos+4]))
#lines(lowess(data3Matmoors[,pos+2]))
lines(pointsDate,(CS_outputDate[,(s*4)+1]),col="blue") #backwards through data 
lines(pointsDate,(CS_outputDate[,(s*6)+1]),col="orange") #forwards through data
lines(pointsDate,(CS_outputDate[,(s*7)+1]),col="green")
}

}

#data3$detectionType<-as.factor(data3$detectionType)
#see freq breakdown of calls 
#cdplot(data3datFrame[,7] ~ data3datFrame[,8], data3datFrame, col=c("cornflowerblue", "orange"), main="Conditional density plot") #meanfreq
#write data to drive
after_model_write(data3) #need to change to vector 

beep(10)

}

###############################################


#EMPLOY MODEL ON FULL DATASETS

################################################

if(runNEW=="y"){
  
whichRun<-"NEW"

resltsTabF <- NULL
MoorInfoMspec<-NULL
for(s in spec){

  #populate global env with species specific variables
  loadSpecVars(s)
    
  #############################
    
  MoorInfo<-makeMoorInfo(NEWmoorings,NEWsf,NEWpath,NEWpath2,NEWsourceFormat,s)
  MoorInfoMspec<-rbind(MoorInfoMspec,MoorInfo)
    
  resltsTabS<-combineDecRaven()
  resltsTabF<-rbind(resltsTabF,resltsTabS)
}

  MoorInfo<-MoorInfoMspec

  resltsTab<-resltsTabF  
  resltsTabF<-NULL
  resltsTabS<-NULL
  resltsVar<-NULL
  
  #Combine and configure spread detectors. 
  DetecTab<-NULL
  for(s in spec){
    
    loadSpecVars(s)
    
    DetecTab<-rbind(DetecTab,process_data())
    #resltsTab<-NULL
  }
  
  DetecTab$detectionType<-0
  
  #print table that has buffer to allow for easily playing shorter GS

DetecTab<-dataConflicts(DetecTab)

DetecTab$combID<-as.factor(paste(substr(DetecTab$Species,1,2),DetecTab$MooringID)) #simplify species ID for correct indexing
dataNEW<- data.matrix(cbind(DetecTab[,c(23,4:7)],as.numeric(DetecTab$RTFb+DetecTab$FileOffsetBegin),as.numeric(DetecTab$detectionType)))
print("extracting features from FFT of each putative call")

dataNEW<-spectral_features(dataNEW)

dataNEW<-data.frame(dataNEW)
DetecTab$combID<-NULL
DetecTab<-cbind(DetecTab,dataNEW[,c(8:ncol(dataNEW))])

#DetecTab<-apply(DetecTab,2,function(x) unlist(x))

DetecTab<-DetecTab[,c(1,4:ncol(DetecTab))]
DetecTab<-data.frame(DetecTab)

modData<-dataArrangeModel(DetecTab)

modData[[1]]$detectionType<-modData[[1]]$detectionType[1]
modData[[1]]$detectionType<-droplevels(modData[[1]]$detectionType)

MoorInfosave<-MoorInfo

inputGT()


#eliminate all but RW for now to see if it generalizes better. 

MoorInfo<-MoorInfosave

GTData<-dataArrangeModel(GTset)

#temporarily exlude BS2 since can't have new factor levels in new data. 
#modData[[1]]<-modData[[1]][which(modData[[1]]$MooringCode!="BS2"),]
#modData[[2]]<-modData[[2]][which(modData[[1]]$Selection %in% modData[[2]]$Selection),]

for(i in 1:ncol(modData[[1]])){

levels(modData[[1]][,i])<-levels(GTData[[1]][,i]) 

}

names(modData[[1]])<-names(GTData[[1]])

CV<-10

#mSpec<-c(spec,unique(substr(GTData[[1]]$detectionType,1,2))[which(!unique(substr(GTData[[1]]$detectionType,1,2)) %in% spec)])

mSpec<-unique(substr(GTData[[1]]$detectionType,1,2))

modelOutput<-runRandomForest(modData[[1]],GTData[[1]])

data3<-cbind(modData[[2]],modelOutput[[1]])
CUTmean<-modelOutput[[2]]

nospec<-ncol(data3)-(length(mSpec)*3)

newSpec<-NULL
for(s in 1:length(mSpec)){
  newSpec<-c(newSpec,(nospec+3*s-2))
}

data3<-data3[,c(1:nospec,newSpec)]  

#put columns excluded from model back in ?

#change data columns from factor to numeric that should be numeric:
for(b in c(1,2,3,4,5,10,11)){
  data3[,b]<-as.numeric(as.character(data3[,b])) 
  
}

for(b in c(13,14)){
  data3[,b]<-as.integer(as.numeric(as.POSIXlt(data3[,b])))
}

#dataPostModel<-formatModelData(data3)

#data3<-dataPostModel[[1]]

AUCadj<-"NA"

after_model_write(data3) #need to change to vector 

stop()



stop()

}
#Define table for later excel file export. 
colClasses = c("character","character","character","character","character","numeric","numeric", "numeric","numeric","numeric","numeric","numeric","character","character","character","character","character","character","character","character","character","character","character","character","numeric","numeric","character")
detecEvalFinal <- read.csv(text="Species, Moorings, Detectors, DetType, RunName, numTP, numFP, numFN, TPR, FPR, TPdivFP,AUCav,CV_TPRthresh,Greatcall_goodcall,Max_modifier_penalty,ZerosAllowed,GroupSize,DownsweepThresh_DownsweepDiff,SkipAllowance,GroupInterval,TimeDiff,TimeDiffself,MinMaxDur,numDetectors,FO,LMS,Notes", colClasses = colClasses)

#commented out 

#MoorTab<-NULL
#MoorVar<-NULL
#for(v in 1:length(unique(DetecTab2$Mooring))){
#  print(paste("Adding interference detector variables to",sort(unique(DetecTab2$Mooring))[v]))   
#  MoorVar<-DetecTab2[which(DetecTab2$Mooring==sort(unique(DetecTab2$Mooring))[v]),]
#  
#  #Define useful comlumns in MoorVar
#  sound.files <- MoorVar[,12]
#  MoorVar <- MoorVar[,c(1:7,13,14:17)]
#  MoorVar<-cbind(sound.files,MoorVar)
#  
#    #compare detections to sources of interference
#    if(interfere=="y"){
#      for(n in 1:max(resltsTabInt$detectorCount)){
#        MoorInt<-resltsTabInt[which(resltsTabInt$Mooring==sort(unique(DetecTab2$Mooring))[v]),]
#        MoorVar[,n+14]<-0
#        colnames(MoorVar)[length(MoorVar)]<-paste(resltsTabInt[which(resltsTabInt$detectorCount==n),10])[1]
#        MoorInt<-MoorInt[which(MoorInt$detectorCount==n),]
#        print(paste("       Compare with detector",MoorInt[1,10]))   
#        for(g in 1:nrow(MoorVar)){
#          hvec <- which(MoorInt$meantime<(MoorVar$meantime[g]+2)&MoorInt$meantime>(MoorVar$meantime[g]-2))
#          if(length(hvec>0)){
#          for(h in min(hvec):max(hvec)){
#            if((MoorInt[h,4]<MoorVar[g,9] & MoorInt[h,5]> MoorVar[g,9])|(MoorInt[h,4]> MoorVar[g,5] & MoorInt[h,5]< MoorVar[g,6])){
#              MoorVar[g,n+14]<-1
#            }
#          }
#        }
#      }
#      }
#    }
  
  #write.table(MoorVar[,2:7],paste(outputpath,runname,"/",unique(DetecTab2$Mooring)[v],"FINAL_Summary_",dettype,"_Ravenformat",".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
  
#  MoorTab<-rbind(MoorTab,MoorVar)
  #MoorVar$Moorpred<-c(MoorVar$Moorpred,predict(data.rf,MoorVar,type="prob"))
#}

#findata<-MoorTab

#make interference columns into factors
#if(length(findata)>17){
#  for(n in 18:length(findata)){
#    findata[,n]<-as.factor(findata[,n])
#  }
#}
#
moorlib<-cbind(seq(1,length(unique(findata$sound.files)),1),as.character(sort(unique(findata$sound.files))),as.character(sort(unique(findata$Mooring))))

findata$sound.files<-as.factor(findata$sound.files)
findataMat<- data.matrix(findata[c(13,4,5,6,7)])
findataMat<-spectral_features(findataMat,moorlib,2,1)

colnames(findataMat)<-c(colnames(findataMat)[1:5],c(letters,strrep(letters,2))[1:(ncol(findataMat)-5)])

#Generate and run a set amount of models from the original GT data. Probabilities are averaged for each mooring. 
if(runProcessGT=="n"){
  recentTab<-file.info(list.files(paste(outputpathfiles,"Processed_GT_data/",sep=""), full.names = T))
  recentPath<-rownames(recentTab)[which.max(recentTab$mtime)]
  data<-read.csv(recentPath) #produces most recently modifed file 
  colnames(data)[1]<-"soundfiles[n]"
  
  recentTab<-file.info(list.files(paste(outputpathfiles,"TPtottab/",sep=""), full.names = T))
  recentPath<-rownames(recentTab)[which.max(recentTab$mtime)]
  TPtottab<-read.csv(recentPath) #produces most recently modifed file 
  
  GTset<-data[,c(9,13:length(data))]
}
GTset$detectionType<-as.factor(GTset$detectionType)
#GTset<-data[,c(7,11:length(data))] #if running full workflow change so that GTset path is reduced to only detectionType and variables

colnames(GTset)<-c(colnames(GTset)[1],c(letters,strrep(letters,2))[1:(ncol(GTset)-1)])

#generate models and apply to whole dataset. Should keep model parameters the same as when you assessed accuracy to have an idea of reliability. 
CUTvec=NULL
for(p in 1:CV){
  print(paste("model",p))
  train<-splitdf(GTset,weight = 2/3)
  #apply model to data
  data.rf<-randomForest(formula=detectionType ~ .,data=train[[1]],mtry=7,na.action=na.roughfix)
  pred<-predict(data.rf,findataMat[,6:ncol(findataMat)],type="prob")
  #assess same model performance on GT data 
  pred2<-predict(data.rf,train[[2]],type="prob")
  ROCRpred<-prediction(pred2[,2],train[[2]]$detectionType)
  perf = performance(ROCRpred, "tpr","fpr")
  
  TPR<-NULL
  TPR <- data.frame(cut=perf@alpha.values[[1]], tpr=perf@y.values[[1]])
  CUT <- max(TPR[which(TPR$tpr>=TPRthresh),1])
  CUTvec<-c(CUTvec,CUT)
  

  if(p==1){
    probstab<-pred[,2]
  }else{
    probstab<-data.frame(probstab,pred[,2])
  }  
}

probmean<-NULL
probstderr<-NULL
n<-NULL

for(x in 1:nrow(probstab)){
probmean[x]<-mean(as.numeric(probstab[x,]))
probstderr[x]<-std.error(as.numeric(probstab[x,]))
n[x]<-length(which(!is.na(probstab[x,2:length(probstab)])))
}

CUTmean<-mean(CUTvec)
CUTstd.err<-std.error(CUTvec)



findataMat<-cbind(findataMat[,1],seq(1:nrow(findataMat)),findataMat[,2:5],0,((findataMat[,4]+findataMat[,5])/2),(findataMat[,5]-findataMat[,4]),((findataMat[,2]+findataMat[,3])/2),findataMat[,6:ncol(findataMat)])
#create columns to be used later 
findataMat<-cbind(findataMat,probmean)
findataMat<-cbind(findataMat,probstderr)
findataMat<-cbind(findataMat,n)

pos<-ncol(findataMat)

findataMat<-adaptive_compare(findataMat,2)

#make value with NA 0
findataMat[!rowSums(!is.finite(findataMat)),]
findataMat[!is.finite(findataMat)] <- 0

findataMat<-context_sim(findataMat)

for(n in unique(findataMat[,1])){
  findataMatmoors<-findataMat[which(findataMat[,1]==n),]
  plot(findataMatmoors[,pos-2], col=findataMatmoors[,7],main=moorlib[which(moorlib[,1]==n),2])
  abline(h=CUTmean,col="red")
  abline(h=0.5,lty=3)
  #lines(lowess(data3moors[,pos+5]))
  #lines(lowess(data3moors[,pos+4]))
  #lines(lowess(data3moors[,pos+2]))
  lines((findataMatmoors[,pos+3]*18),col="blue") #backwards through data 
  lines((findataMatmoors[,pos+1]*18),col="orange") #forwards through data
  #lines(((pmax(findataMatmoors[,pos+1],findataMatmoors[,pos+3])*6)),col="green")
}


AUCadj<-NA
after_model_write(findataMat,moorlib,2)
#apply models and average probabilities. ASSUMPTIONS for stats to make sense: Full mooring is run, and all HG data is included in GT. 
MoorVar1<-NULL
run<-"no"
if(run=="yes"){
for(v in 1:length(unique(findataMat[,1]))){
  
  MoorVar1<-as.data.frame(findataMat[which(findataMat[,1]==sort(unique(findataMat[,1]))[v]),2])
  
  #MoorVar1$detectionType<-ifelse(((MoorVar1$probstderr<CUTstd.err)&(MoorVar1$probstmean>CUTmean)),"RFselected","RFrejected")
  MoorVar1[,7]<-ifelse(MoorVar1[,pos-2]>CUTmean,"RFselected","RFrejected")
  MoorVar1<-MoorVar1[which(MoorVar1[,7]=="RFselected"),]
  
  numTP<-TPtottab[which(TPtottab$MoorCor==sort(unique(findataMat[,1]))[v]),1]*TPRthresh
  numTPtruth<-TPtottab[which(TPtottab$MoorCor==sort(unique(findataMat[,1]))[v]),2]
  detTotal<-nrow(MoorVar1)
  numFP<-detTotal-numTP
  numFN<-numTPtruth-numTP
  
  TPR <- numTP/(numTP+numFN)
  FPR <- numFP/(numTP+numFP)
  TPdivFP<- numTP/numFP
  
  #save stats and parameters to excel file
  detecEvalFinal<-read.csv(paste(outputpath,"DetectorRunLog.csv",sep=""))
  detecEval<-detecEvalFinal[0,]
  if(dettype=="spread"|dettype=="combined"){
    detecEval[1,]<-c(spec,paste("full",moorlib[which(moorlib[,1]==sort(unique(findataMat[,1]))[v]),2]),paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,detTotal,numTP,numFP,numFN,TPR,FPR,TPdivFP,NA,NA,NA,NA,paste(allowedZeros,collapse=","),paste(grpsize,collapse=","),paste(downsweepCompMod,downsweepCompAdjust,sep=","),paste(detskip,collapse=","),paste(groupInt,collapse=","),NA,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")
  }else{
    detecEval[1,]<-c(spec,paste("full",moorlib[which(moorlib[,1]==sort(unique(findataMat[,1]))[v]),2]),paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,detTotal,numTP,numFP,numFN,TPR,FPR,TPdivFP,NA,NA,NA,NA,NA,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")   
  }
  
  detecEval2<-read.csv(paste(outputpath,"DetectorRunLog.csv",sep=""))
  detecEvalFinal<-rbind(detecEval2,detecEval)
  write.csv(detecEvalFinal,paste(outputpath,"DetectorRunLog.csv",sep=""),row.names=FALSE)
  
  MoorVar1<-MoorVar1[,c(2:8)]
  
  MoorVar2<-findataMat[which(findataMat[,1]==sort(unique(findataMat[,1]))[v]),]
  
  MoorVar3<-data.frame(MoorVar2$Selection,MoorVar2$FileOffsetBegin,MoorVar2$FileOffsetEnd,MoorVar2$`Low Freq (Hz)`,MoorVar2$`High Freq (Hz)`,MoorVar2$sound.files,MoorVar2$File,MoorVar2$probmean,MoorVar2$probstderr,MoorVar2$probn)
  colnames(MoorVar3)<-c("Selection","FileOffsetBegin","FileOffsetEnd","Low Freq (Hz)","High Freq (Hz)","Mooring","File","probs","probsstderr","probsn")
  
  write.table(MoorVar1,paste(outputpath,runname,"/",sub(".wav", "", moorlib[which(moorlib[,1]==sort(unique(findataMat[,1]))[v]),2]),"FINAL_Model_Applied_Ravenformat",".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
  write.table(MoorVar2,paste(outputpath,runname,"/",sub(" .wav", "", moorlib[which(moorlib[,1]==sort(unique(findataMat[,1]))[v]),2]),"FINAL_Model_Applied_probs",".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
  
}
  

#After we save and report the stats, could make a section here for analyzing probabilites and variance. 
  
#no avg
colfunc <- colorRampPalette(c("red", "green"))

plot(as.numeric(probmean),probstderr, col = ifelse(as.numeric(probstderr) > CUTstd.err,'red','green'))
plot(as.numeric(probmean),probstderr, col = ifelse(as.numeric(probmean) < CUTmean,'red','green'))

#this looks like cleanest portion of data- but how to subset to this while keeping a known TPR? Even if it takes a after the fact analysis, should explore only taking the "tail" of the data
plot(as.numeric(probmean),probstderr, col = ifelse(((as.numeric(probmean) < CUTmean)|(as.numeric(probstderr)>CUTstd.err)),'red','green'))

cor.test(as.numeric(probmean),probstderr)

}


