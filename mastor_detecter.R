#Evaluate and score custom BLED detector(s) in R. Sox must be installed and pathed to accordingly

###############################################


#RUN RAVEN DETECTORS


################################################

#install.packages("flightcallr")install.packages("randomForest")install.packages("seewave")install.packages("tuneR")install.packages("plotrix")install.packages("aod")install.packages("ggplot2")install.packages("usdm")install.packages("ROCR")install.packages("e1071")  install.packages("caret")  

library(e1071)  
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

duration_store<- function(){
  durVar<-NULL
  durVar<-data.frame(SFfp=character(),
                     SFsh=character(),
                     Duration=numeric(), 
                     CumDur=numeric(),
                     CombSF=character(), 
                     Mooring=character(),
                     stringsAsFactors=FALSE) 
for(f in 1:length(sound_filesfullpath)){
  durVar[f,1]<-sound_filesfullpath[f]
  durVar[f,2]<-sound_files[f]
  audio<-readWave(sound_filesfullpath[f], header=TRUE)
  durVar[f,3]<-round(audio$samples / audio$sample.rate, 2)
}
  durVar[,4]<-cumsum(durVar$Duration)
  durVar[,5]<-paste(m,"_files_entire",bigFile_breaks[b],".wav",sep="")
  durVar[,6]<-m

durTab<-rbind(durTab,durVar)
return(durTab)
}

spectral_features<- function(specdata,whichRun){
  if(whichRun==1){
    if(whiten=="y"){
      specpath<-paste(startcombpath,"/",spec,"/Bbandp",LMS*100,"x_FO",FO,"/",sep="")
    }else{
      specpath<-paste(startcombpath,"/",spec,"/No_whiten","/",sep="")
        }
  }else{
    if(whiten=="y" & moorType=="HG"){
      specpath<-paste(startcombpath,"/",spec,"/Entire_Bbandp",LMS*100,"x_FO",FO,"/",sep="")
    }else if(whiten=="y" & moorType!="HG"){
      specpath<-paste(startcombpath,"/Entire_full_Bbandp",LMS*100,"x_FO",FO,"/",sep="")
    }else if(whiten=="n" & moorType=="HG"){
      specpath<-paste(startcombpath,"/",spec,"/Entire_No_whiten","/",sep="")
    }else{
      specpath<-paste(startcombpath,"/Entire_full_No_whiten","/",sep="")
    }
  }

  
print("extracting spectral parameters")
for(z in 1:nrow(specdata)){
  foo <- readWave(paste(specpath,specdata[z,1],sep=""),specdata[z,5],specdata[z,6],units="seconds")
  foo.spec <- spec(foo, plot=F, PSD=T,ylim=c(yminn,ymaxx))
  foo.specprop <- specprop(foo.spec)
  foo.meanspec = meanspec(foo, plot=FALSE, ovlp=90)#not sure what ovlp parameter does but initially set to 90
  foo.autoc = autoc(foo, plot=F)
  foo.dfreq = dfreq(foo, plot=F, ovlp=90)
  specdata$rugosity[z] = rugo(foo@left / max(foo@left)) 
  specdata$crest.factor[z] = crest(foo)$C
  foo.env = seewave:::env(foo, plot=F) 
  specdata$temporal.entropy[z] = th(foo.env)
  specdata$shannon.entropy[z] = sh(foo.spec)
  specdata$spectral.flatness.measure[z] = sfm(foo.spec)
  specdata$spectrum.roughness[z] = roughness(foo.meanspec[,2])
  specdata$autoc.mean[z] = mean(foo.autoc[,2], na.rm=T)
  specdata$autoc.median[z] = median(foo.autoc[,2], na.rm=T)
  specdata$autoc.se[z] = std.error(foo.autoc[,2], na.rm=T)
  specdata$dfreq.mean[z] = mean(foo.dfreq[,2], na.rm=T)
  specdata$dfreq.se[z] = std.error(foo.dfreq[,2], na.rm=T)
  specdata$specprop.mean[z] = foo.specprop$mean[1]
  specdata$specprop.sd[z] = foo.specprop$sd[1]
  specdata$specprop.sem[z] = foo.specprop$sem[1]
  specdata$specprop.median[z] = foo.specprop$median[1]
  specdata$specprop.mode[z] = foo.specprop$mode[1]
  specdata$specprop.Q25[z] = foo.specprop$Q25[1]
  specdata$specprop.Q75[z] = foo.specprop$Q75[1]
  specdata$specprop.IQR[z] = foo.specprop$IQR[1]
  specdata$specprop.cent[z] = foo.specprop$cent[1]
  specdata$specprop.skewness[z] = foo.specprop$skewness[1]
  specdata$specprop.kurtosis[z] = foo.specprop$kurtosis[1]
  specdata$specprop.sfm[z] = foo.specprop$sfm[1]
  specdata$specprop.sh[z] = foo.specprop$sh[1]

}
  return(specdata)
}

process_data<-function(whichRun){
#Combine and configure spread detectors. 
l=1
DetecTab<-NULL
if(dettype=="spread"|dettype=="combined"){
  resltsTabspr<-resltsTab[which(resltsTab$detectorType=="spread"),]
  RTSVar<-unique(substr(resltsTabspr$detector,1,3))
  for(d in 1:length(RTSVar)){
    print(paste("Combining spread for",RTSVar[d]))
    resltsTSPVd<-resltsTabspr[which(substr(resltsTabspr$detector,1,3)==RTSVar[d]),]
    for(e in unique(resltsTSPVd$Mooring)){
      resltsTSPV<-resltsTSPVd[which(resltsTSPVd$Mooring==e),]
      print(paste("    ",e))
      for(f in 1:length(unique(resltsTSPV$bottom.freq))){
        resltsTSPV[resltsTSPV$bottom.freq==unique(resltsTSPV$bottom.freq)[f],13]<-f
      }
      colnames(resltsTSPV)[13] <- "detectorRank"
      resltsTSPV$detectorRank<-as.numeric(resltsTSPV$detectorRank)
      resltsTSPV$group[1]<-1
      
      #make sure these values are what they say they are- encountered some behaviors to suggest there was some behind the scenes decimals earlier. This may not be necessary. 
      resltsTSPV$start<-as.numeric(as.character(resltsTSPV$start))
      resltsTSPV$end<-as.numeric(as.character(resltsTSPV$end))
      resltsTSPV$meantime<-(resltsTSPV$start+resltsTSPV$end)/2
      
      #need to order chronologically. 
      resltsTSPV<-resltsTSPV[order(resltsTSPV$meantime,resltsTSPV$bottom.freq),]
      
      #assign groups based on groupInt value
      f<-1
      print("assigning group values")
      for(z in 1:(nrow(resltsTSPV)-1)){
        if(resltsTSPV[z,15]+groupInt[d]>=resltsTSPV[z+1,15]){
          resltsTSPV[z+1,14]<-f
        }else{
          f<-f+1
          resltsTSPV[z+1,14]<-f
        }
      }
      
      #if the meantime is the same take only the lowest # box. 
      #resltsTSPV$remove<-0
      #for(g in 1:(nrow(resltsTSPV)-1)){
      #  if(resltsTSPV[g+1,15]==resltsTSPV[g,15]){
      #    resltsTSPV[g+1,16]<-1}}
      #resltsTSPV<- subset(resltsTSPV,remove==0) #this is working as intended- looks like R truncates the values after 4 digits but does calculate with the full values. 
      #resltsTSPV$remove<-NULL
      
      #remove groups based on grpsize value
      removegrp <- table(resltsTSPV$group)
      resltsTSPV <- subset(resltsTSPV, group %in% names(removegrp[removegrp > (grpsize[d]-1)]))
      #maxgrp<-max(resltsTSPV[,14])
      
      #section for new algorithm, to precede previous RM method. Picks best sequence and subsets data. 
      print(paste("calculating best runs for each group"))
      for(f in unique(resltsTSPV[,14])){
        #print(paste("calculating run for",f,"of",maxgrp))
        groupdat<- subset(resltsTSPV,group==f)
        grpvec<-groupdat[,13]
        colClasses = c("numeric","numeric","numeric","numeric","numeric")
        runsum<- read.csv(text="start, ones, zeros, length,skip", colClasses = colClasses)
        
        for(g in 1:(nrow(groupdat)-(grpsize[d]-1))){
          RM<-groupdat[g,13]
          groupdat[,15+g]<-99
          groupdat[g,15+g]<-2
          skipvec<-0
          for(h in g:(nrow(groupdat)-1)){
            rsltvec0s<-rle(groupdat[,15+g])
            if(any(rsltvec0s$lengths[rsltvec0s$values==0]>allowedZeros[d])){
              break
            }  
            if(RM<grpvec[h+1]&RM+(detskip[d]+1)>grpvec[h+1]&groupdat[h,15]!=groupdat[h+1,15]){
              groupdat[h+1,15+g]<-1
              skipvec<-c(skipvec,(grpvec[h+1]-RM))
              RM<-grpvec[h+1]
            }else if(groupdat[h,15]!=groupdat[h+1,15]){
              groupdat[h+1,15+g]<-0
            }
            if(groupdat[h,15]==groupdat[h+1,15]&groupdat[h,15+g]==0&RM<grpvec[h+1]&RM+(detskip[d]+1)>grpvec[h+1]){
              groupdat[h+1,15+g]<-1
              skipvec<-c(skipvec,(grpvec[h+1]-RM))
              RM<-grpvec[h+1]
            }else if(groupdat[h,15]==groupdat[h+1,15]){
              groupdat[h+1,15+g]<-98
            }
          }
          runsum[g,1]<-g
          runsum[g,2]<-sum(groupdat[,15+g]==1)
          runsum[g,3]<-sum(groupdat[,15+g]==0)
          runsum[g,4]<-(sum(groupdat[,15+g]==1)+sum(groupdat[,15+g]==0)+1)
          runsum[g,5]<-max(skipvec)
        }
        runsum<-runsum[which(runsum[,2]==max(runsum[,2])),] #choose w most ones
        runsum<-runsum[which(runsum[,3]==min(runsum[,3])),] #choose w least 0s
        runsum<-runsum[which(runsum[,5]==min(runsum[,5])),] #choose w smallest maximum skip (most gradual)
        runsum<-runsum[which(runsum[,4]==min(runsum[,4])),] #choose w least length
        runsum<-runsum[1,1] #choose first one
        
        #groupdat<-groupdat[runsum[1,1]:as.numeric((runsum[1,1]+runsum[4]-1)),]
        
        groupdat<-groupdat[,c(1:15,15+runsum)]
        groupdat<-groupdat[which(groupdat[,16]==2|groupdat[,16]==1),]
        groupdat<-groupdat[,c(1:15)]
        
        resltsTSPV<- subset(resltsTSPV,group!=f)
        if(nrow(groupdat)> (grpsize[d]-1)){ #only rbind if group is above minimum groupdat size after processing 
        resltsTSPV<-rbind(resltsTSPV,groupdat)
        }
      }
      
      if(nrow(resltsTSPV)==0){
        write.table("FINAL There were no detections",paste(outputpath,runname,"/",e,"/FINAL_Summary_spread_",substr(resltsTSPVd$detector[1],1,3),"_",length(detectorsspr[[d]]),"dnum_","_",d,".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE,col.names=FALSE)
      }else{
        
        colClasses = c("numeric", "character","numeric","numeric", "numeric","numeric","numeric","numeric","character","character", "numeric","character")
        resltsTSPVFinal <- read.csv(text="Selection,View,Channel,Begin Time (s),End Time (s),Low Freq (Hz),High Freq (Hz), DetectorCount, DetectorName, DetectorType, numDetectors, Mooring", colClasses = colClasses)
        colnames(resltsTSPVFinal)<-c("Selection","View","Channel","Begin Time (s)","End Time (s)","Low Freq (Hz)","High Freq (Hz)","DetectorCount", "DetectorName", "DetectorType","numDetectors","Mooring")
        
        
        p=1
        for(j in unique(resltsTSPV$group)){
          grpminfreq <- min(resltsTSPV[resltsTSPV$group==j,6])
          grpmaxfreq <- max(resltsTSPV[resltsTSPV$group==j,7])
          grpstarttime <- min(resltsTSPV[resltsTSPV$group==j,4])
          grpendtime <- max(resltsTSPV[resltsTSPV$group==j,5])
          
          resltsTSPVFinal[p,1]<-p
          resltsTSPVFinal[p,4]<-grpstarttime
          resltsTSPVFinal[p,5]<-grpendtime
          resltsTSPVFinal[p,6]<-grpminfreq
          resltsTSPVFinal[p,7]<-grpmaxfreq
          resltsTSPVFinal[p,8]<-l
          resltsTSPVFinal[p,9]<-substr(resltsTSPV$detector[1],1,3)
          resltsTSPVFinal[p,10]<-"spread"
          resltsTSPVFinal[p,11]<-length(detectorsspr[[d]])
          resltsTSPVFinal[p,12]<-e
          
          p<-p+1
        }
        
        resltsTSPVFinal$View<-"Spectrogram 1"
        resltsTSPVFinal$Channel<-1
        
        DetecTab<- rbind(DetecTab,resltsTSPVFinal)
        
        resltsTSPVFinal<- resltsTSPVFinal[,1:7]
        write.table(resltsTSPVFinal,paste(outputpath,runname,"/",e,"FINAL_Summary_spread_",substr(resltsTSPVd$detector[1],1,3),"_",length(detectorsspr[[d]]),"dnum_","_",d,".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
      }
    }
    l<-l+1}
}


if(dettype=="combined"){
}else{
  l=1
}

#Configure single detectors. 
if(dettype=="single"|dettype=="combined"){
  resltsTabsin<-resltsTab[which(resltsTab$detectorType=="single"),]
  for(d in 1:length(unique(resltsTabsin$detector))){
    print(paste("Configuring single detector",unique(resltsTabsin$detector)[d]))
    resltsTSGVd<-resltsTabsin[which(resltsTabsin$detector==unique(resltsTabsin$detector)[d]),]
    for(e in unique(resltsTSGVd$Mooring)){
      resltsTSGV<-resltsTSGVd[which(resltsTSGVd$Mooring==e),]
      print(paste("     ",e))
      if(nrow(resltsTSGV)==0){
        write.table("FINAL There were no detections",paste(outputpath,runname,"/",e,"/FINAL_Summary_single_",resltsTSGVd$detector[1],"_","_",d,".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE,col.names=FALSE)
      }else{
        
        colClasses = c("numeric", "character","numeric","numeric", "numeric","numeric","numeric","numeric","character","character", "numeric","character")
        resltsTSGVFinal <- read.csv(text="Selection,View,Channel,Begin Time (s),End Time (s),Low Freq (Hz),High Freq (Hz), DetectorCount, DetectorName, DetectorType, numDetectors, Mooring", colClasses = colClasses)[1:nrow(resltsTSGV), ]
        colnames(resltsTSGVFinal)<-c("Selection","View","Channel","Begin Time (s)","End Time (s)","Low Freq (Hz)","High Freq (Hz)","DetectorCount", "DetectorName", "DetectorType","numDetectors","Mooring")
        
        resltsTSGVFinal$Selection<-seq(1:nrow(resltsTSGV))
        resltsTSGVFinal$View<-"Spectrogram 1"
        resltsTSGVFinal$Channel<-1
        resltsTSGVFinal[,4]<-resltsTSGV$start
        resltsTSGVFinal[,5]<-resltsTSGV$end
        resltsTSGVFinal[,6]<-resltsTSGV$bottom.freq
        resltsTSGVFinal[,7]<-resltsTSGV$top.freq
        resltsTSGVFinal$DetectorCount<-l
        resltsTSGVFinal$DetectorName<- resltsTSGV$detector[1]
        resltsTSGVFinal$DetectorType<-"single"
        resltsTSGVFinal$numDetectors<-1
        resltsTSGVFinal$Mooring<-e
        
        DetecTab<- rbind(DetecTab,resltsTSGVFinal)
        
        resltsTSGVFinal<- resltsTSGVFinal[,1:7]
        write.table(resltsTSGVFinal,paste(outputpath,runname,"/",e,"FINAL_Summary_single_",resltsTSGVd$detector[1],"_",d,".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
      }
    }
    l=l+1}
}

############

#now need to average detections between detectors. set 
DetecTab$meantime<-(DetecTab[,4]+DetecTab[,5])/2
#need to order chronologically
DetecTab<-DetecTab[order(DetecTab$meantime),]
DetecTab$meanfreq<-(DetecTab[,6]+DetecTab[,7])/2
DetecTab$UniqueID<-seq(1:nrow(DetecTab))
DetecTab$remove<-0

#average detections within detectors using timediffself parameter
n=0
for(o in unique(DetecTab$Mooring)){
  Tab<-DetecTab[which(DetecTab$Mooring==o),]
  print(paste("For mooring",o))
  for(p in 1:max(Tab$DetectorCount)){
    r=0
    AvgTab<- Tab[which(Tab$DetectorCount==p),]
    print(paste("       Average within detector",AvgTab$DetectorName[1]))
    newrow<-AvgTab[0,]
    IDvec<-NULL
    if(nrow(AvgTab)>1){
    for(q in 1:(nrow(AvgTab)-1)){
      if(AvgTab[q+1,13]<=(AvgTab[q,13]+timediffself)){
        
        newdat<-AvgTab[0,]
        newdat[1:2,]<-AvgTab[c(q,q+1),]
        IDvec<-c(IDvec,newdat$UniqueID)
        s<-mean(newdat[,4])
        e<-mean(newdat[,5])
        h<-mean(newdat[,6])
        l<-mean(newdat[,7])
        mt<-(s+e)/2
        mf<-(h+l)/2
        unqID<-newdat[1,15]
        newrow[r+1,]<-data.frame(newdat[1,1],newdat[1,2],newdat[1,3],s,e,h,l,newdat[1,8],newdat[1,9],newdat[1,10],newdat[1,11],newdat[1,12],mt,mf,unqID,0, stringsAsFactors = FALSE)
        r=r+1
      }
    }
    }
    if(r>0){
      AvgTab<-subset(AvgTab,!(AvgTab$UniqueID %in% IDvec))
      AvgTab<-rbind(AvgTab,newrow)
      AvgTab<-AvgTab[order(AvgTab$meantime),]
      Tab<- Tab[-which(Tab$DetectorCount==p),]
      Tab<-rbind(Tab,AvgTab)
      n=n+1
    }
    
  }
  if(n>0){
    DetecTab<-DetecTab[-which(DetecTab$Mooring==o),]
    DetecTab<-rbind(DetecTab,Tab)
  }
}


DetecTab<-DetecTab[order(DetecTab$meantime),]

#average between detectors
AvgDet<-DetecTab[0,]
DetecTab2<-DetecTab[0,]
for(w in unique(DetecTab$Mooring)){
  print(paste("For mooring",w))
  CompareDet<-DetecTab[which(DetecTab$Mooring==w),]
  if(detnum>1){
    CDvar<-unique(CompareDet$DetectorCount)
    for(x in 1:(length(CDvar)-1)){
      i = CompareDet[which(CompareDet$DetectorCount==CDvar[x]),]
      j = CompareDet[which(CompareDet$DetectorCount==CDvar[x+1]),]
      print(paste("      Average detectors",i[1,9],"and",j[1,9]))
      for(y in 1:nrow(i)){
        jvec <- which(j$meantime<(i$meantime+2)&j$meantime>(i$meantime-2))
        if(length(jvec)>0){
        for(z in min(jvec):max(jvec)){
          if(((((i[y,13]-j[z,13])<=timediff) & (i[y,13]>=j[z,13])) | (((j[z,13]-i[y,13])<=timediff) & (j[z,13]>=i[y,13]))) & ((((i[y,14]-j[z,14])<=freqdiff) & (i[y,14]>=j[z,14])) | (((j[z,14]-i[y,14])<=freqdiff) & (j[z,14]>=i[y,14])))){
            CompareDet[which(i[y,15]==CompareDet$UniqueID),16]<-1
            CompareDet[which(j[z,15]==CompareDet$UniqueID),16]<-1
            meanS<-mean(c(i[y,4],j[z,4]))
            meanE<-mean(c(i[y,5],j[z,5]))
            meanL<-mean(c(i[y,6],j[z,6]))
            meanH<-mean(c(i[y,7],j[z,7]))
            
            AvgDet2<-data.frame(99,as.character("Spectrogram 1"),1,meanS,meanE,meanL,meanH,j[1,8],as.character(paste(i[1,9],"+",j[1,9])),as.character(paste(i[1,10],"+",j[1,10])),(j[1,11]+i[1,11]),as.character(w),mean(c(meanS,meanE)),mean(c(meanL,meanH)),as.integer(99),0)
            names(AvgDet2)<-colnames(AvgDet)
            #make new dataframe be second detector ID so it will loop properly
            
            AvgDet<- rbind(AvgDet,AvgDet2)  
            
          }
        }
        }
      }
      CompareDet<-CompareDet[-which(CompareDet$remove==1),]
      
      CompareDet[which(CompareDet$DetectorCount==i[1,8]),8]<-j[1,8]
      CompareDet[which(CompareDet$DetectorName==i[1,9] |CompareDet$DetectorName==j[1,9]),9]<-as.character(paste(i[1,9],"+",j[1,9]))
      CompareDet[which(CompareDet$DetectorType==i[1,10] |CompareDet$DetectorType==j[1,10]),10]<-as.character(paste(i[1,10],"+",j[1,10]))
      CompareDet[which(CompareDet$numDetector==i[1,11] |CompareDet$numDetector==j[1,11]),11]<-(j[1,11]+i[1,11])
      
      CompareDet<-rbind(CompareDet,AvgDet)
      CompareDet<-CompareDet[order(CompareDet$meantime),]
      
      CompareDet$UniqueID<-seq(1:nrow(CompareDet))
      AvgDet<-DetecTab[0,]
    }
  }
  DetecTab2<-rbind(DetecTab2,CompareDet)
}

DetecTab2<-DetecTab2[order(DetecTab2$meantime),]

#average detections within combined detector using timediffself parameter (3x)
for(a in 1:3){
n=0
for(o in unique(DetecTab2$Mooring)){
  Tab<-DetecTab2[which(DetecTab2$Mooring==o),]
  print(paste("For mooring",o))   
  for(p in unique(Tab$DetectorCount)){
    r=0
    AvgTab<- Tab[which(Tab$DetectorCount==p),]
    print(paste("       Average combined detector for time number",a))
    
    newrow<-AvgTab[0,]
    IDvec<-NULL
    if(nrow(AvgTab)>1){
    for(q in 1:(nrow(AvgTab)-1)){
      if(AvgTab[q+1,13]<=(AvgTab[q,13]+timediffself)){
        
        newdat<-AvgTab[0,]
        newdat[1:2,]<-AvgTab[c(q,q+1),]
        IDvec<-c(IDvec,newdat$UniqueID)
        s<-mean(newdat[,4])
        e<-mean(newdat[,5])
        h<-mean(newdat[,6])
        l<-mean(newdat[,7])
        mt<-(s+e)/2
        mf<-(h+l)/2
        unqID<-newdat[1,15]
        newrow[r+1,]<-data.frame(newdat[1,1],newdat[1,2],newdat[1,3],s,e,h,l,newdat[1,8],newdat[1,9],newdat[1,10],newdat[1,11],newdat[1,12],mt,mf,unqID,0, stringsAsFactors = FALSE)
        r=r+1
      }
    }
    }
    if(r>0){
      AvgTab<-subset(AvgTab,!(AvgTab$UniqueID %in% IDvec))
      AvgTab<-rbind(AvgTab,newrow)
      AvgTab<-AvgTab[order(AvgTab$meantime),]
      Tab<- Tab[-which(Tab$DetectorCount==p),]
      Tab<-rbind(Tab,AvgTab)
      n=n+1
    }
    
  }
  if(n>0){
    DetecTab2<-DetecTab2[-which(DetecTab2$Mooring==o),]
    DetecTab2<-rbind(DetecTab2,Tab)
  }
}
}

#remove detections that do not fit min/max duration parameters 
DetecTab2$UniqueID<-NULL
for(d in unique(DetecTab2$Mooring)){
  Tab<-DetecTab2[which(DetecTab2$Mooring==d),]
  for(f in 1:nrow(Tab)){
    if((Tab[f,5]-Tab[f,4])>Maxdur|(Tab[f,5]-Tab[f,4])<Mindur){
      Tab[f,15]<-1
    }
  }
  DetecTab2<-DetecTab2[-which(DetecTab2$Mooring==d),]
  DetecTab2<-rbind(DetecTab2,Tab)
}

DetecTab2<-DetecTab2[which(DetecTab2$remove==0),]

DetecTab2$remove<-NULL

#add information on original sound files and calculate time since file start
if(whichRun==2){
DetecTab2$File<-""
DetecTab2$FileStartSec<-0
DetecTab2$FileOffsetBegin<-0
DetecTab2$FileOffsetEnd<-0
for(w in unique(DetecTab2$Mooring)){
  print(paste("calculate file ID and begin time and end time relative to file for mooring",w))
  DetecVar<-DetecTab2[which(DetecTab2$Mooring==w),]
  durVar<-durTab[which(durTab$CombSF==w),]
  for(c in 1:nrow(DetecVar)){
    DetecVar$File[c]<-as.character(durVar[findInterval(DetecVar[c,13], c(0,durVar$CumDur)),2])
    DetecVar$FileStartSec[c]<-durVar[findInterval(DetecVar[c,13], c(0,durVar$CumDur)),4]
  }
  DetecVar$FileStartSec<-DetecVar$FileStartSec-DetecVar$FileStartSec[1]
  DetecVar$FileOffsetBegin<-DetecVar$`Begin Time (s)`-DetecVar$FileStartSec
  DetecVar$FileOffsetEnd<-DetecVar$`End Time (s)`-DetecVar$FileStartSec  
  DetecTab2<-DetecTab2[-which(DetecTab2$Mooring==w),]
  DetecTab2<-rbind(DetecTab2,DetecVar)
  
}

DetecTab2$FileStartSec<-NULL
}
DetecTab2$Selection<-seq(1,nrow(DetecTab2))

return(DetecTab2)
}

#paths
startcombpath<-"E:/Combined_sound_files/"
BLEDpath<-"C:/Users/danby456/Raven Pro 1.5/Presets/Detector/Band Limited Energy Detector/"
ravenpath<-"C:/Users/danby456/Raven Pro 1.5"
outputpath<-"E:/DetectorRunOutput/"

#moorings completed 
allmooringsGT<- c("BS15_AU_02a","BS14_AU_04","AW12_AU_BS3","BS13_AU_04","BS16_AU_02a","BS15_AU_02b","AW14_AU_BS3") #add as complete GTs 
allmooringsSF<-list()#list sound file range for comleted GT of each mooring 
allmooringsSF[[1]]<-c(1,104)
allmooringsSF[[2]]<-c(1,179)
allmooringsSF[[3]]<-c(1,217)
allmooringsSF[[4]]<-c(1,304)
allmooringsSF[[5]]<-c(1,175)
allmooringsSF[[6]]<-c(1,62)
allmooringsSF[[7]]<-c(1,160)

#############################

MooringsDat<-rbind(allmooringsGT,matrix(unlist(allmooringsSF), nrow=length(unlist(allmooringsSF[1]))))
colnames(MooringsDat)<-c(allmooringsGT)
MooringsDat<-MooringsDat[,order(colnames(MooringsDat))] 

################Script function

#enter the run name:
runname<- "Decent detector test"

#Run type: all (all) or specific (spf) moorings to run
runtype<-"all"

#enter the detector type: "spread" or "single" or "combined". Can run and combine any combination of spread and single detectors that will be averaged after returning their detections. 
dettype<- "spread" 

#enter the type of mooring you'd like to analyze data: high graded (HG) or on full mooring (FULL)
moorType<-"FULL"

#Enter the name of the species you'd like to evaluate (RW,GS):
spec <- "RW"

#compare detections with pulses and fin/mooring noise, other sources of intereference y or n
interfere<-"n"

interfereVec<-c(dir(BLEDpath)[6])



if(dettype=="spread"|dettype=="combined"){
#make a list of detectors you wish to run. Must correspond with those of same name already in BLED folder in Raven. 
detectorsspr<-list()
detectorsspr[[1]] <- dir(BLEDpath)[20:37] #add more spreads with notation detectorspr[[x]]<-... #15-32
#detectorsspr[[2]] <- dir(BLEDpath)[3:14]
detectorssprshort<- detectorsspr
}

if(dettype=="single"|dettype=="combined"){
detectorssin <- c(dir(BLEDpath)[1]
                  ) #list single detectors to run 
detectorssinshort<- detectorssin
}
##################Sampling rate of sound files

samplingRate<-16384 #hz

##########################max min length parameters (applies on R final detections, can also change in Raven to change initial box size)

Maxdur<-3.5
Mindur<-0.2

############################Combine detector  parameters
timediffself<-1

#multiple detectors
freqdiff<-100
timediff<-1


############################Whiten parameters (need to have done this in Raven previously)

#Pre whiten data?(y or no)
whiten<-"y"
FO<-100 #filter order
LMS<-.10 #LMS step size

############################Spread parameters. must be same length as number of spread detectors you are running
#p7 good ones: 2,1,3,0.75
#p9 working ones: 3,2,3,.25
#p10 good ones: 3,2,4,0.5
#(SPREAD) enter the desired smallest sequence size for detection. 
grpsize<-c(4)

#(SPREAD) allowed consecutive descending boxes allowed to still constitute an ascending sequence. Will end sequence after the maximum has been exceeded
allowedZeros<-c(2)

#(SPREAD) threshold of how many detectors at most can be skipped to be counted as sequential increase. 
detskip<-c(5)

#(SPREAD) max time distance for detectors to be considered in like group 
groupInt<-c(0.35)

############################
runname<-paste(runname,gsub("\\D","",Sys.time()),sep="_")
dir.create(paste(outputpath,runname,sep=""))

if(runtype=="all"){
moorings<- colnames(MooringsDat)
#SF<-allmooringsSF
}else{
  allmooringsGT<- c("AW12_AU_BS3") #add as complete GTs 
  allmooringsSF<-list()#list sound file range for comleted GT of each mooring 
  allmooringsSF[[1]]<-c(1,217)
 # allmooringsSF[[2]]<-c(1,96)
  
  MooringsDat<-rbind(allmooringsGT,matrix(unlist(allmooringsSF), nrow=length(unlist(allmooringsSF[1]))))
  colnames(MooringsDat)<-c(allmooringsGT)
  if(ncol(MooringsDat)>1){
  MooringsDat<-MooringsDat[,order(colnames(MooringsDat))] 
  }else{
  }
  moorings<-colnames(MooringsDat)
}


detlist<-NULL
detlist2<-NULL
if(dettype=="spread"|dettype=="combined"){
  for(n in 1:length(detectorsspr)){
    detlist<-c(detlist,length(detectorsspr[[n]]))
    sonlydetlist<-detlist
    detlist2<-c(detlist2,substr(detectorssprshort[[n]][1],1,3))
    sonlydetlist2<-detlist2
  }
}

if(dettype=="single"|dettype=="combined"){
    detlist<- c(detlist,length(detectorssin))
    detlist2<-c(detlist2,detectorssinshort)
}

if(dettype=="single"){
  sonlydetlist<-0
  sonlydetlist2<-NA
}

detnum<-if(dettype=="spread"){
    length(detectorsspr)
  }else if(dettype=="single"){
      length(detectorssin)
    }else if(dettype=="combined"){
        length(detectorssin)+length(detectorsspr)
    }

###########################make txt file of params for run:
colClasses = c("numeric", "character","character")
ParamSum <- read.csv(text="Parameter,Value,Description",colClasses = colClasses)
colnames(ParamSum)<-c("Parameter","Value","Description")

ParamSum[1,1]<-"Run name:"
ParamSum[1,2]<- runname
ParamSum[1,3]<-""
ParamSum[2,1]<-"Run type:"
ParamSum[2,2]<- runtype
ParamSum[2,3]<-"All available moorings ('all') or specific moorings ('spf')"
ParamSum[3,1]<-"Detector type:"
ParamSum[3,2]<- dettype
ParamSum[3,3]<-"Type of detector(s) used in run. spread, single, or combined (both)"
ParamSum[3,1]<-"Number/name of detectors ran:"
ParamSum[3,2]<- paste(detnum,paste(detlist2,collapse=" "))
ParamSum[3,3]<-"Number of detectors used in run"
ParamSum[4,1]<-"Species:"
ParamSum[4,2]<- spec
ParamSum[4,3]<-"Call type that detector(s) will look for (RW=right whale,GS=gunshot etc.)"
ParamSum[5,1]<-"Time threshold for combination:"
ParamSum[5,2]<- paste(paste(timediff,"s",sep=""),paste(timediffself,"s",sep=""),sep=",")
ParamSum[5,3]<-"Maximum second difference between detection mean time to be considered a combined detection. Self and between detectors"
ParamSum[6,1]<-"Frequency threshold for combination:"
ParamSum[6,2]<- paste(freqdiff,"Hz",sep="")
ParamSum[6,3]<-"Maximum Hz difference between detection mean freq to be considered a combined detection"
ParamSum[7,1]<-"Min/Max duration for final detections"
ParamSum[7,2]<- paste(Mindur,Maxdur)
ParamSum[7,3]<- " "
ParamSum[8,1]<-"Number/name spread detectors ran:"
ParamSum[8,2]<- paste(length(sonlydetlist),paste(sonlydetlist2,collapse="/"))
ParamSum[8,3]<- " "


if(dettype=="spread"|dettype=="combined"){
  
#different table for spread parameters
colClasses = c("numeric", "character","character")
ParamSum2 <- read.csv(text="Parameter,Value,Description",colClasses = colClasses)
colnames(ParamSum2)<-c("Parameter","Value","Description")

ParamSum2[1,1]<-"Spread group size:"
ParamSum2[1,2]<- paste(grpsize,collapse = ",")
ParamSum2[1,3]<-"Minimum amount of detections needed to be considered a detection group"
ParamSum2[2,1]<-"Zeros:"
ParamSum2[2,2]<- paste(allowedZeros,collapse = ",")
ParamSum2[2,3]<-"Maximum number of consecutive descending boxes allowed in an ascending sequence"
ParamSum2[3,1]<-"Skip threshold:"
ParamSum2[3,2]<- paste(detskip,collapse = ",")
ParamSum2[3,3]<-"Maximum amount of skips between detector rank to count towards ascending sequence distinction"
ParamSum2[4,1]<-"Time threshold for group:"
ParamSum2[4,2]<- paste(groupInt,collapse = ",")
ParamSum2[4,3]<-"Maximum second difference between detection mean time to be considered a detection group"
ParamSum2[5,1]<-"Detectors ran in each spread:"
ParamSum2[5,2]<- paste(detlist,collapse = ",")
ParamSum2[5,3]<- " "

  ParamSum<-rbind(ParamSum,ParamSum2)
}


if(whiten=="y"){

colClasses = c("numeric", "character","character")
ParamSum3 <- read.csv(text="Parameter,Value,Description",colClasses = colClasses)
colnames(ParamSum3)<-c("Parameter","Value","Description")

ParamSum3[1,1]<-"Whiten Filter Order (FO)"
ParamSum3[1,2]<- FO
ParamSum3[1,3]<-" "
ParamSum3[2,1]<-"Whiten LMS step size (x10^-9)" 
ParamSum3[2,2]<- LMS
ParamSum3[2,3]<-" " 

ParamSum<-rbind(ParamSum,ParamSum3)

}


write.table(ParamSum,paste(outputpath,runname,"/","Params_",dettype,"_",runname,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE)

######################################
#set FO and LMS to NA if no whiten
if(whiten=="n"){
  FO<-NA
  LMS<-NA
}

#path to ground truth table
GT<-list()
for(f in 1:length(moorings)){
  GT[[f]] <- read.delim(paste("E:/Selection tables/",moorings[f],"Sum/",moorings[f],"_All.txt",sep=""))
  GT[[f]] <- GT[[f]][GT[[f]]$View=="Spectrogram 1",]
}

if(dettype=="spread"|dettype=="combined"){
for(i in 1:length(detectorsspr)){
  for(j in 1:length(detectorsspr[[i]])){
    detectorsspr[[i]][j]<-paste(BLEDpath,detectorsspr[[i]][j],sep="")
  }
}
}
if(dettype=="single"|dettype=="combined"){
  for(k in 1:length(detectorssin)){
    detectorssin[k]<-paste(BLEDpath,detectorssin[k],sep="")
  }  
}

#run sound files:
resltsTab <- NULL
resltsTabInt<- NULL
durTab<-NULL
for(m in moorings){
  if(whiten=="n"){
  whiten2<-"No_whiten"
  sound_files <- dir(paste("E:/Datasets/",m,"/",spec,"_ONLY_yesUnion",sep = ""))[MooringsDat[2,colnames(MooringsDat)==m]:MooringsDat[3,colnames(MooringsDat)==m]] #based on amount analyzed in GT set
  sound_filesfullpath <- paste("E:/Datasets/",m,"/",spec,"_ONLY_yesUnion/",sound_files,sep = "")
  #too ineffecient to run sound files one by one, so check to see if combined file exists and if not combine them. 
  combSound<-paste(startcombpath,spec,"/",whiten2,"/",m,"_files",MooringsDat[2,colnames(MooringsDat)==m],"-",MooringsDat[3,colnames(MooringsDat)==m],".wav",sep="")
  if(file.exists(combSound)){
  }else{
    dir.create(paste(startcombpath,spec,sep=""))
    dir.create(paste(startcombpath,spec,"/",whiten2,"/",sep=""))
    sox_alt(paste(noquote(paste(paste(sound_filesfullpath[MooringsDat[2,colnames(MooringsDat)==m]:MooringsDat[3,colnames(MooringsDat)==m]],collapse=" ")," ",combSound,sep=""))),exename="sox.exe",path2exe="E:\\Accessory\\sox-14-4-2")
  }
  
  
  }else{
  whiten2 <- paste("/Bbandp",100*LMS,"x_","FO",FO,"/",sep = "")
  }

  combname<- paste(m,"_files",MooringsDat[2,colnames(MooringsDat)==m],"-",MooringsDat[3,colnames(MooringsDat)==m],".wav",sep="")
  
#run pulse and fin/mooring detector, if selected:
if(interfere=="y"){
  for(i in interfereVec){
    print(paste("Running detector for",m))
    resltVarInt <- raven_batch_detec(raven.path = ravenpath, sound.files = combname, path = paste(startcombpath,spec,"/",whiten2,sep=""),detector = "Band Limited Energy Detector",dpreset=i,vpreset="RW_Upcalls")
    resltVarInt$Mooring<-m
    resltVarInt$detector<-i
    resltVarInt$detectorType<-"intereference"
    resltVarInt$detectorCount<-which(interfereVec==i)
    if(is.null(nrow(resltVarInt))==FALSE){
    resltsTabInt<-rbind(resltsTabInt,resltVarInt)
    }
    resltVarInt<-NULL
  }
  
}
#run detector(s)
if(dettype=="spread"|dettype=="combined"){
  for(q in 1:length(detectorssprshort)){
    for(r in detectorssprshort[[q]]){
      print(paste("Running detector for",m))
      resltVar <- raven_batch_detec(raven.path = ravenpath, sound.files = combname, path = paste(startcombpath,spec,"/",whiten2,sep=""),detector = "Band Limited Energy Detector",dpreset=r,vpreset="RW_Upcalls")
      resltVar$Mooring<-m
      resltVar$detector<-r
      resltVar$detectorType<-"spread"
      resltVar$detectorCount<-q
      if(is.null(nrow(resltVar))==FALSE){
      resltsTab<- rbind(resltsTab,resltVar)
      }
      resltVar<-NULL 
    }
    }
  }
 
if(dettype=="single"|dettype=="combined"){
  for(n in detectorssinshort){
    print(paste("Running detector for",m))
    resltVar <- raven_batch_detec(raven.path = ravenpath, sound.files =  combname, path = paste(startcombpath,spec,"/",whiten2,sep=""),detector = "Band Limited Energy Detector",dpreset=n,vpreset ="RW_Upcalls")
    resltVar$Mooring<-m
    resltVar$detector<-n
    resltVar$detectorType<-"single"
    resltVar$detectorCount<-which(detectorssinshort==n)
    if(is.null(nrow(resltVar))==FALSE){
    resltsTab<- rbind(resltsTab,resltVar)
    }
    resltVar<-NULL
    }
  }
}

#Combine and configure spread detectors. 
DetecTab2<-process_data(1)

DetecTab2$detectionType<-0

#Define table for later excel file export. 
colClasses = c("character","character","character","character","character","numeric","numeric","numeric", "numeric","numeric","numeric","character","character","character","character","character","character","character","character","numeric","numeric","character")
detecEvalFinal <- read.csv(text="Species, Moorings, Detectors, DetType, RunName, numTP, numFP, numFN, TPhitRate, TPR, TPdivFP, ZerosAllowed,GroupSize,SkipAllowance,GroupInterval,TimeDiff,TimeDiffself,MinMaxDur,numDetectors,FO,LMS,Notes", colClasses = colClasses)

  GTtot=0
  GTtot2=0
  TPtot=0
  MoorCor=0
for(v in 1:length(unique(DetecTab2$Mooring))){
  print(paste("Comparing ground truth of",sort(unique(DetecTab2$Mooring))[v],"with final detector"))   
  MoorVar<-DetecTab2[which(DetecTab2$Mooring==sort(unique(DetecTab2$Mooring))[v]),]
  MoorVar$Selection<-seq(1:nrow(MoorVar))
  #this table is mostly buggy and useless in its stats 
  write.csv(MoorVar,paste(outputpath,runname,"/",MoorVar[1,12],"_Summary_",dettype,"_Info",".csv",sep=""),quote=FALSE,row.names=FALSE)
  #useable table to evaluate just results of combined detectors. 
  write.table(MoorVar[,1:7],paste(outputpath,runname,"/",MoorVar[1,12],"_Summary_",dettype,"_Ravenformat",".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
  
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
    gvec <- which(GT[[v]]$meantime<(MoorVar$meantime[h]+2)&GT[[v]]$meantime>(MoorVar$meantime[h]-2))
    if(length(gvec)>0){
    for(g in min(gvec):max(gvec)){
      if((MoorVar[h,13]>GT[[v]][g,4]) & (MoorVar[h,13]<GT[[v]][g,5])){
        OutputCompare[p,]<-MoorVar[h,c(1:7,13,15)]
        OutputCompare[p,9]<-"TP"
        p=p+1
      }
    }
  }
  }
  #Identify and add FPs. if selection in MoorVar row does not match that in Output compare, add it to Output compare under designation FP.  
  if(nrow(OutputCompare)>0){
  OutputCompare <- rbind(OutputCompare,MoorVar[-which(MoorVar$Selection %in% OutputCompare$Selection),c(1:7,13,15)])
  OutputCompare[which(OutputCompare$detectionType!="TP"),9]<-"FP"
  }
  
  #Add rows where GT meantime was in between 
  p=1
  for(h in 1:nrow(GT[[v]])){
    gvec <- which(MoorVar$meantime<(GT[[v]]$meantime[h]+2)&MoorVar$meantime>(GT[[v]]$meantime[h]-2))
    if(length(gvec)>0){
    for(g in min(gvec):max(gvec)){
      if(GT[[v]][h,8]>MoorVar[g,4] & GT[[v]][h,8]<MoorVar[g,5]){
        OutputCompare2[p,]<-GT[[v]][h,]
        OutputCompare2[p,9]<-"TP truth"
        p=p+1
      }
    }
  }
  }
  
  #Identify and add FNs. if selection in GT row does not match that in OutputCompare2, add it to Output compare under designation FN.  
  if(nrow(OutputCompare2)>0){
  OutputCompare2 <- rbind(OutputCompare2,GT[[v]][-which(GT[[v]]$Selection %in% OutputCompare2$Selection),])
  OutputCompare2[which(OutputCompare2$detectionType!="TP truth"),9]<-"FN"
  
  #Combine tables and remove GT TPs from dataset. 
  OutputCompare<-rbind(OutputCompare,OutputCompare2)
  OutputCompare$meantime<-as.numeric(OutputCompare$meantime)
  OutputCompare<-OutputCompare[order(OutputCompare$meantime),]
  OutputCompare[which(OutputCompare$detectionType=="TP truth"),9]<-"x"
  OutputCompare <- subset(OutputCompare,detectionType!="x")
  OutputCompare$Selection<-seq(1:nrow(OutputCompare))
  
  #compare detections to sources of interference
  if(interfere=="y"){
    for(n in 1:max(resltsTabInt$detectorCount)){
      MoorInt<-resltsTabInt[which(resltsTabInt$Mooring==sort(unique(DetecTab2$Mooring))[v]),]
      OutputCompare[,n+9]<-0
      colnames(OutputCompare)[length(OutputCompare)]<-paste(resltsTabInt[which(resltsTabInt$detectorCount==n),10])[1]
      MoorInt<-MoorInt[which(MoorInt$detectorCount==n),]
      print(paste("       Compare with detector",MoorInt[1,10]))   
      for(g in 1:nrow(OutputCompare)){
        hvec <- which(MoorInt$meantime<(OutputCompare$meantime[g]+2)&MoorInt$meantime>(OutputCompare$meantime[g]-2))
        if(length(hvec)>0){
        for(h in min(hvec):max(hvec)){
          if((MoorInt[h,4]<OutputCompare[g,8] & MoorInt[h,5]>OutputCompare[g,8])|(MoorInt[h,4]>OutputCompare[g,4] & MoorInt[h,5]<OutputCompare[g,5])){
            OutputCompare[g,n+9]<-1
          }
        }
        }
        }
      }
    }
  
  #Ready table for Raven and save. 
  colnames(OutputCompare)[8]<-'TP/FP/FN'
  OutputCompareRav<- OutputCompare[,-8]
  OutputCompareRav<- OutputCompareRav[,1:8]
  
  OutputCompareRav$Selection<-seq(1:nrow(OutputCompareRav))
  write.table(OutputCompareRav,paste(outputpath,runname,"/",MoorVar[1,12],OutputCompare$Mooring[1],"_TPFPFN_Tab_Ravenformat.txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)

  OutputCompareRF<- OutputCompare[,-8]
  OutputCompareRF$Selection<-seq(1:nrow(OutputCompareRF))
  write.table(OutputCompareRF,paste(outputpath,runname,"/",MoorVar[1,12],OutputCompare$Mooring[1],"_TPFPFN_Tab_RF.txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
  
  OutputCompare<- OutputCompare[,-8]
  OutputCompare<- OutputCompare[,1:8]

  }else{
    write.table("There were no true positive detections for this mooring",paste(outputpath,runname,"/",MoorVar[1,12],OutputCompare$Mooring[1],"_TPFPFN_Tab_Ravenformat.txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE,col.names=FALSE)
    
  }
  
  
  #Make summary table of statistics for table comparison. 
  numTP <- nrow(OutputCompare[which(OutputCompare[,8]=="TP"),])
  numFN <- nrow(OutputCompare[which(OutputCompare[,8]=="FN"),])
  numFP <- nrow(OutputCompare[which(OutputCompare[,8]=="FP"),])
  numTPtruth<- nrow(GT[[v]])

  #store this for comparison with full mooring later
  TPtot<-c(TPtot,numTP)
  MoorCor<-c(MoorCor,sort(unique(DetecTab2$Mooring))[v])
  GTtot2<-c(GTtot2,GTtot)
  
  TPhitRate <- numTP/numTPtruth*100
  TPR <- numTP/(numTP+numFN)
  TPdivFP<- numTP/numFP

  #save stats and parameters to excel file
  detecEval<-detecEvalFinal[0,]
  if(dettype=="spread"|dettype=="combined"){
    detecEval[1,]<-c(spec,sort(unique(DetecTab2$Mooring))[v],paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,numTP,numFP,numFN,TPhitRate,TPR,TPdivFP,paste(allowedZeros,collapse=","),paste(grpsize,collapse=","),paste(detskip,collapse=","),paste(groupInt,collapse=","),timediff,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")
  }else{
    detecEval[1,]<-c(spec,sort(unique(DetecTab2$Mooring))[v],paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,numTP,numFP,numFN,TPhitRate,TPR,TPdivFP,NA,NA,NA,NA,timediff,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")   
  }
  detecEvalFinal <- rbind(detecEvalFinal,detecEval)
  
  png(paste(outputpath,runname,"/",MoorVar[1,12],"_Distribution.png",sep =""), width = 600, height = 300)
  hist(MoorVar$meantime,main=paste(MoorVar[1,12],"Detections"))
  dev.off()
  
  png(paste(outputpath,runname,"/",MoorVar[1,12],"_GTDistribution.png",sep =""), width = 600, height = 300)
  hist(GT[[v]]$meantime,main=paste(MoorVar[1,12],"Ground Truth Detections"))
  dev.off()
}

#Make summary table of whole run statistics for table comparison. 
numTP <- sum(as.numeric(detecEvalFinal[,6]))
numFN <- sum(as.numeric(detecEvalFinal[,8]))
numFP <- sum(as.numeric(detecEvalFinal[,7]))
numTPtruth<- GTtot

TPhitRate <- numTP/numTPtruth*100
TPR <- numTP/(numTP+numFN)
TPdivFP<- numTP/numFP

#save stats and parameters to excel file
detecEval<-detecEvalFinal[0,]
if(dettype=="spread"|dettype=="combined"){
  detecEval[1,]<-c(spec,"all",paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,numTP,numFP,numFN,TPhitRate,TPR,TPdivFP,paste(allowedZeros,collapse=","),paste(grpsize,collapse=","),paste(detskip,collapse=","),paste(groupInt,collapse=","),timediff,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")
  }else{
detecEval[1,]<-c(spec,"all",paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,numTP,numFP,numFN,TPhitRate,TPR,TPdivFP,NA,NA,NA,NA,timediff,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")   
  }
detecEvalFinal <- rbind(detecEvalFinal,detecEval)
 
detecEval2<-read.csv(paste(outputpath,"DetectorRunLog.csv",sep=""))
detecEvalFinal<-rbind(detecEval2,detecEvalFinal)

beep(10)

write.csv(detecEvalFinal,paste(outputpath,"DetectorRunLog.csv",sep=""),row.names=FALSE)

###################
#MAYBE TEMPORARY- Save dataset for next steps so don't have to rerun after crash

#write.csv(DetecTab2,paste(outputpath,"DetecTab2.csv",sep=""),row.names=FALSE)


###################

###############################################


#RUN RANDOM FOREST MODEL 


################################################
runname<-runname
spec<-spec
#Which data would you like to evaluate?
#species
yminn<-0 #for spec plotting. Should be the same as detector window preset
ymaxx<-1000 #" "

detfiles<-list.files(paste("E:/DetectorRunOutput/",runname,sep=""),pattern = "RF")  

#extract mooring names from moorings used in run
mooringpat=NULL
for(n in 1:length(detfiles)){
  mooringpat<-c(mooringpat,substr(detfiles[n],1,11))
}

if(whiten=="y"){
  soundfile<-(paste("Bbandp",LMS*100,"x_FO",FO,sep=""))
}else{
  soundfile<-"No_whiten"
}

#only choose soundfiles that match those used in run
soundfiles<-NULL
for(n in 1:length(dir(paste("E:/Combined_sound_files/",spec,"/",soundfile,sep="")))){
  if(substr(dir(paste("E:/Combined_sound_files/",spec,"/",soundfile,sep=""))[n],1,11) %in% mooringpat){
    soundfiles<-c(soundfiles,dir(paste("E:/Combined_sound_files/",spec,"/",soundfile,sep=""))[n])
  }
}

#make sure they are in same order
soundfiles<-sort(soundfiles)
detfiles<-sort(detfiles)


data=NULL
for(n in 1:length(detfiles)){
  data2<-read.csv(paste("E:/DetectorRunOutput/",runname,"/",detfiles[n],sep=""), sep = "\t")
  data2<-cbind(soundfiles[n],data2)
  data<-rbind(data,data2)
}

#translate response to binary 
data$detectionType<- as.character(data$detectionType)
data[which(data$detectionType=="TP"),9]<-1
data[which(data$detectionType=="FP"),9]<-0
data[which(data$detectionType=="FN"),9]<-2

#remove FN from data
data<-data[which(data$detectionType==0|data$detectionType==1),]
data$detectionType<-as.numeric(data$detectionType)

#make interference columns into factors
if(length(data)>9){
  for(n in 10:length(data)){
    data[,n]<-as.factor(data[,n])
  }
}
#######1 mooring test######
#data<-data[which(data$`soundfiles[n]`=="BS15_AU_02a_files1-104.wav"),]

data<-splitdf(data,weight = 1/2)[[1]]

data<-spectral_features(data,1)

data2<-data[,9:length(data)]
data2$detectionType<-as.factor(data2$detectionType)
#names(data2)[1]<-"Mooring"

my.xval = list()
my.xval$predictions = list()
my.xval$labels = list()

#number of iterations 
CV=150

#set desired TPR threshold
TPRthresh<-.95

AUC_avg<-c()
f=1
CUTvec=NULL
#how many cross validation runs you want
for(p in 1:CV){
  print(paste("model",p))
  train<-splitdf(data2,weight = 2/3)
  data.rf<-randomForest(formula=detectionType ~ .,data=train[[1]],mtry=7)
  pred<-predict(data.rf,train[[2]],type="prob")
  ROCRpred<-prediction(pred[,2],train[[2]]$detectionType)
  auc.perf = performance(ROCRpred, measure = "auc",plot=F)
  AUC_avg<-c(AUC_avg,as.numeric(auc.perf@y.values))
  
  my.xval$predictions[[p]] = ROCRpred@predictions[[1]]
  my.xval$labels[[p]] = ROCRpred@labels[[1]]
  
#  prob.perf = performance(ROCRpred, "tpr","fpr")
  
#  TPR<-NULL
#  TPR <- data.frame(cut=prob.perf@alpha.values[[1]], tpr=prob.perf@y.values[[1]])
#  CUT <- max(TPR[which(TPR$tpr>=TPRthresh),1])
#  CUTvec<-c(CUTvec,CUT)
  
#  if(f==1){
#    probstab<-pred[,2]
#  }else{
#    probstab<-data.frame(probstab,pred[,2])
#  }  
#  f=f+1
}

#first round of graphs: maybe could merge this section and the next, but mainly copied this code so keep seperate for now 
pp = my.xval$predictions
ll = my.xval$labels
predd = prediction(pp, ll)
perff = performance(predd, "tpr", "fpr")

#no avg
plot(perff, xaxs="i", yaxs="i",main=paste("All",CV," cross validation runs"))
abline(a=0, b= 1)

#avg. Don't really know what the points on the line mean. 
plot(perff, avg = "vertical", spread.estimate = "stddev",spread.scale=2, xaxs="i", yaxs="i", 
     #show.spread.at=c(.05,.075,.1,.125,.15,.2,.3),
     lwd = 2, main = paste("Vertical avg w/ std dev\n"))
plot(perff, avg = "threshold",  xaxs="i", yaxs="i", spread.scale=2,
      lwd = 2, main = paste("Threshold avg"),colorize=T)
abline(a=0, b= 1)
print(mean(AUC_avg))

varImpPlot(data.rf,  
           sort = 27,
           n.var=27,
           main="Top 10 - Variable Importance")

#save last model
save(data.rf, file = paste("E:/DetectorRunOutput/",runname,"/an_example_model.rda",sep=""))

##Graph std error and probability rates, with true detection included 
#probmean<-NULL
#probstderr<-NULL
#for(x in 1:nrow(probstab)){
#  probmean[x]<-mean(as.numeric(probstab[x,]))
#  probstderr[x]<-std.error(as.numeric(probstab[x,]))
#}

#CUTmean<-mean(CUTvec)
#CUTstd.err<-std.error(CUTvec)

#data3<- data2

##assuming $detection type is already in this data NOTE not 
#data3$probmean<-probmean
#data3$probstderr<-probstderr
#data3$probn<-CV

##no avg
#colfunc <- colorRampPalette(c("red", "green"))

#plot(as.numeric(probmean),probstderr, col = ifelse(data3$detectionType==1,'green','red'))
#abline(h=CUTstd.err)
#abline(v=CUTmean)

#this looks like cleanest portion of data- but how to subset to this while keeping a known TPR? Even if it takes a after the fact analysis, should explore only taking the "tail" of the data
#plot(as.numeric(probmean),probstderr, col = ifelse(((as.numeric(probmean) < CUTmean)|(as.numeric(probstderr)>CUTstd.err)),'red','green'))

#cor.test(as.numeric(probmean),probstderr)





#run 1 full data: can achieve .95 TPR with FPR of .7 84% of TPs total accounted for. With FPR of a little of .5, can get TPR of .9, 80% of TPs total accounted for. 
#about the same for second .25 sampling of data. 
#.95 tpr .7 fpr seems more consistently on curve line than .9/.5. 




###############################################


#EMPLOY MODEL ON FULL DATASETS

################################################

allDataPath<-"E:/Datasets"
allMoorings<-dir(allDataPath)[1] #Just AW12_AU_BS3 right now, need fully analyzed GT to test on full mooring

if(whiten!="y"){
fileSizeInt<-345
}else{
fileSizeInt<-(345*3) #whitened files are smaller so still under 6 gigs. 
}

if(moorType=="HG"){
  sfpath<-paste("E:/Datasets/",dir(allDataPath)[1],"/",spec,"_ONLY_yesUnion",sep = "")
}else{
  sfpath<-paste("E:/Full_datasets/",dir(allDataPath)[1],sep = "")
}

#
#run sound files:
resltsTab <- NULL
resltsTabInt<- NULL
durTab<-NULL
for(m in allMoorings){
  sound_files <- dir(sfpath) #
  #make 300 increment break points for sound files. SoX and RRaven can't handle full sound files. 
  bigFile_breaks<-c(seq(1,length(sound_files),fileSizeInt),length(sound_files)) #[sample.int(58,size=2,replace=F)] #last index for run test. 

  #if end of file happens to end on last break make sure its not redundant
  if(bigFile_breaks[length(bigFile_breaks)]==bigFile_breaks[length(bigFile_breaks)-1]){
    bigFile_breaks<-bigFile_breaks[1:(length(bigFile_breaks)-1)]}
  for(b in 1:(length(bigFile_breaks)-1)){
    sound_files <- dir(sfpath)[bigFile_breaks[b]:bigFile_breaks[b+1]]
    sound_filesfullpath <- paste(sfpath,"/",sound_files,sep = "")
    if(whiten=="n" & moorType=="HG"){
      whiten2<-"Entire_No_whiten"
      
      combSound<-paste(startcombpath,spec,"/",whiten2,"/",m,"_files_entire",bigFile_breaks[b],".wav",sep="")
      if(file.exists(combSound)){
        durTab <-read.csv(paste(startcombpath,spec,"/",whiten2,"/SFiles_and_durations.csv",sep=""))  
      }else{
        dir.create(paste(startcombpath,spec,sep=""))
        dir.create(paste(startcombpath,spec,"/",whiten2,"/",sep=""))
        print(paste("Creating file ",m,bigFile_breaks[b],sep=""))
        sox_alt(paste(noquote(paste(paste(sound_filesfullpath,collapse=" ")," ",combSound,sep=""))),exename="sox.exe",path2exe="E:\\Accessory\\sox-14-4-2")
        durTab<-duration_store()
      }
      
      filePath<- paste(startcombpath,spec,whiten2,sep="")
      
    }else if(whiten=="n" & moorType!="HG"){
      whiten2<-"Entire_full_No_whiten"
      
      combSound<-paste(startcombpath,"/",whiten2,"/",m,"_files_entire",bigFile_breaks[b],".wav",sep="")
      if(file.exists(combSound)){
        durTab <-read.csv(paste(startcombpath,whiten2,"/SFiles_and_durations.csv",sep=""))   
      }else{
        dir.create(paste(startcombpath,"/",whiten2,"/",sep=""))
        print(paste("Creating file ",m,bigFile_breaks[b],sep=""))
        sox_alt(paste(noquote(paste(paste(sound_filesfullpath,collapse=" ")," ",combSound,sep=""))),exename="sox.exe",path2exe="E:\\Accessory\\sox-14-4-2")
        durTab<-duration_store()
      }
      
      filePath<- paste(startcombpath,whiten2,sep="")
    }

    if(whiten=="y" & moorType=="HG"){
      whiten2 <- paste("Entire_Bbandp",100*LMS,"x_","FO",FO,sep = "")
      filePath<- paste(startcombpath,spec,whiten2,sep="")
      durTab <-read.csv(paste(startcombpath,spec,"/",whiten2,"/SFiles_and_durations.csv",sep=""))  
    }else if(whiten=="y" & moorType!="HG"){
      whiten2 <- paste("Entire_full_Bbandp",100*LMS,"x_","FO",FO,sep = "")
      filePath<- paste(startcombpath,whiten2,sep="")
      durTab <-read.csv(paste(startcombpath,whiten2,"/SFiles_and_durations.csv",sep=""))  
    }
  }
  
  #run pulse and fin/mooring detector, if selected:
  if(interfere=="y"){
    for(b in bigFile_breaks[1:length(bigFile_breaks)-1]){
      combname<- paste(m,"_files_entire",b,".wav",sep="")
      for(i in interfereVec){
        print(paste("Running detector for",combname))
        resltVarInt <- raven_batch_detec(raven.path = ravenpath, sound.files = combname, path =filePath,detector = "Band Limited Energy Detector",dpreset=i,vpreset="RW_Upcalls")
        resltVarInt$Mooring<-paste(m,"_files_entire",b,".wav",sep="")
        resltVarInt$detector<-i
        resltVarInt$detectorType<-"intereference"
        resltVarInt$detectorCount<-which(interfereVec==i)
        if(is.null(nrow(resltVar))==resltVarInt){
        resltsTabInt<-rbind(resltsTabInt,resltVarInt)
        }
        resltVarInt<-NULL
        }
      }
    }
  #run detector(s)
  if(dettype=="spread"|dettype=="combined"){
    for(b in bigFile_breaks[1:length(bigFile_breaks)-1]){
      combname<- paste(m,"_files_entire",b,".wav",sep="")
      for(q in 1:length(detectorssprshort)){
        for(r in detectorssprshort[[q]]){
          print(paste("Running detector for",combname))
          resltVar <- raven_batch_detec(raven.path = ravenpath, sound.files = combname, path = filePath,detector = "Band Limited Energy Detector",dpreset=r,vpreset="RW_Upcalls")
          resltVar$Mooring<-paste(m,"_files_entire",b,".wav",sep="")
          resltVar$detector<-r
          resltVar$detectorType<-"spread"
          resltVar$detectorCount<-q
          if(is.null(nrow(resltVar))==FALSE){
          resltsTab<- rbind(resltsTab,resltVar)
          }
          resltVar<-NULL 
        }
      }  
    }
  }
  
  if(dettype=="single"|dettype=="combined"){
    for(b in bigFile_breaks[1:length(bigFile_breaks)-1]){
      combname<- paste(m,"_files_entire",b,".wav",sep="")
      for(n in detectorssinshort){
        print(paste("Running detector for",combname))
        resltVar <- raven_batch_detec(raven.path = ravenpath, sound.files =  combname, path = filePath,detector = "Band Limited Energy Detector",dpreset=n,vpreset ="RW_Upcalls")
        resltVar$Mooring<-paste(m,"_files_entire",b,".wav",sep="")
        resltVar$detector<-n
        resltVar$detectorType<-"single"
        resltVar$detectorCount<-which(detectorssinshort==n)
        if(is.null(nrow(resltVar))==FALSE){
        resltsTab<- rbind(resltsTab,resltVar)
        }
        resltVar<-NULL
      }
    }  
  }
}

#write durTab to file. 1st time run will set but will not modify durTab after in any case so no need for conditional
write.csv(durTab,paste(filePath,"/SFiles_and_durations.csv",sep=""),row.names = F)

DetecTab2<-process_data(2)

#Define table for later excel file export. 
colClasses = c("character","character","character","character","character","numeric","numeric","numeric", "numeric","numeric","numeric","character","character","character","character","character","character","character","character","numeric","numeric","character")
detecEvalFinal <- read.csv(text="Species, Moorings, Detectors, DetType, RunName, numTP, numFP, numFN, TPhitRate, TPR, TPdivFP, ZerosAllowed,GroupSize,SkipAllowance,GroupInterval,TimeDiff,TimeDiffself,MinMaxDur,numDetectors,FO,LMS,Notes", colClasses = colClasses)


MoorTab<-NULL
MoorVar<-NULL
for(v in 1:length(unique(DetecTab2$Mooring))){
  print(paste("Adding interference detector variables to",sort(unique(DetecTab2$Mooring))[v]))   
  MoorVar<-DetecTab2[which(DetecTab2$Mooring==sort(unique(DetecTab2$Mooring))[v]),]
  
  #Define useful comlumns in MoorVar
  sound.files <- MoorVar[,12]
  MoorVar <- MoorVar[,c(1:7,13,14:17)]
  MoorVar<-cbind(sound.files,MoorVar)
  
    #compare detections to sources of interference
    if(interfere=="y"){
      for(n in 1:max(resltsTabInt$detectorCount)){
        MoorInt<-resltsTabInt[which(resltsTabInt$Mooring==sort(unique(DetecTab2$Mooring))[v]),]
        MoorVar[,n+14]<-0
        colnames(MoorVar)[length(MoorVar)]<-paste(resltsTabInt[which(resltsTabInt$detectorCount==n),10])[1]
        MoorInt<-MoorInt[which(MoorInt$detectorCount==n),]
        print(paste("       Compare with detector",MoorInt[1,10]))   
        for(g in 1:nrow(MoorVar)){
          hvec <- which(MoorInt$meantime<(MoorVar$meantime[g]+2)&MoorInt$meantime>(MoorVar$meantime[g]-2))
          if(length(hvec>0)){
          for(h in min(hvec):max(hvec)){
            if((MoorInt[h,4]<MoorVar[g,9] & MoorInt[h,5]> MoorVar[g,9])|(MoorInt[h,4]> MoorVar[g,5] & MoorInt[h,5]< MoorVar[g,6])){
              MoorVar[g,n+14]<-1
            }
          }
        }
      }
      }
    }
  
  #write.table(MoorVar[,2:7],paste(outputpath,runname,"/",unique(DetecTab2$Mooring)[v],"FINAL_Summary_",dettype,"_Ravenformat",".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
  
  MoorTab<-rbind(MoorTab,MoorVar)
  #MoorVar$Moorpred<-c(MoorVar$Moorpred,predict(data.rf,MoorVar,type="prob"))
}

findata<-MoorTab

#create columns to be used later 
findata$probmean<-0
findata$probstderr<-0
findata$probn<-0
findata$detectionType<-0

#make interference columns into factors
if(length(findata)>17){
  for(n in 18:length(findata)){
    findata[,n]<-as.factor(findata[,n])
  }
}

findata<-spectral_features(findata,2)

#Generate and run a set amount of models from the original GT data. Probabilities are averaged for each mooring. 
data2$detectionType<-as.factor(data2$detectionType)

AUC_avg<-c()
#generate models and apply to whole dataset. Should keep model parameters the same as when you assessed accuracy to have an idea of reliability. 
f=1
CUTvec=NULL
for(p in 1:CV){
  print(paste("model",p))
  train<-splitdf(data2,weight = 2/3)
  #apply model to data
  data.rf<-randomForest(formula=detectionType ~ .,data=train[[1]],mtry=7)
  pred<-predict(data.rf,findata[,17:length(findata)],type="prob")
  #assess same model performance on GT data 
  pred2<-predict(data.rf,train[[2]],type="prob")
  ROCRpred<-prediction(pred2[,2],train[[2]]$detectionType)
  perf = performance(ROCRpred, "tpr","fpr")
  
  TPR<-NULL
  TPR <- data.frame(cut=perf@alpha.values[[1]], tpr=perf@y.values[[1]])
  CUT <- max(TPR[which(TPR$tpr>=TPRthresh),1])
  CUTvec<-c(CUTvec,CUT)
  

  if(f==1){
    probstab<-pred[,2]
  }else{
    probstab<-data.frame(probstab,pred[,2])
  }  
  f=f+1
}

probmean<-NULL
probstderr<-NULL
for(x in 1:nrow(probstab)){
probmean[x]<-mean(as.numeric(probstab[x,]))
probstderr[x]<-std.error(as.numeric(probstab[x,]))
}

CUTmean<-mean(CUTvec)
CUTstd.err<-std.error(CUTvec)

findata$probmean<-probmean
findata$probstderr<-probstderr
findata$probn<-CV
findata[,1]<-substr(findata$sound.files,1,11)
TPtottab<-data.frame(TPtot,GTtot,MoorCor)

#apply models and average probabilities. ASSUMPTIONS for stats to make sense: Full mooring is run, and all HG data is included in GT. 
MoorVar<-NULL
for(v in 1:length(unique(findata$sound.files))){
  
  MoorVar1<-findata[which(findata$sound.files==sort(unique(findata$sound.files))[v]),]
  
  #MoorVar1$detectionType<-ifelse(((MoorVar1$probstderr<CUTstd.err)&(MoorVar1$probstmean>CUTmean)),"RFselected","RFrejected")
  MoorVar1$detectionType<-ifelse(MoorVar1$probmean>CUTmean,"RFselected","RFrejected")
  MoorVar1<-MoorVar1[which(MoorVar1$detectionType=="RFselected"),]
  
  numTP<-TPtottab[which(TPtottab$MoorCor==sort(unique(findata$sound.files))[v]),1]*TPRthresh
  numTPtruth<-TPtottab[which(TPtottab$MoorCor==sort(unique(findata$sound.files))[v]),2]
  detTotal<-nrow(MoorVar1)
  numFP<-detTotal-numTP
  numFN<-numTPtruth-numTP
  
  TPhitRate <- numTP/numTPtruth*100
  TPR <- numTP/(numTP+numFN)
  TPdivFP<- numTP/numFP
  
  #save stats and parameters to excel file
  detecEval<-detecEvalFinal[0,]
  if(dettype=="spread"|dettype=="combined"){
    detecEval[1,]<-c(spec,paste("full",sort(unique(findata$sound.files))[v]),paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,numTP,numFP,numFN,TPhitRate,TPR,TPdivFP,paste(allowedZeros,collapse=","),paste(grpsize,collapse=","),paste(detskip,collapse=","),paste(groupInt,collapse=","),timediff,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")
  }else{
    detecEval[1,]<-c(spec,paste("full",sort(unique(findata$sound.files))[v]),paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,numTP,numFP,numFN,TPhitRate,TPR,TPdivFP,NA,NA,NA,NA,timediff,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")   
  }
  
  detecEval2<-read.csv(paste(outputpath,"DetectorRunLog.csv",sep=""))
  detecEvalFinal<-rbind(detecEval2,detecEval)
  write.csv(detecEvalFinal,paste(outputpath,"DetectorRunLog.csv",sep=""),row.names=FALSE)
  
  MoorVar1<-MoorVar1[,c(2:8)]
  
  MoorVar2<-findata[which(findata$sound.files==sort(unique(findata$sound.files))[v]),]
  
  MoorVar3<-data.frame(MoorVar2$Selection,MoorVar2$FileOffsetBegin,MoorVar2$FileOffsetEnd,MoorVar2$`Low Freq (Hz)`,MoorVar2$`High Freq (Hz)`,MoorVar2$sound.files,MoorVar2$File,MoorVar2$probmean,MoorVar2$probstderr,MoorVar2$probn)
  colnames(MoorVar3)<-c("Selection","FileOffsetBegin","FileOffsetEnd","Low Freq (Hz)","High Freq (Hz)","Mooring","File","probs","probsstderr","probsn")
  
  write.table(MoorVar1,paste(outputpath,runname,"/",sub(".wav", "", sort(unique(findata$sound.files))[v]),"FINAL_Model_Applied_Ravenformat",".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
  write.table(MoorVar2,paste(outputpath,runname,"/",sub(" .wav", "", sort(unique(findata$sound.files))[v]),"FINAL_Model_Applied_probs",".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
  
}

#After we save and report the stats, could make a section here for analyzing probabilites and variance. 
  
#no avg
colfunc <- colorRampPalette(c("red", "green"))

plot(as.numeric(probmean),probstderr, col = ifelse(as.numeric(probstderr) > CUTstd.err,'red','green'))
plot(as.numeric(probmean),probstderr, col = ifelse(as.numeric(probmean) < CUTmean,'red','green'))

#this looks like cleanest portion of data- but how to subset to this while keeping a known TPR? Even if it takes a after the fact analysis, should explore only taking the "tail" of the data
plot(as.numeric(probmean),probstderr, col = ifelse(((as.numeric(probmean) < CUTmean)|(as.numeric(probstderr)>CUTstd.err)),'red','green'))

cor.test(as.numeric(probmean),probstderr)


  