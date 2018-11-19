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
library(stringi)
library(signal)

freqstat.normalize<- function(freqstat,lowFreq,highFreq){
  newstat<-((freqstat-lowFreq)/(highFreq-lowFreq))
  return(newstat)
}

frd_wrblr_int <- function(spc, fsmooth = 0.1, threshold = 10, wn = "hanning", bp = NULL,sr=NULL)
{
  
  # get frequency windows length for smoothing
  step <- sr/512/1000
  
  fsmooth <- fsmooth/step
  
  # number of samples
  n <- nrow(spc)
  
  # smoothing parameter
  FWL <- fsmooth - 1
  
  # smooth 
  z <- apply(as.matrix(1:(n - FWL)), 1, function(y) sum(spc[y:(y + FWL), 2]))
  zf <- seq(min(spc[,1]), max(spc[,1]), length.out = length(z))
  
  # make minimum amplitude 0
  z <- z - min(z)
  z[z < 0] <- 0
  
  # normalize amplitude from 0 to 1
  z <- z/max(z)
  
  meanpeakf <- zf[which.max(z)] + (step / 2)
  
  # return low and high freq
  return(meanpeakf) 
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

context_sim <-function(sdata){
  #context simulator- add or subtract % points based on how good neighboring calls were. Only useful for full mooring dataset. 
  for(w in 1:length(unique(sdata[,1]))){
    datVar<-sdata[which(sdata[,1]==unique(sdata[,1])[w]),]
    datVar$probrollscore<-0
    datVar$probmean2<-datVar$probmean
    for(n in 1:(nrow(datVar)-1)){
      if(datVar$probmean2[n]>=greatcallThresh){
        datVar$probrollscore[n+1]<-maxBonus
        if(datVar$probrollscore[n+1]>maxBonus){
          datVar$probrollscore[n+1]<-maxBonus
        }
      }else if(datVar$probmean2[n]>=(-maxPenalty+CUTmean)&datVar$probmean2[n]<greatcallThresh){
        datVar$probrollscore[n+1]<-datVar$probrollscore[n]+(datVar$probmean2[n]*goodcallBonus)
        if(datVar$probrollscore[n+1]>maxBonus){
          datVar$probrollscore[n+1]<-maxBonus
        }
      }else{
        datVar$probrollscore[n+1]<-datVar$probrollscore[n]+badcallPenalty
        if(datVar$probrollscore[n+1]<maxPenalty){
          datVar$probrollscore[n+1]<-maxPenalty
        }
      }
      if(datVar$probmean2[n]+datVar$probrollscore[n]>0 & datVar$probmean2[n]+datVar$probrollscore[n]<1){
        datVar$probmean2[n]<-datVar$probmean2[n]+datVar$probrollscore[n]
      }else if(datVar$probmean2[n]+datVar$probrollscore[n]<0){
        datVar$probmean2[n]<-0
      }else if(datVar$probmean2[n]+datVar$probrollscore[n]>1){
        datVar$probmean2[n]<-1
      }
    }
  #same but backwards through data 
  datVar$probrollscore<-0
  datVar$probmean3<-datVar$probmean
  for(n in (nrow(datVar)-1):1){
    if(datVar$probmean3[n]>=greatcallThresh){
      datVar$probrollscore[n+1]<-maxBonus
      if(datVar$probrollscore[n+1]>maxBonus){
        datVar$probrollscore[n+1]<-maxBonus
      }
    }else if(datVar$probmean3[n]>(-maxPenalty)&datVar$probmean3[n]<greatcallThresh){
      datVar$probrollscore[n+1]<-datVar$probrollscore[n]+(datVar$probmean3[n]*goodcallBonus)
      if(datVar$probrollscore[n+1]>maxBonus){
        datVar$probrollscore[n+1]<-maxBonus
      }
    }else{
      datVar$probrollscore[n+1]<-datVar$probrollscore[n]+badcallPenalty
      if(datVar$probrollscore[n+1]<maxPenalty){
        datVar$probrollscore[n+1]<-maxPenalty
      }
    }
    if(datVar$probmean3[n]+datVar$probrollscore[n]>0 & datVar$probmean3[n]+datVar$probrollscore[n]<1){
      datVar$probmean3[n]<-datVar$probmean3[n]+datVar$probrollscore[n]
    }else if(datVar$probmean3[n]+datVar$probrollscore[n]<0){
      datVar$probmean3[n]<-0
    }else if(datVar$probmean3[n]+datVar$probrollscore[n]>1){
      datVar$probmean3[n]<-1
    }
  }
  
  datVar$probmean<-(datVar$probmean2+datVar$probmean3)/2 #average
  #datVar$probmean<-pmax(datVar$probmean2,datVar$probmean3) #max
  #datVar$probmean<-pmin(datVar$probmean2,datVar$probmean3) #min
  
  sdata<-sdata[which(sdata[,1]!=unique(sdata[,1])[w]),]
  datVar<-datVar[,c(1:length(sdata))]
  sdata<-rbind(sdata,datVar)
  }
  return(sdata)
}

after_model_write <-function(mdata,finaldatrun){
  mdata[,1]<-substr(mdata[,1],1,11)
  MoorVar1<-NULL
  for(v in 1:length(unique(mdata[,1]))){
    if(finaldatrun==1){
      name<-"RFA"
    }else{
      name<-"full"
    }
    MoorVar1<-mdata[which(mdata[,1]==sort(unique(mdata[,1]))[v]),]
    
    MoorVar1<-MoorVar1[which(MoorVar1$probmean>CUTmean),]
    
    numTPtruth<-TPtottab[which(stri_detect_fixed(sort(unique(mdata[,1]))[v],TPtottab$MoorCor)),2]
    
    if(finaldatrun==2){
    numTP<-TPtottab[which(stri_detect_fixed(sort(unique(mdata[,1]))[v],TPtottab$MoorCor)),1]*TPRthresh
    detTotal<-nrow(MoorVar1)
    numFP<-detTotal-numTP
    numFN<-numTPtruth-numTP
    }else{
      #use real TP values for GT data. 
      numTP<-sum(as.numeric(as.character(MoorVar1$detectionType)))
      numFP<-(nrow(MoorVar1)-numTP)
      numFN<-numTPtruth-numTP
    }
    TPhitRate <- numTP/numTPtruth*100
    TPR <- numTP/(numTP+numFN)
    TPdivFP<- numTP/numFP
    
    #save stats and parameters to excel file
    detecEval<-detecEvalFinal[0,]
    detecEval[,2]<-as.character(detecEval[,2])
    detecEval[,13]<-as.character(detecEval[,13])
    if(dettype=="spread"|dettype=="combined"){
      detecEval[1,]<-c(spec,paste(name,sort(unique(mdata[,1]))[v]),paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,numTP,numFP,numFN,TPhitRate,TPR,TPdivFP,AUCadj,paste(CV,TPRthresh,sep=","),paste(greatcallThresh,-maxPenalty,sep=","),paste(maxBonus,goodcallBonus,badcallPenalty,sep=","), paste(allowedZeros,collapse=","),paste(grpsize,collapse=","),paste(downsweepCompMod,downsweepCompAdjust,sep=","),paste(detskip,collapse=","),paste(groupInt,collapse=","),timediff,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")
    }else{
      detecEval[1,]<-c(spec,paste(name,sort(unique(mdata[,1]))[v]),paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,numTP,numFP,numFN,TPhitRate,TPR,TPdivFP,AUCadj,paste(CV,TPRthresh,sep=","),paste(greatcallThresh,-maxPenalty,sep=","),paste(maxBonus,goodcallBonus,badcallPenalty,sep=","), paste(allowedZeros,collapse=","),paste(grpsize,collapse=","),paste(downsweepCompMod,downsweepCompAdjust,sep=","),paste(detskip,collapse=","),paste(groupInt,collapse=","),timediff,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")
    }
    
    detecEval2<-read.csv(paste(outputpath,"DetectorRunLog.csv",sep=""))
    detecEvalFinal<-rbind(detecEval2,detecEval)
    write.csv(detecEvalFinal,paste(outputpath,"DetectorRunLog.csv",sep=""),row.names=FALSE)

    if(finaldatrun==2){    
    MoorVar1<-MoorVar1[,c(2:8)]
    
    write.table(MoorVar1,paste(outputpath,runname,"/",sub(".wav", "", sort(unique(mdata[,1]))[v]),"FINAL_Model_Applied_Ravenformat",".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
    
    MoorVar2<-mdata[which(mdata[,1]==sort(unique(mdata[,1])))[v],]
    MoorVar3<-data.frame(MoorVar2$Selection,MoorVar2$FileOffsetBegin,MoorVar2$FileOffsetEnd,MoorVar2$`Low Freq (Hz)`,MoorVar2$`High Freq (Hz)`,MoorVar2$sound.files,MoorVar2$File,MoorVar2$probmean,MoorVar2$probstderr,MoorVar2$probn)
    colnames(MoorVar3)<-c("Selection","FileOffsetBegin","FileOffsetEnd","Low Freq (Hz)","High Freq (Hz)","Mooring","File","probs","probsstderr","probsn")
    write.table(MoorVar3,paste(outputpath,runname,"/",sub(" .wav", "", sort(unique(mdata[,1]))[v]),"FINAL_Model_Applied_probs",".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
    
    }else{
      RavenExport<-data.frame(MoorVar1$Selection)
      RavenExport[,2]<-"Spectrogram 1"
      RavenExport[,3]<-1
      RavenExport[,4]<-MoorVar1$Begin.Time..s.
      RavenExport[,5]<-MoorVar1$End.Time..s.
      RavenExport[,6]<-MoorVar1$Low.Freq..Hz.
      RavenExport[,7]<-MoorVar1$High.Freq..Hz.
      RavenExport[,8]<-as.numeric(as.character(MoorVar1$detectionType))
      RavenExport[,9]<-as.character(MoorVar1$probmean)
      RavenExport[,10]<-as.character(MoorVar1$probstderr)
      RavenExport[,11]<-as.character(MoorVar1$n)
      
      colnames(RavenExport)<-c("Selection","View","Channel","Begin Time (s)","End Time (s)","Low Freq (Hz)","High Freq (Hz)","detectionType","probs","err","n")
      
      RavenExport[which(RavenExport[,8]==1),8]<-"TP"
      RavenExport[which(RavenExport[,8]==0),8]<-"FP"
      
      write.table(RavenExport,paste(outputpath,runname,"/",sub(" .wav", "", sort(unique(mdata[,1]))[v]),"Model_Applied_probs",".txt",sep=""),quote=FALSE,sep = "\t",row.names=FALSE)
        
      }
  }
  numTPtruth<-sum(TPtottab$GTtot)
  mdata<-mdata[which(mdata$probmean>CUTmean),]
  
  if(finaldatrun==1){
  numTP<-sum(as.numeric(as.character(mdata$detectionType)))
  numFP<-(nrow(mdata)-numTP)
  numFN<-numTPtruth-numTP
  }else{
    
  }
  TPhitRate <- numTP/numTPtruth*100
  TPR <- numTP/(numTP+numFN)
  TPdivFP<- numTP/numFP
  
  detecEval<-detecEvalFinal[0,]
  detecEval[,2]<-as.character(detecEval[,2])
  detecEval[,13]<-as.character(detecEval[,13])
  if(dettype=="spread"|dettype=="combined"){
    detecEval[1,]<-c(spec,paste(name,"all"),paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,numTP,numFP,numFN,TPhitRate,TPR,TPdivFP,AUCadj,paste(CV,TPRthresh,sep=","),paste(greatcallThresh,-maxPenalty,sep=","),paste(maxBonus,goodcallBonus,badcallPenalty,sep=","), paste(allowedZeros,collapse=","),paste(grpsize,collapse=","),paste(downsweepCompMod,downsweepCompAdjust,sep=","),paste(detskip,collapse=","),paste(groupInt,collapse=","),timediff,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")
  }else{
    detecEval[1,]<-c(spec,paste(name,"all"),paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,numTP,numFP,numFN,TPhitRate,TPR,TPdivFP,AUCadj,paste(CV,TPRthresh,sep=","),paste(greatcallThresh,-maxPenalty,sep=","),paste(maxBonus,goodcallBonus,badcallPenalty,sep=","), paste(allowedZeros,collapse=","),paste(grpsize,collapse=","),paste(downsweepCompMod,downsweepCompAdjust,sep=","),paste(detskip,collapse=","),paste(groupInt,collapse=","),timediff,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")
  }
  
  
  detecEval2<-read.csv(paste(outputpath,"DetectorRunLog.csv",sep=""))
  detecEvalFinal<-rbind(detecEval2,detecEval)
  write.csv(detecEvalFinal,paste(outputpath,"DetectorRunLog.csv",sep=""),row.names=FALSE)
}

adaptive_compare<-function(Compdata,specfeatrun){
  for(a in 1:2){#go through twice in case there are mulitple boxes close to one another. 
  for(o in unique(Compdata[,1])){
    CompVar<-Compdata[which(Compdata[,1]==o),]
    CompVar<-CompVar[order(CompVar$meantime),]
    print(paste("For mooring",o))
    r=0
    n=0
    newrow<-CompVar[0,]
    IDvec<-NULL
    for(q in 1:(nrow(CompVar)-1)){
      if(CompVar$meantime[q+1]<=(CompVar$meantime[q]+timediffself)){
        if(CompVar$probmean[q+1]+probdist<CompVar$probmean[q]|CompVar$probmean[q+1]-probdist>CompVar$probmean[q]){#take only the best one
         # print(q)
          newdat<-CompVar[0,]
          newdat[1:2,]<-CompVar[c(q,q+1),]
          newdat<-newdat[order(newdat$probmean),]
          IDvec<-c(IDvec,newdat$Selection)
          newrow[r+1,]<-newdat[2,]
          r=r+1
        }else{
          #print(q)
          newdat<-CompVar[0,]
          newdat[1:2,]<-CompVar[c(q,q+1),]
          IDvec<-c(IDvec,newdat$Selection)
          s<-as.numeric(min(newdat[,3]))
          e<-as.numeric(max(newdat[,4]))
          l<-as.numeric(min(newdat[,5]))
          h<-as.numeric(max(newdat[,6]))
          dt<-max(as.numeric(as.character(newdat[,7])))
          mt<-(s+e)/2
          mf<-(h+l)/2
          fr<-(h-l)
          
          newrow[r+1,]<-data.frame(newdat[1,1],newdat[1,2],s,e,l,h,dt,mf,fr,mt)
          newrow[r+1,]<-spectral_features(newrow[r+1,],specfeatrun)
          
          newrow[r+1,]$probmean<-mean(newdat$probmean)
          newrow[r+1,]$probstderr<-mean(newdat$probstderr)
          newrow[r+1,]$n<-mean(newdat$n)
                                   
          r=r+1
        }
          
        }

      }
      if(r>0){
        CompVar<-subset(CompVar,!(CompVar$Selection %in% IDvec))
        CompVar<-rbind(CompVar,newrow)
        CompVar<-CompVar[order(CompVar$meantime),]
        n=n+1
      }
    if(n>0){
      Compdata<-Compdata[-which(Compdata[,1]==o),]
      Compdata<-rbind(Compdata,CompVar)      
    }

  }
  }
 return(Compdata) 
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
  print(z)
  #store reused calculations to avoid indexing 
  Low<-specdata$Low.Freq..Hz.[z]
  High<-specdata$High.Freq..Hz.[z]
  Start<-specdata$Begin.Time..s.[z]
  End<-  specdata$End.Time..s.[z]
  
  foo <- foo2 <- readWave(paste(specpath,specdata[z,1],sep=""),Start,End,units="seconds")
  sample_rate<-foo@samp.rate
  foo<-ffilter(foo,from=Low,to=High,output="Wave")
  foo.spec <- spec(foo,plot=F, PSD=T)
  #foo.spec <- foo.spec[which(foo.spec[,1]<(High/1000)&foo.spec[,1]>(Low/1000)),]#,ylim=c(specdata$Low.Freq..Hz.[z],specdata$High.Freq..Hz.[z])
  foo.specprop <- specprop(foo.spec) #
  #spectro(foo) #could do image analysis on this guy 
  foo.meanspec = meanspec(foo, plot=F,ovlp=90)#not sure what ovlp parameter does but initially set to 90 #
  foo.meanspec2 = meanspec(foo2, plot=F,ovlp=90)
  #foo.meanspec2<-foo.meanspec[which(foo.meanspec[,1]<High&foo.meanspec[,1]>Low),] #only returns like 2 numbers- probably not getting a lot out of this if we freq limit
  #foo.meanspec.db = meanspec(foo, plot=F,ovlp=90,dB="max0",flim=c(specdata$Low.Freq..Hz.[z]/1000,specdata$High.Freq..Hz.[z]/1000))#not sure what ovlp parameter does but initially set to 90 #,flim=c(specdata$Low.Freq..Hz.[z]/1000,specdata$High.Freq..Hz.[z]/1000)
  foo.autoc = autoc(foo, plot=F) #
  foo.dfreq = dfreq(foo, plot=F, ovlp=90) #tried bandpass argument, limited dfreq to only 2 different values for some reason. Seemed wrong. 
  Startdom<-foo.dfreq[,2][1]
  Enddom<-foo.dfreq[,2][length(foo.dfreq[,2])]
  Mindom <- min(foo.dfreq, na.rm = TRUE)
  Maxdom <- max(foo.dfreq, na.rm = TRUE)
  Dfrange <- Maxdom - Mindom
  specdata$rugosity[z] = rugo(foo@left / max(foo@left)) 
  specdata$crest.factor[z] = crest(foo)$C
  foo.env = seewave:::env(foo, plot=F) 
  specdata$temporal.entropy[z] = th(foo.env)
  specdata$shannon.entropy[z] = sh(foo.spec)
  specdata$spectrum.roughness[z] = roughness(foo.meanspec[,2])
  specdata$autoc.mean[z] = freqstat.normalize(mean(foo.autoc[,2], na.rm=T),Low,High)
  specdata$autoc.median[z] = freqstat.normalize(median(foo.autoc[,2], na.rm=T),Low,High)
  specdata$autoc.se[z] = std.error(foo.autoc[,2], na.rm=T)
  specdata$dfreq.mean[z] = freqstat.normalize(mean(foo.dfreq[,2], na.rm=T),Low,High)
  specdata$dfreq.se[z] = std.error(foo.dfreq[,2], na.rm=T)
  specdata$specprop.mean[z] = freqstat.normalize(foo.specprop$mean[1],Low,High)
  specdata$specprop.sd[z] = foo.specprop$sd[1]
  specdata$specprop.sem[z] = foo.specprop$sem[1]
  specdata$specprop.median[z] = freqstat.normalize(foo.specprop$median[1],Low,High)
  specdata$specprop.mode[z] = freqstat.normalize(foo.specprop$mode[1],Low,High)
  specdata$specprop.Q25[z] = foo.specprop$Q25[1]
  specdata$specprop.Q75[z] = foo.specprop$Q75[1]
  specdata$specprop.IQR[z] = foo.specprop$IQR[1]
  specdata$specprop.cent[z] = foo.specprop$cent[1]
  specdata$specprop.skewness[z] = foo.specprop$skewness[1]
  specdata$specprop.kurtosis[z] = foo.specprop$kurtosis[1]
  specdata$specprop.sfm[z] = foo.specprop$sfm[1]
  specdata$specprop.sh[z] = foo.specprop$sh[1]
  specdata$specprop.prec[z] = foo.specprop$prec[1]
  specdata$amp.env.median[z] = M(foo)
  specdata$total.entropy[z] = H(foo)
 # specdata$resonant.qual.fact[z]<-Q(foo.meanspec.db,plot=T)$Q #0s introduced
  #warbler params
  specdata$time.ent[z]<-th(foo.env)[1]
  specdata$modindx[z]<- (sum(sapply(2:length(foo.dfreq[,2]), function(j) abs(foo.dfreq[,2][j] - foo.dfreq[,2][j - 1])))/(Dfrange))
  specdata$startdom[z]<-freqstat.normalize(Startdom,Low,High)
  specdata$enddom[z]<-freqstat.normalize(Enddom,Low,High)
  specdata$mindom[z]<-freqstat.normalize(Mindom,Low,High)
  specdata$maxdom[z]<-freqstat.normalize(Maxdom,Low,High)
  specdata$dfrange[z]<-Dfrange
  specdata$dfslope[z]<-((Enddom-Startdom)/(End-Start))
  specdata$meanpeakf[z]<-frd_wrblr_int(foo.meanspec2,sr=sample_rate)

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
        resltsTSPV[resltsTSPV$bottom.freq==sort(unique(resltsTSPV$bottom.freq))[f],13]<-f
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
    
      #remove groups based on grpsize value
      removegrp <- table(resltsTSPV$group)
      resltsTSPV <- subset(resltsTSPV, group %in% names(removegrp[removegrp > (grpsize[d]-1)]))
      #maxgrp<-max(resltsTSPV[,14])
      
      resltsTSPV[1,]<-seq(1:nrow(resltsTSPV))
      #updated algorithm, optimized for performance. avoids r bind
      print(paste("calculating best runs for each group"))
      unwantedSelections<-c()
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
            rsltvec0s<-rle(groupdat[,15+g][! groupdat[,15+g] %in% 98])
            if(any(rsltvec0s$lengths[rsltvec0s$values==0]>allowedZeros[d])){
              break
            }  
            if(RM<grpvec[h+1]&RM+(detskip[d]+1)>grpvec[h+1]&groupdat[h,15]!=groupdat[h+1,15]){
              groupdat[h+1,15+g]<-1
              skipvec<-c(skipvec,(grpvec[h+1]-RM))
              RM<-grpvec[h+1]
            }else if(groupdat[h,15]!=groupdat[h+1,15]&RM+(detskip[d]+1)>grpvec[h+1]&RM-(detskip[d]+1)<=grpvec[h+1]){
              groupdat[h+1,15+g]<-0
            }else if(groupdat[h,15]!=groupdat[h+1,15]&(RM+(detskip[d]+1)<=grpvec[h+1]|RM-(detskip[d]+1)>grpvec[h+1])){
              groupdat[h+1,15+g]<-98
            }
            if(groupdat[h,15]==groupdat[h+1,15]&groupdat[h,15+g]==0&RM<grpvec[h+1]&RM+(detskip[d]+1)>grpvec[h+1]){
              groupdat[h+1,15+g]<-1
              groupdat[h,15+g]<-98
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
        runsum<-runsum[1,] #choose first one
        
        #if run is less than 33% of boxes, build a downsweep. If the downsweep has equal or more ones disqualify it.
        kill="n"
        if((runsum[,2]+1)<grpsize[d]){
          kill="y"
        }
        if(((runsum[,2]+1)*downsweepCompMod)<nrow(groupdat)&kill=="n"){
          groupdat2<- subset(resltsTSPV,group==f)
          groupdat2<-groupdat2[order(groupdat2$meantime,-groupdat2$bottom.freq),]#reverse the order it counts stacks detections
          grpvec2<-groupdat2[,13]
          colClasses = c("numeric","numeric","numeric","numeric","numeric")
          runsum2<- read.csv(text="start, ones, zeros, length,skip", colClasses = colClasses)
          for(g in 1:(nrow(groupdat2)-(grpsize[d]-1))){
            RM2<-groupdat2[g,13]
            groupdat2[,15+g]<-99
            groupdat2[g,15+g]<-2
            skipvec2<-0
            for(h in g:(nrow(groupdat2)-1)){
              rsltvec0s2<-rle(groupdat2[,15+g][! groupdat2[,15+g] %in% 98])
              if(any(rsltvec0s2$lengths[rsltvec0s2$values==0]>allowedZeros[d])){
                break
              }  
              if(RM2>grpvec2[h+1]&RM2-(detskip[d]+1)<grpvec2[h+1]&groupdat2[h,15]!=groupdat2[h+1,15]){
                groupdat2[h+1,15+g]<-1
                skipvec2<-c(skipvec2,(grpvec2[h+1]-RM2))
                RM2<-grpvec2[h+1]
              }else if(groupdat2[h,15]!=groupdat2[h+1,15]&RM2-(detskip[d]+1)<grpvec2[h+1]&RM2+(detskip[d]+1)>=grpvec2[h+1]){
                groupdat2[h+1,15+g]<-0
              }else if(groupdat2[h,15]!=groupdat2[h+1,15]&(RM2-(detskip[d]+1)>=grpvec2[h+1])|RM2+(detskip[d]+1)<grpvec2[h+1]){
                groupdat2[h+1,15+g]<-98
              }
              if(groupdat2[h,15]==groupdat2[h+1,15]&(groupdat2[h,15+g]==0|groupdat2[h,15+g]==98)&RM2>grpvec2[h+1]&(RM2-(detskip[d]+1))<grpvec2[h+1]){
                groupdat2[h+1,15+g]<-1
                groupdat2[h,15+g]<-98
                skipvec2<-c(skipvec2,(grpvec2[h+1]-RM2))
                RM2<-grpvec2[h+1]
              }else if(groupdat2[h,15]==groupdat2[h+1,15]){
                groupdat2[h+1,15+g]<-98
              }
            }
            runsum2[g,1]<-g
            runsum2[g,2]<-sum(groupdat2[,15+g]==1)
            runsum2[g,3]<-sum(groupdat2[,15+g]==0)
            runsum2[g,4]<-(sum(groupdat2[,15+g]==1)+sum(groupdat2[,15+g]==0)+1)
            runsum2[g,5]<-max(-skipvec2)
          }
          runsum2<-runsum2[which(runsum2[,2]==max(runsum2[,2])),] #choose w most ones
          runsum2<-runsum2[which(runsum2[,3]==min(runsum2[,3])),] #choose w least 0s
          runsum2<-runsum2[which(runsum2[,5]==min(runsum2[,5])),] #choose w smallest maximum skip (most gradual)
          runsum2<-runsum2[which(runsum2[,4]==min(runsum2[,4])),] #choose w least length
          runsum2<-runsum2[1,] #choose first one
          
          downsweep<-groupdat2[,c(1:15,15+runsum2[,1])]
          downsweep<-downsweep[which(downsweep[,16]==2|downsweep[,16]==1),]
          downsweep<-downsweep[,c(1:15)]
          upsweep<-groupdat[,c(1:15,15+runsum[,1])]
          upsweep<-upsweep[which(upsweep[,16]==2|upsweep[,16]==1),]
          upsweep<-upsweep[,c(1:15)]
          if(mean(downsweep[,15])<mean(upsweep[,15])&downsweep[1,13]<upsweep[nrow(upsweep),13]&downsweep[nrow(downsweep),15]==upsweep[1,15]){
            groupdat<-rbind(downsweep,upsweep[c(2:nrow(upsweep)),])
            kill="s"
            unwantedSelections<-c(unwantedSelections,downsweep[which(downsweep[,16]!=2&downsweep[,16]!=1),1])
            unwantedSelections<-c(unwantedSelections,upsweep[which(upsweep[,16]!=1),1])
          }else if(runsum2[,2]>=runsum[,2]+downsweepCompAdjust){
            kill="y"
          }else{
            kill="n"
          }
        }
        
        if(kill=="n"){
          groupdat<-groupdat[,c(1:15,15+runsum[,1])]
          unwantedSelections<-c(unwantedSelections,groupdat[which(groupdat[,16]!=2&groupdat[,16]!=1),1])
        }else if(kill=="y"){
          unwantedSelections<-c(unwantedSelections,groupdat[,1])
        }
      }
      
      resltsTSPV<-resltsTSPV[-which(resltsTSPV[,1] %in% unwantedSelections),]
      
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
drivepath<-"E:/"
startcombpath<-paste(drivepath,"Combined_sound_files/",sep="")
BLEDpath<-"C:/Users/danby456/Raven Pro 1.5/Presets/Detector/Band Limited Energy Detector/"
ravenpath<-"C:/Users/danby456/Raven Pro 1.5"
outputpath<-paste(drivepath,"DetectorRunOutput/",sep="")
processedGTpath<-paste(drivepath,"Processed_GT_data/",sep="")

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

##########sections to run
runGT<-"y"
runTestModel<-"n" #run model on GT data
runNewData<-"n" #run on data that has not been ground truthed. 

#enter the run name:
runname<- "break apart script test "

#Run type: all (all) or specific (spf) moorings to run
runtype<-"all"

#enter the detector type: "spread" or "single" or "combined". Can run and combine any combination of spread and single detectors that will be averaged after returning their detections. 
dettype<- "spread" 

#for newData: enter the type of mooring you'd like to analyze data: high graded (HG) or on full mooring (FULL)
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
detectorssin <- c(dir(BLEDpath)[2]
                  ) #list single detectors to run 
detectorssinshort<- detectorssin
}
##################Sampling rate of sound files

samplingRate<-16384 #hz

##########################max min length parameters (applies on R final detections, can also change in Raven to change initial box size)

Maxdur<-3.5
Mindur<-0.2

############################Combine detector  parameters
timediffself<-1.25
probdist<-.2 #how apart the probabilities can be before only choosing the best one. If within probdist of each other, combine them and average probability

#multiple detectors
freqdiff<-100
timediff<-1

#############################compare with downsweeps parameters
downsweepCompMod<-2
downsweepCompAdjust<-(3)

############################Whiten parameters (need to have done this in Raven previously)

#Pre whiten data?(y or no)
whiten<-"y"
FO<-100 #filter order
LMS<-.10 #LMS step size

###########################Model paramaters
CV=30 #number of cross validation models to create (more consistent results with greater number of models)
TPRthresh=.95 #estimated number of TP's to retain from total that initially go into model. Results can vary if data does not resemble GT data. 

#context sim parameters
greatcallThresh<-0.8 #level at which rolling score will reset to greatcallLevel
maxBonus<-0.05 #Maximum bonus level. Same as reset value when greatcallThresh is reached. 

goodcallBonus<-0.1 #increase factor for every call within resonable range (-maxPenalty to greatcallThresh)

maxPenalty<-(-0.35) #same abs value as lowest value for "good calls" that give bonus 
badcallPenalty<-(-0.001) #constant decrease every detection below threshold

############################Spread parameters. must be same length as number of spread detectors you are running
#p7 good ones: 2,1,3,0.75
#p9 working ones: 3,2,3,.25
#p10 good ones: 3,2,4,0.5
#(SPREAD) enter the desired smallest sequence size for detection. 
grpsize<-c(4)

#(SPREAD) allowed consecutive descending boxes allowed to still constitute an ascending sequence. Will end sequence after the maximum has been exceeded
allowedZeros<-c(2)

#(SPREAD) threshold of how many detectors at most can be skipped to be counted as sequential increase. 
detskip<-c(3)

#(SPREAD) max time distance for detectors to be considered in like group 
groupInt<-c(0.45)

############################file combine parameters
fileCombinesize<-345

downsample<-"y"
setSampRate<-4000

############################

runname<-paste(runname,gsub("\\D","",Sys.time()),sep="_")
dir.create(paste(outputpath,runname,sep=""))

if(runtype=="all"){
moorings<- colnames(MooringsDat)
#SF<-allmooringsSF
}else{
  allmooringsGT<- c("BS13_AU_04") #add as complete GTs 
  allmooringsSF<-list()#list sound file range for comleted GT of each mooring 
  allmooringsSF[[1]]<-c(1,304)
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

if(spec=="RW"){
  ravenView<-"RW_Upcalls"
}else if(spec=="GS"){
  ravenView<-"RW_GS"
}

if(runGT=="y"){
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
    resltVarInt <- raven_batch_detec(raven.path = ravenpath, sound.files = combname, path = paste(startcombpath,spec,"/",whiten2,sep=""),detector = "Band Limited Energy Detector",dpreset=i,vpreset=ravenView)
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
      resltVar <- raven_batch_detec(raven.path = ravenpath, sound.files = combname, path = paste(startcombpath,spec,"/",whiten2,sep=""),detector = "Band Limited Energy Detector",dpreset=r,vpreset=ravenView)
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
    resltVar <- raven_batch_detec(raven.path = ravenpath, sound.files =  combname, path = paste(startcombpath,spec,"/",whiten2,sep=""),detector = "Band Limited Energy Detector",dpreset=n,vpreset =ravenView)
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
colClasses = c("character","character","character","character","character","numeric","numeric","numeric", "numeric","numeric","numeric","numeric","character","character","character","character","character","character","character","character","character","character","character","character","numeric","numeric","character")
detecEvalFinal <- read.csv(text="Species, Moorings, Detectors, DetType, RunName, numTP, numFP, numFN, TPhitRate, TPR, TPdivFP,AUCav,CV_TPRthresh,Greatcall_goodcall,Max_modifier_penalty,ZerosAllowed,GroupSize,DownsweepThresh_DownsweepDiff,SkipAllowance,GroupInterval,TimeDiff,TimeDiffself,MinMaxDur,numDetectors,FO,LMS,Notes", colClasses = colClasses)

  GTtot=0
  GTtot2=NULL
  TPtot=NULL
  MoorCor=NULL
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
  GTtot2<-c(GTtot2,nrow(GT[[v]]))
  
  TPhitRate <- numTP/numTPtruth*100
  TPR <- numTP/(numTP+numFN)
  TPdivFP<- numTP/numFP

  #save stats and parameters to excel file
  detecEval<-detecEvalFinal[0,]
  if(dettype=="spread"|dettype=="combined"){
    detecEval[1,]<-c(spec,sort(unique(DetecTab2$Mooring))[v],paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,numTP,numFP,numFN,TPhitRate,TPR,TPdivFP,NA,NA,NA,NA,paste(allowedZeros,collapse=","),paste(grpsize,collapse=","),paste(downsweepCompMod,downsweepCompAdjust,sep=","),paste(detskip,collapse=","),paste(groupInt,collapse=","),timediff,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")
  }else{
    detecEval[1,]<-c(spec,sort(unique(DetecTab2$Mooring))[v],paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,numTP,numFP,numFN,TPhitRate,TPR,TPdivFP,NA,NA,NA,NA,NA,NA,NA,NA,NA,timediff,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")   
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
  detecEval[1,]<-c(spec,"all",paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,numTP,numFP,numFN,TPhitRate,TPR,TPdivFP,NA,NA,NA,NA,paste(allowedZeros,collapse=","),paste(grpsize,collapse=","),paste(downsweepCompMod,downsweepCompAdjust,sep=","),paste(detskip,collapse=","),paste(groupInt,collapse=","),timediff,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")
  }else{
detecEval[1,]<-c(spec,"all",paste(detnum,paste(detlist2,collapse="+"),sep=";"),dettype,runname,numTP,numFP,numFN,TPhitRate,TPR,TPdivFP,NA,NA,NA,NA,NA,NA,NA,NA,NA,timediff,timediffself,paste(Mindur,Maxdur,sep=","),as.character(paste(detnum,sum(detlist),sep=";")),FO,LMS," ")   
  }
detecEvalFinal <- rbind(detecEvalFinal,detecEval)
 
detecEval2<-read.csv(paste(outputpath,"DetectorRunLog.csv",sep=""))
detecEvalFinal<-rbind(detecEval2,detecEvalFinal)

write.csv(detecEvalFinal,paste(outputpath,"DetectorRunLog.csv",sep=""),row.names=FALSE)

beep(10)
###################
#MAYBE TEMPORARY- Save dataset for next steps so don't have to rerun after crash

#write.csv(DetecTab2,paste(outputpath,"DetecTab2.csv",sep=""),row.names=FALSE)


###################

###############################################


#RUN RANDOM FOREST MODEL 


#################

################################################
runname<-runname
spec<-spec
#Which data would you like to evaluate?
#species
#yminn<-0 #for spec plotting. Should be the same as detector window preset
#ymaxx<-1000 #" "

#define this table to compare counts after running model
TPtottab<-data.frame(TPtot,GTtot2,MoorCor)

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


#######1 mooring test######
#data<-data[which(data$`soundfiles[n]`=="BS15_AU_02a_files1-104.wav"),]

#add frequency stats to data 
data$meanfreq<- (data$Low.Freq..Hz.+data$High.Freq..Hz.)/2
data$freqrange<- (data$High.Freq..Hz.-data$Low.Freq..Hz.)
data$meantime<- (data$Begin.Time..s.+data$End.Time..s.)/2

#make interference columns into factors
if(length(data)>12){
  for(n in 13:length(data)){
    data[,n]<-as.factor(data[,n])
  }
}
data$Selection<-seq(1,nrow(data))

#data<-splitdf(data,weight = 1/4)[[1]]

data<-spectral_features(data,1)
write.csv(data,paste(processedGTpath,runname,"_processedGT.csv"),row.names = F)

}else{
  recentTab<-file.info(list.files(processedGTpath, full.names = T))
  recentPath<-rownames(recentTab)[which.max(recentTab$mtime)]
  data<-read.csv(recentPath) #produces most recently modifed file 
}

if(runTestModel=="y"){
  
data2<-data[,c(1,2,5,6,7,8,9:length(data))]
data2$detectionType<-as.factor(data2$detectionType)
#names(data2)[1]<-"Mooring"

my.xval = list()
my.xval$predictions = list()
my.xval$labels = list()

AUC_avg<-c()
f=1
CUTvec=NULL
#how many cross validation runs you want
for(p in 1:CV){
  print(paste("model",p))
  train<-splitdf(data2,weight = 2/3)
  data.rf<-randomForest(formula=detectionType ~ . -Selection -`soundfiles[n]`-meantime -Begin.Time..s. -End.Time..s. -Low.Freq..Hz. -High.Freq..Hz.,data=train[[1]],mtry=7,na.action=na.roughfix)
  pred<-predict(data.rf,train[[2]],type="prob")
  pred<-cbind(pred,train[[2]]$Selection)
  ROCRpred<-prediction(pred[,2],train[[2]]$detectionType)
  auc.perf = performance(ROCRpred, measure = "auc",plot=F)
  AUC_avg<-c(AUC_avg,as.numeric(auc.perf@y.values))
  
  my.xval$predictions[[p]] = ROCRpred@predictions[[1]]
  my.xval$labels[[p]] = ROCRpred@labels[[1]]
  
  prob.perf = performance(ROCRpred, "tpr","fpr")
  
  TPR<-NULL
  TPR <- data.frame(cut=prob.perf@alpha.values[[1]], tpr=prob.perf@y.values[[1]])
  CUT <- max(TPR[which(TPR$tpr>=TPRthresh),1])
  CUTvec<-c(CUTvec,CUT)
  
  if(f==1){
    probstab<-data.frame(data2$Selection)
    probstab[,f+1]<-NA
    for(n in 1:nrow(pred)){
      probstab[which(probstab$data2.Selection==pred[n,3]),f+1]<-pred[n,2]
    }
  }else{
    probstab[,f+1]<-NA
    for(n in 1:nrow(pred)){
      probstab[which(probstab$data2.Selection==pred[n,3]),f+1]<-pred[n,2]
    }
  }  
  f=f+1
}

##Graph std error and probability rates, with true detection included 
probmean<-NULL
probstderr<-NULL
n<-NULL
for(x in 1:nrow(probstab)){
  probmean[x]<-mean(as.numeric(probstab[x,2:length(probstab)]),na.rm=TRUE)
  probstderr[x]<-std.error(as.numeric(probstab[x,2:length(probstab)]),na.rm=TRUE)
  n[x]<-length(which(is.na(probstab[x,2:length(probstab)])))
}

CUTmean<-mean(CUTvec)
CUTstd.err<-std.error(CUTvec)

data3<- data2
data4<- data2

##assuming $detection type is already in this data NOTE not 
data3$probmean<-probmean
data3$probstderr<-probstderr
data3$n<-n

######################
#adaptively combine detections based on probability
data3<-adaptive_compare(data3,1)

#simulate context over time using probability scores 
data3<-context_sim(data3)

#number of TPs in data3
finTPs<-sum(as.numeric(as.character(data3[which(data3$probmean>CUTmean),]$detectionType)))

finFPs<-(nrow(data3[which(data3$probmean>CUTmean),])-finTPs)
finRat<-finTPs/finFPs

pp2<-data3$probmean
ll2<-data3$detectionType
predd2<-prediction(pp2,ll2)
perff2<-performance(predd2,"tpr","fpr")

varImpPlot(data.rf,  
           sort = TRUE,
           main="Variable Importance")

#with permutations on probs
plot(perff2, avg = "threshold",  xaxs="i", yaxs="i", spread.scale=2,
     lwd = 2, main = paste("Threshold avg"),colorize=T)
abline(a=0, b= 1)
auc.perf = performance(predd2, measure = "auc",plot=F)
print(auc.perf@y.values)

AUCadj<-auc.perf@y.values

plot(data3$probmean,data3$probstderr, col = ifelse(data3$detectionType==1,'blue','red'),cex=0.25)
abline(v=CUTmean)

#see freq breakdown of calls 
cdplot(data3$detectionType ~ data3$meanfreq, data3, col=c("cornflowerblue", "orange"), main="Conditional density plot")
cdplot(data3$detectionType ~ data3$freqrange, data3, col=c("cornflowerblue", "orange"), main="Conditional density plot")
cdplot(data3$detectionType ~ data3$specprop.mode, data3, col=c("cornflowerblue", "orange"), main="Conditional density plot")
cdplot(data3$detectionType ~ data3$meanpeakf, data3, col=c("cornflowerblue", "orange"), main="Conditional density plot")
cdplot(data3$detectionType ~ data3$High.Freq..Hz., data3, col=c("cornflowerblue", "orange"), main="Conditional density plot")
cdplot(data3$detectionType ~ data3$Low.Freq..Hz., data3, col=c("cornflowerblue", "orange"), main="Conditional density plot")


#write data to drive
after_model_write(data3,1)

beep(10)



#comparison dataset to data3#####################################
data4$probmean<-probmean
data4$probstderr<-probstderr
data4$n<-n

#number of TPs in data4
sum(as.numeric(as.character(data4[which(data4$probmean>CUTmean),]$detectionType)))

#cor.test(as.numeric(probmean),probstderr)
pp = my.xval$predictions
ll = my.xval$labels
predd = prediction(pp, ll)
perff = performance(predd, "tpr", "fpr")

#no avg
#plot(perff, xaxs="i", yaxs="i",main=paste("All",CV," cross validation runs"))
#abline(a=0, b= 1)

#avg. Don't really know what the points on the line mean. Without manipulations on probs
#plot(perff, avg = "vertical", spread.estimate = "stddev",spread.scale=2, xaxs="i", yaxs="i", 
     #show.spread.at=c(.05,.075,.1,.125,.15,.2,.3),
#     lwd = 2, main = paste("Vertical avg w/ std devn"))
#plot(perff, avg = "threshold",  xaxs="i", yaxs="i", spread.scale=2,
#     lwd = 2, main = paste("Threshold avg"),colorize=T)
#abline(a=0, b= 1)
#print(mean(AUC_avg))

#save last model
save(data.rf, file = paste("E:/DetectorRunOutput/",runname,"/an_example_model.rda",sep=""))

}

###############################################


#EMPLOY MODEL ON FULL DATASETS

################################################

if(runNewData=="y"){

allDataPath<-"E:/Datasets"
allMoorings<-dir(allDataPath)[1] #Just AW12_AU_BS3 right now, need fully analyzed GT to test on full mooring

if(whiten!="y"){
fileSizeInt<-(fileCombinesize*4)
}else{
fileSizeInt<-(fileCombinesize*4*3) #whitened files are smaller so still under 6 gigs. 
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

#decimate dataset. need to add txt doc to say whether decimated already or not or will continue to decimate. right now, decimating 1000 files per hour (17,000 files is AW15_AU_BS3 mooring)
if(downsample=="y"){
  ravenView<-paste(ravenView,"_",setSampRate,"DS",sep="")
    if(!file.exists(paste(sfpath,"downsamp_to",setSampRate,sep="_"))){
      dir.create(paste(sfpath,"downsamp_to",setSampRate,sep="_"))
      sfpath2<-sfpath
      sfpath<-paste(sfpath,"downsamp_to",setSampRate,sep="_")
      for(z in dir(sfpath2)){
        wav<-readWave(paste(sfpath2,"/",z,sep=""),unit="sample")
        wav.samp.rate<-wav@samp.rate
        wav<-downsample(wav,setSampRate)
        writeWave(wav, filename=paste(sfpath,"/",z,sep=""), extensible = FALSE)
      }
      write.table(paste("This data has been decimated by factor of",decimationFactor),paste(sfpath,"/decimationStatus"),quote=FALSE,sep = "\t",row.names=FALSE,col.names=FALSE)
      
    }
}

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
        sox_alt(paste(noquote(paste(paste(sound_filesfullpath[1],collapse=" ")," ",combSound,sep=""))),exename="sox.exe",path2exe="E:\\Accessory\\sox-14-4-2")
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
        resltVarInt <- raven_batch_detec(raven.path = ravenpath, sound.files = combname, path =filePath,detector = "Band Limited Energy Detector",dpreset=i,vpreset=ravenView)
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
          resltVar <- raven_batch_detec(raven.path = ravenpath, sound.files = combname, path = filePath,detector = "Band Limited Energy Detector",dpreset=r,vpreset=ravenView)
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
        resltVar <- raven_batch_detec(raven.path = ravenpath, sound.files =  combname, path = filePath,detector = "Band Limited Energy Detector",dpreset=n,vpreset =ravenView)
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
colClasses = c("character","character","character","character","character","numeric","numeric","numeric", "numeric","numeric","numeric","numeric","numeric","character","character","character","character","character","character","character","character","character","character","character","character","numeric","numeric","character")
detecEvalFinal <- read.csv(text="Species, Moorings, Detectors, DetType, RunName, numTP, numFP, numFN, TPhitRate, TPR, TPdivFP,AUCav,CV_TPRthresh,Greatcall_goodcall,Max_modifier_penalty,ZerosAllowed,GroupSize,DownsweepThresh_DownsweepDiff,SkipAllowance,GroupInterval,TimeDiff,TimeDiffself,MinMaxDur,numDetectors,FO,LMS,Notes", colClasses = colClasses)


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
#combine or eliminate detections within timediffself parameter
findata<-adaptive_compare(findata,2)

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

}
  