#install.packages("flightcallr")
#install.packages("randomForest")
#install.packages("seewave")
#install.packages("tuneR")
#install.packages("plotrix")
#install.packages("aod")
#install.packages("ggplot2")
#install.packages("usdm")
#install.packages("ROCR")

library(randomForest)
library(seewave)
library(tuneR)
library(plotrix)
library(aod)
library(ggplot2)
library(usdm)
library(flightcallr)
library(ROCR)


#Which data would you like to evaluate?
#species
samplingRate<-16384 #hz
yminn<-0 #for spec plotting. Should be the same as detector window preset
ymaxx<-1000 #" "

spec<-"RW"

runname<-"full run no minmax test_20181006144510"



##########
detfiles<-list.files(paste("E:/DetectorRunOutput/",runname,sep=""),pattern = "TPFPFN")  

#extract mooring names from moorings used in run
mooringpat=NULL
for(n in 1:length(detfiles)){
  mooringpat<-c(mooringpat,substr(detfiles[n],1,11))
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

data<-splitdf(data,weight = 1/4)[[1]]

for(z in 1:nrow(data)){
  foo <- readWave(paste("E:/Combined_sound_files/",spec,"/",soundfile,"/",data[z,1],sep=""),data[z,5],data[z,6],units="seconds")
  foo.spec <- spec(foo, plot=F, PSD=T,ylim=c(yminn,ymaxx))
  foo.specprop <- specprop(foo.spec)
  foo.meanspec = meanspec(foo, plot=FALSE, ovlp=90)#not sure what ovlp parameter does but initially set to 90
  foo.autoc = autoc(foo, plot=F)
  foo.dfreq = dfreq(foo, plot=F, ovlp=90)
  data$rugosity[z] = rugo(foo@left / max(foo@left)) #no idea how these @s work
  data$crest.factor[z] = crest(foo)$C
  foo.env = seewave:::env(foo, plot=F) 
  data$temporal.entropy[z] = th(foo.env)
  data$shannon.entropy[z] = sh(foo.spec)
  data$spectral.flatness.measure[z] = sfm(foo.spec)
  data$spectrum.roughness[z] = roughness(foo.meanspec[,2])
  data$autoc.mean[z] = mean(foo.autoc[,2], na.rm=T)
  data$autoc.median[z] = median(foo.autoc[,2], na.rm=T)
  data$autoc.se[z] = std.error(foo.autoc[,2], na.rm=T)
  data$dfreq.mean[z] = mean(foo.dfreq[,2], na.rm=T)
  data$dfreq.se[z] = std.error(foo.dfreq[,2], na.rm=T)
  data$specprop.mean[z] = foo.specprop$mean[1]
  data$specprop.sd[z] = foo.specprop$sd[1]
  data$specprop.sem[z] = foo.specprop$sem[1]
  data$specprop.median[z] = foo.specprop$median[1]
  data$specprop.mode[z] = foo.specprop$mode[1]
  data$specprop.Q25[z] = foo.specprop$Q25[1]
  data$specprop.Q75[z] = foo.specprop$Q75[1]
  data$specprop.IQR[z] = foo.specprop$IQR[1]
  data$specprop.cent[z] = foo.specprop$cent[1]
  data$specprop.skewness[z] = foo.specprop$skewness[1]
  data$specprop.kurtosis[z] = foo.specprop$kurtosis[1]
  data$specprop.sfm[z] = foo.specprop$sfm[1]
  data$specprop.sh[z] = foo.specprop$sh[1]
  print(paste("done with",z))
}

data<-data[,9:length(data)]

#data1<-data.frame(scale(data))

#test for collinearity
vif(data1)
#first<-vif(data)[which(vif(data)$VIF>4),1]
#second<-vif(data)[which(vif(data)$VIF>4),1]
#third<-vif(data)[which(vif(data)$VIF>4)
#fourth<-vif(data)[which(vif(data)$VIF>4),1]
new<-vif(data)[which(vif(data)$VIF>4),1]
#vif >4
#of original params: take one of (8,9),(13,16,18,19,21 keep: ),(14,20),(23,24),(6,24),(5,25). just take the first for now, remove the rest. 

#eigen(cor(data))$values
#kappa(cor(data),exact=T)
#cor(data)


#before going too deep into trying to tune the glm, give random forest a shot.
data2$detectionType<-factor(data2$detectionType)
data2<-data[,1:(length(data)-2)]
train<-splitdf(data2,weight = 2/3)


data.rf<-randomForest(x=train[[1]],formula=detectionType ~ .)
data.rf




#pairs(data2, upper.panel = NULL)


mylogit<-glm(formula=detectionType~rugosity+crest.factor+temporal.entropy+shannon.entropy+spectral.flatness.measure+
               spectrum.roughness+autoc.mean+autoc.se+dfreq.mean+dfreq.se+specprop.mean+specprop.sd+specprop.sem+
               specprop.mode,family="binomial",data=train[[1]])
summary(mylogit)

pred<-predict(mylogit,train[[2]],type="response")
model_pred_det<-rep("0",nrow(train[[2]]))
model_pred_det[pred>0.135]<-"1"
tab<-table(model_pred_det,train[[2]]$detectionType)
print(tab)
1-sum(diag(tab))/sum(tab) #15.5% misclassification with pred>.5 (first 6000 rows)


ROCRpred<-prediction(as.numeric(model_pred_det),train[[2]]$detectionType)

roc.perf = performance(ROCRpred, measure = "tpr", x.measure = "fpr")
plot(roc.perf)
abline(a=0, b= 1)

#going with this approach, should include site and season as fixed effects, or include site/season (mooring name) as a random effect and run a glmm
######################################




