#random forest model to seperate FPs and TPs from Raven output

#install.packages("flightcallr")
#install.packages("randomForest")
#install.packages("seewave")
#install.packages("tuneR")
#install.packages("plotrix")
#install.packages("aod")
#install.packages("ggplot2")
#install.packages("usdm")
#install.packages("ROCR")
#install.packages("e1071")  
#install.packages("caret")  

library(e1071)  
library(caret)  
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


#Which data would you like to evaluate?
#species
samplingRate<-16384 #hz
yminn<-0 #for spec plotting. Should be the same as detector window preset
ymaxx<-1000 #" "

spec<-"RW"

runname<-"mooring and pulse full test_20181010145447"



##########
detfiles<-list.files(paste("E:/DetectorRunOutput/",runname,sep=""),pattern = "RF")  

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

#make interference columns into factors
if(length(data)>8){
  for(n in 9:length(data)){
    data[,n]<-as.factor(data[,n])
  }
}
#######1 mooring test######
#data<-data[which(data$`soundfiles[n]`=="BS15_AU_02a_files1-104.wav"),]

data<-splitdf(data,weight = 1/3)[[1]]

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

data2<-data[,9:length(data)]
colnames(data2)[1]<-"File"

#data1<-data.frame(scale(data))

#test for collinearity
#vif(data1)
#first<-vif(data)[which(vif(data)$VIF>4),1]
#second<-vif(data)[which(vif(data)$VIF>4),1]
#third<-vif(data)[which(vif(data)$VIF>4)
#fourth<-vif(data)[which(vif(data)$VIF>4),1]
#new<-vif(data)[which(vif(data)$VIF>4),1]
#vif >4
#of original params: take one of (8,9),(13,16,18,19,21 keep: ),(14,20),(23,24),(6,24),(5,25). just take the first for now, remove the rest. 

#eigen(cor(data))$values
#kappa(cor(data),exact=T)
#cor(data)

#########################################
#before going too deep into trying to tune the glm, give random forest a shot.
#test individual moorings to see if it makes it better

data2$detectionType<-as.factor(data2$detectionType)
train<-splitdf(data2,weight = 2/3)

#test using file as variable
#randfile<- sample(1:7, 4, replace=F)
#train<-list()
#train[[1]]<-data2[which(data2$File %in% unique(data2$File)[randfile]),]
#train[[2]]<-data2[which(data2$File %in% unique(data2$File)[-randfile]),]
#test out what happens if I remove most of the 0s. Did not seem to have any positive effect. 
#train2<-splitdf(train[[1]][which(train[[1]]$detectionType==0),],weight=1/6)[[1]]
#train3<-train[[1]][which(train[[1]]$detectionType==1),]
#train2<-rbind(train2,train3)
#train3<-NULL


mtry_avg<-c()
mtry_count<-c()
for(i in seq(100,3000,100)){
  print(paste("ntree value",i))
  AUC_avg<-c()
  for(p in 1:20){
    print(paste("model",p))
    train<-splitdf(data2,weight = 2/3)
    data.rf<-randomForest(formula=detectionType ~ .,data=train[[1]],mtry=7,ntree=i)
    pred<-predict(data.rf,train[[2]],type="prob")
    ROCRpred<-prediction(pred[,2],train[[2]]$detectionType)
    auc.perf = performance(ROCRpred, measure = "auc",plot=F)
    AUC_avg<-c(AUC_avg,as.numeric(auc.perf@y.values))
      }
  mtry_avg<-c(mtry_avg,mean(AUC_avg))
  mtry_count<-c(mtry_count,i)
}

plot(mtry_count,mtry_avg)
lm(formula = mtry_avg ~ mtry_count)
abline(8.23e-01,  3.48e-07 )
#appears that mtry at 7 is consistently best (when ntree is not specified). 
#when mtry is 7, looks like ntree is pretty much neutral from 400 onwards...? Seems like a very low # of trees is worse (~100), but it is the same after that. Say 500 for ease of computation. 

#cross validation:
#data.rf.cv<-rfcv(train[[1]][,2:length(train[[1]])],train[[1]]$detectionType,20,step=0.8,formula=detectionType ~ .)
#data.rf.cv

#smallest error at step .8, nfold 20, at 7th iteration. Only 6 variables really used for rf?? 
#top 6 variables: specprop.mode,autoc.median,temporal.entropy,spectrum.roughness,specprop.sem,specprop.kurtosis

#test<-data.rf.cv$predicted[7]




#test against test data
pred<-predict(data.rf,train[[2]],type="prob")

ROCRpred<-prediction(pred[,2],train[[2]]$detectionType)

roc.perf = performance(ROCRpred, measure = "tpr", x.measure = "fpr")
plot(roc.perf)
abline(a=0, b= 1)

#acc.perf = performance(ROCRpred, measure = "acc")
#plot(acc.perf)

auc.perf = performance(ROCRpred, measure = "auc",plot=F)
auc.perf@y.values

varImpPlot(data.rf,  
           sort = 27,
           n.var=27,
           main="Top 10 - Variable Importance")

#did better on only BS15_AU_02a- at .35 pred cutoff, achieved AUC score of .87,.82,.85,.87,.85... (small sample size)  
                                #larger sample size: .82,.86,.82,.82,.85,.87,.88 (all done without specifying mtry)
                    #AW12_AU_BS3- full samp size auc: .77,.79,.76,.78,.78,78,.79... correctly gets about 60% TPs. Overall including Raven BLED would be about .6*.83 = .5 of all calls (including manual review of yeses produced here, which adds another 50% to look at beyond TPs). # is also mooring specific, would be worse with a more general detector. 
                    #BS13_AU_04- full samp size auc: .71,(from here on with better looking auc curve, no cutoff..),.79,.799,.77...
                    #BS15_AU_02a with no cutoff: .94,.94.,.9,.9,.9,.92,.9,...
                    #full moorings AUC with no cutoff: .82,.815,.8,.815,.83,.8... about .85 sensitivity for 50% FP rate. meaning on average for full analysis, could get .8 from BLED, .85 from RF for a total of 68% of upcalls, and would have to manually review an equivalent amount of FPs as TPs.  
                    #full moorings AUC with half of 0s taken out:.8,.82
                    #full moorings AUC with only 1/6 of original no's:.8,.815,.8
                    #OOB % error rate: at 500 trees: 11.06,11.13
                                      #at 100 trees: 11.41%
                                      #at 300 trees: 11.32%
                                      #at 1000 trees: 11.16%
                                      #at 10000 trees: 11.09% Does not appear to be making a big difference

#After reducing it to 6 vars instead on 24 based on rfcv() analysis...:
#run with top 6 vars: 0.8339438, .8,.82,.8 - fairly similar to with all 24. Good to know if need to reduce computation time, otherwise will include all vars for after I introduce fin and mooring detectors.  
                                                
#test: use only BS15_AU_02a to train model, and use it to predict rest of moorings. result: AUC .73,.72... This could be an option, if you just want the "best" calls, would involve more manual analysis. 

#test: use binary fin and mooring noise detector as categorical variables. All other variables included. AUC: .845, .85, .815, .827, .81, .83. May be a slight improvement, mooring and fin detectors could be made better as well. With .83 AUC run could get .9 from random forest with 50% FP rate. 
#   The variable importance suggests that the mooring/fin and pulses detectors were the worst predictors, visaully by a pretty clear margin. 

#test: include file as a variable, randomly compare 3 moorings (train) to 4 (test). AUC: .75, .79, .75, .75,.76 (with 4 in train now) .7,.81,.75  ... seems to do nothing, as expected, when variable isn't repeated from the training to the test set. Could be informative if I choose to include more GT data that matches in part with each mooring. 

#could try with mooring, or with year- variables that represent a better 
############################################################

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


