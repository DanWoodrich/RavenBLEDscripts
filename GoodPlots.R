#ROCnoColor
plot(perff2, avg = "threshold",  xaxs="i", yaxs="i", print.cutoffs.at=c(0.35,0.4,0.5,0.6,0.7,0.8,0.9,0.99),
     cutoff.label.function = function(x) return(paste("              ",x)),
     lwd = 2, main = paste("ROC for Average Upsweep Random Forest Classifier"),colorize=F)
abline(a=0, b= 1)

#giniscore
varImpPlot_AVG(giniAv,  
               sort = TRUE,
               main="Variable Importance Average Upsweep Random Forest Classifier")

plot(dataSPEC[which(dataSPEC$detectionType==paste(mSpec[s],0)),probIndex],dataSPEC[which(dataSPEC$detectionType==paste(mSpec[s],0)),varIndex], col = "red",cex=0.25,xlab="Average probability of true detection",ylab="Average standard error",main="Average Upsweep Random Forest Classifier Distribution Negatives")
abline(v=CUTmeanspec)
legend("topright", 
       legend = c("Negatives"), 
       col = c("red"), 
       pch = c(15), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1))

plot(dataSPEC[which(dataSPEC$detectionType==paste(mSpec[s],1)),probIndex],dataSPEC[which(dataSPEC$detectionType==paste(mSpec[s],1)),varIndex], col = "blue",cex=0.25,xlab="Average probability of true detection",ylab="Average standard error",main="Average Upsweep Random Forest Classifier Distribution Positives")
abline(v=CUTmeanspec)
legend("topright", 
       legend = c("Positives"), 
       col = c("blue"), 
       pch = c(15), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1))

plot(dataSPEC[,probIndex],dataSPEC[,varIndex], col = ifelse(dataSPEC$detectionType==paste(mSpec[s],1),'blue','red'),cex=0.25,xlab="Average probability of true detection",ylab="Average standard error",main="Average Upsweep Random Forest Classifier Distribution All")
abline(v=CUTmeanspec)
plot.new()
legend("right", 
       legend = c("Negatives", "Positives"), 
       col = c("red","blue"), 
       pch = c(15,15), 
       bty = "n", 
       pt.cex = 10, 
       cex = 6, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))

#plot of probabilities after context sim:
par(mgp=c(2.5,1,0)) 
par(cex.main=2,cex.lab=2,cex.axis=1.5)
#betternames<- c("AW15_AU_BS2","AW12_AU_BS3 1","BS14_AU_04","AW14_AU_BS3 1","AW15_AU_BS3","AW14_AU_BS3 2","BS12_AU_02a","BS12_AU_02b","AL16_AU_BS3","AW12_AU_BS3 2","BS13_AU_04")
betternames<-c("BS15_AU_02a","BS14_AU_04","AW12_AU_BS3","BS13_AU_04","BS16_AU_02a","BS15_AU_02b","AW14_AU_BS3","AL16_AU_BS1","BS13_AU_02a","BS16_AU_05","BS15_AU_04","AW14_AU_BS2")
for(m in unique(dataSPEC[,6])){
  data3Matmoors<-dataSPEC[which(dataSPEC[,6]==m),]
  CS_outputDate<-CS_output[which(CS_output[,1] %in% as.numeric(data3Matmoors[,1])),]
  pointsDate<-data3Matmoors$Begin.Time..s.
  #pointsDate<-as.POSIXlt((as.numeric(data3Matmoors$RTFb)+data3Matmoors$FileOffsetBegin), origin="1970-01-01",tz="UTC")
  plot(x=pointsDate,y=data3Matmoors[,nospec+1], col=c("red","blue")[as.numeric(data3Matmoors$detectionType)],cex=0.5,pch=19,xlab="Detection index (ordered by increasing time from DS start)",ylab="Average probability")
  title(betternames[which(unique(dataSPEC[,6])==m)],line=0.6,cex=2)
  abline(h=CUTmean[s],col="blue",lwd=2)
  #abline(h=0.5,lty=3)
  #lines(lowess(data3Matmoors[,pos+5]))
  #lines(lowess(data3Matmoors[,pos+4]))
  #lines(lowess(data3Matmoors[,pos+2]))
  #lines(pointsDate,(CS_outputDate[,5]),col="blue") #backwards through data 
  #lines(pointsDate,(CS_outputDate[,7]),col="orange") #forwards through data
  #lines(pointsDate,(CS_outputDate[,8]),col="green")
}

