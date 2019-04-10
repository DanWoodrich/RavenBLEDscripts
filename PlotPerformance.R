#plot boxplots showing difference in GT data and NEW data

data<-read.csv("E:/RavenBLEDscripts/Data/Paper/AddDataSummaryGS.csv")

dataGT<-data[which(data$Type=="GT"),]
dataNEW<-data[which(data$Type=="NEW"),]
dataALL<-data[which(data$Type=="ALL"),]


boxplot(dataGT$AUC ~ dataGT$Test, data = dataGT, lwd = 2, ylab = 'AUC',ylim = c(min(dataNEW$AUC)-0.01,1),xlab="Experiment",main="")
stripchart(dataALL$AUC ~ dataALL$Test, vertical = TRUE, data = dataALL, 
           add = TRUE, pch = 12, col = 'blue',lwd=3)
stripchart(dataNEW$AUC ~ dataNEW$Test, vertical = TRUE, data = dataNEW, 
            add = TRUE, pch = 18, col = 'red',lwd=4)



boxplot(dataGT$TPR ~ dataGT$Test, data = dataGT, lwd = 2, ylab = 'TPR',ylim = c(min(dataNEW$TPR)-0.01,1),xlab="Experiment")
stripchart(dataALL$TPR ~ dataALL$Test, vertical = TRUE, data = dataALL, 
           add = TRUE, pch = 12, col = 'blue',lwd=3)
stripchart(dataNEW$TPR ~ dataNEW$Test, vertical = TRUE, data = dataNEW, 
            add = TRUE, pch = 18, col = 'red',lwd=4)



boxplot(dataGT$FPR ~ dataGT$Test, data = dataGT, lwd = 2, ylab = 'FPR',ylim = c(min(dataNEW$FPR)-0.01,1),xlab="Experiment")
stripchart(dataALL$FPR ~ dataALL$Test, vertical = TRUE, data = dataALL, 
           add = TRUE, pch = 12, col = 'blue',lwd=3)
stripchart(dataNEW$FPR ~ dataNEW$Test, vertical = TRUE, data = dataNEW, 
            add = TRUE, pch = 18, col = 'red',lwd=4)


