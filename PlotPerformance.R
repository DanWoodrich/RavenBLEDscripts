#plot boxplots showing difference in GT data and NEW data

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

data<-read.csv("E:/DetectorRunOutput/AddDataSummaryGS_automate.csv")

dataGT<-data[which(data$Type=="GT"),]
dataNEW<-data[which(data$Type=="NEW"),]
dataALL<-data[which(data$Type=="ALL"),]

boxplot(dataGT$AUC ~ dataGT$Test, data = dataGT, lwd = 2, ylab = 'AUC',ylim = c(min(dataNEW$AUC)-0.01,1),xlab="# data segments in labelled data",main="Upcall classifier - scatter", cex.lab=1.2)
stripchart(dataNEW$AUC ~ dataNEW$Test, vertical = TRUE, data = dataNEW, 
            add = TRUE, pch = 18, col = 'red',lwd=4,at=2:(length(unique(dataNEW$Test))+1))
legend("bottomright", 
       legend = c("Labelled data", "Holdout data"), 
       col = c("black","red"), 
       pch = c(15,18), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))

boxplot(dataGT$AUC ~ dataGT$Test, data = dataGT, lwd = 2, ylab = 'AUC',ylim = c(min(dataNEW$AUC)-0.01,1),xlab="# data segments in labelled data",main="Upcall classifier - boxplot", cex.lab=1.2)
boxplot(dataNEW$AUC ~ dataNEW$Test, data = dataNEW, col=makeTransparent("white",0),boxcol=makeTransparent("red",99),medcol=makeTransparent("red",99),staplecol=makeTransparent("red",99),whiskcol=makeTransparent("red",99),outcol=makeTransparent("red",99), xaxt='n',lwd = 2,add=TRUE,at=2:(length(unique(dataNEW$Test))+1))
boxplot(dataNEW$AUC ~ dataNEW$Test, data = dataNEW, col=makeTransparent("white",0),boxcol=makeTransparent("red",99),medcol=makeTransparent("red",99),staplecol=makeTransparent("red",99),whiskcol=makeTransparent("red",99),outcol=makeTransparent("red",99), xaxt='n',lwd = 2,add=TRUE,at=2:(length(unique(dataNEW$Test))+1))
boxplot(dataNEW$AUC ~ dataNEW$Test, data = dataNEW, col=makeTransparent("white",0),boxcol=makeTransparent("red",99),medcol=makeTransparent("red",99),staplecol=makeTransparent("red",99),whiskcol=makeTransparent("red",99),outcol=makeTransparent("red",99), xaxt='n',lwd = 2,add=TRUE,at=2:(length(unique(dataNEW$Test))+1))
legend("bottomright", 
       legend = c("Labelled data", "Holdout data"), 
       col = c("black","red"), 
       pch = c(15,15), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))

#line<-lm(dataNEW$AUC ~ dataNEW$Test)
#abline(a=line$coefficients[1],b=line$coefficients[2],col="red")

boxplot(dataGT$FPR ~ dataGT$Test, data = dataGT, lwd = 2, ylab = 'FPR',ylim = c(min(dataNEW$FPR)-0.01,1),xlab="# data segments in labelled data",main="Upcall classifier - scatter",cex.lab=1.2)
stripchart(dataNEW$FPR ~ dataNEW$Test, vertical = TRUE, data = dataNEW, 
           add = TRUE, pch = 18, col = 'red',lwd=4,at=2:(length(unique(dataNEW$Test))+1))
legend("topright", 
       legend = c("Labelled data", "Holdout data"), 
       col = c("black","red"), 
       pch = c(15,18), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))

boxplot(dataGT$FPR ~ dataGT$Test, data = dataGT, lwd = 2, ylab = 'FPR',ylim = c(min(dataNEW$FPR)-0.01,1),xlab="# data segments in labelled data",main="Upcall classifier - boxplot", cex.lab=1.2)
boxplot(dataNEW$FPR ~ dataNEW$Test, data = dataNEW, col=makeTransparent("white",0),boxcol=makeTransparent("red",99),medcol=makeTransparent("red",99),staplecol=makeTransparent("red",99),whiskcol=makeTransparent("red",99),outcol=makeTransparent("red",99), xaxt='n',lwd = 2,add=TRUE,at=2:(length(unique(dataNEW$Test))+1))
boxplot(dataNEW$FPR ~ dataNEW$Test, data = dataNEW, col=makeTransparent("white",0),boxcol=makeTransparent("red",99),medcol=makeTransparent("red",99),staplecol=makeTransparent("red",99),whiskcol=makeTransparent("red",99),outcol=makeTransparent("red",99), xaxt='n',lwd = 2,add=TRUE,at=2:(length(unique(dataNEW$Test))+1))
boxplot(dataNEW$FPR ~ dataNEW$Test, data = dataNEW, col=makeTransparent("white",0),boxcol=makeTransparent("red",99),medcol=makeTransparent("red",99),staplecol=makeTransparent("red",99),whiskcol=makeTransparent("red",99),outcol=makeTransparent("red",99), xaxt='n',lwd = 2,add=TRUE,at=2:(length(unique(dataNEW$Test))+1))
legend("topright", 
       legend = c("Labelled data", "Holdout data"), 
       col = c("black","red"), 
       pch = c(15,15), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))








boxplot(dataNEW$FPR ~ dataNEW$Test, data = dataNEW, col=makeTransparent("white",0),boxcol=makeTransparent("red",99),medcol=makeTransparent("red",99),staplecol=makeTransparent("red",99),whiskcol=makeTransparent("red",99),outcol=makeTransparent("red",99), xaxt='n',lwd = 2,add=TRUE,at=2:(length(unique(dataNEW$Test))+1))

line<-lm(dataNEW$FPR ~ dataNEW$Test)
abline(a=line$coefficients[1],b=line$coefficients[2],col="red")


boxplot(dataGT$TPR ~ dataGT$Test, data = dataGT, lwd = 2, ylab = 'TPR',ylim = c(min(dataNEW$TPR)-0.01,1.3),xlab="# data segments in GT set",main="GS")
#stripchart(dataALL$TPR ~ dataALL$Test, vertical = TRUE, data = dataALL, 
#           add = TRUE, pch = 12, col = 'blue',lwd=3)
#stripchart(dataNEW$TPR ~ dataNEW$Test, vertical = TRUE, data = dataNEW, 
#            add = TRUE, pch = 18, col = 'red',lwd=4)
boxplot(dataNEW$TPR ~ dataNEW$Test, data = dataNEW, col=makeTransparent("white",0),boxcol=makeTransparent("red",99),medcol=makeTransparent("red",99),staplecol=makeTransparent("red",99),whiskcol=makeTransparent("red",99),outcol=makeTransparent("red",99), xaxt='n',lwd = 2,add=TRUE,at=2:(length(unique(dataNEW$Test))+1))




boxplot(dataALL$DeltGini ~ dataALL$Test, data = dataALL, lwd = 2, ylab = 'Delta Gini score',xlab="# data segments in GT set",main="GS")


stripchart(dataALL$FPR ~ dataALL$Test, vertical = TRUE, data = dataALL, 
           add = TRUE, pch = 12, col = 'blue',lwd=3)
stripchart(dataNEW$FPR ~ dataNEW$Test, vertical = TRUE, data = dataNEW, 
            add = TRUE, pch = 18, col = 'red',lwd=4)


