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

boxplot(dataGT$AUC ~ dataGT$Test, data = dataGT, lwd = 2, ylab = 'AUC',ylim = c(min(dataNEW$AUC)-0.01,1),xlab="# moorings in GT set",main="GS")
#stripchart(dataALL$AUC ~ dataALL$Test, vertical = TRUE, data = dataALL, 
#           add = TRUE, pch = 12, col = 'blue',lwd=3)
#stripchart(dataNEW$AUC ~ dataNEW$Test, vertical = TRUE, data = dataNEW, 
   #         add = TRUE, pch = 18, col = 'red',lwd=4)
boxplot(dataNEW$AUC ~ dataNEW$Test, data = dataNEW, col=makeTransparent("white",0),boxcol=makeTransparent("red",99),medcol=makeTransparent("red",99),staplecol=makeTransparent("red",99),whiskcol=makeTransparent("red",99),outcol=makeTransparent("red",99), xaxt='n',lwd = 2,add=TRUE,at=2:(length(unique(dataNEW$Test))+1))


boxplot(dataGT$TPR ~ dataGT$Test, data = dataGT, lwd = 2, ylab = 'TPR',ylim = c(min(dataNEW$TPR)-0.01,1.3),xlab="# moorings in GT set",main="GS")
#stripchart(dataALL$TPR ~ dataALL$Test, vertical = TRUE, data = dataALL, 
#           add = TRUE, pch = 12, col = 'blue',lwd=3)
#stripchart(dataNEW$TPR ~ dataNEW$Test, vertical = TRUE, data = dataNEW, 
#            add = TRUE, pch = 18, col = 'red',lwd=4)
boxplot(dataNEW$TPR ~ dataNEW$Test, data = dataNEW, col=makeTransparent("white",0),boxcol=makeTransparent("red",99),medcol=makeTransparent("red",99),staplecol=makeTransparent("red",99),whiskcol=makeTransparent("red",99),outcol=makeTransparent("red",99), xaxt='n',lwd = 2,add=TRUE,at=2:(length(unique(dataNEW$Test))+1))



boxplot(dataGT$FPR ~ dataGT$Test, data = dataGT, lwd = 2, ylab = 'FPR',ylim = c(min(dataNEW$FPR)-0.01,1),xlab="# moorings in GT set",main="GS")
boxplot(dataNEW$FPR ~ dataNEW$Test, data = dataNEW, col=makeTransparent("white",0),boxcol=makeTransparent("red",99),medcol=makeTransparent("red",99),staplecol=makeTransparent("red",99),whiskcol=makeTransparent("red",99),outcol=makeTransparent("red",99), xaxt='n',lwd = 2,add=TRUE,at=2:(length(unique(dataNEW$Test))+1))

boxplot(dataALL$DeltGini ~ dataALL$Test, data = dataALL, lwd = 2, ylab = 'Delta Gini score',xlab="# moorings in GT set",main="GS")


stripchart(dataALL$FPR ~ dataALL$Test, vertical = TRUE, data = dataALL, 
           add = TRUE, pch = 12, col = 'blue',lwd=3)
stripchart(dataNEW$FPR ~ dataNEW$Test, vertical = TRUE, data = dataNEW, 
            add = TRUE, pch = 18, col = 'red',lwd=4)


