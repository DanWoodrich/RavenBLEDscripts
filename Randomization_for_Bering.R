#choose randomized data 

#for RW

assTab<-NULL
for(n in 1:20){
  
assignment<-NULL

species<-c("RW","GS")

species<-sample(species,1)

assignment<-c(species)

Years<-c(12,13,14,15,16)

assignment<-c(assignment,sample(Years,1))

Moorings<-c("02","BS3","04","BS2","05","BS1","08")

assignment<-c(assignment,sample(Moorings,1))

if(species=="RW"){
  SFs<-100:200
}else if(species=="GS"){
  SFs<-50:100
}

assignment<-c(assignment,sample(SFs,1))

extent<-0:100

assignment<-c(assignment,sample(extent,1))

assTab<-rbind(assTab,t(data.frame(assignment)))

}

colnames(assTab)<-c("spec","year","Mooring","num_sf","per_through_dat_start")

write.csv(assTab,"E:/Other/",row.names=FALSE)
