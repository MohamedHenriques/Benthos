setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("ggplot2","viridis","effects","RColorBrewer","xlsx","psych")
lapply(packs,require,character.only=T)

db_mh<-read.xlsx("D:/Work/FCUL/Doutoramento/Capitulos/Exclosure_experiments/Data/Cores_ID/DBMH/Cores_DB_all_20201105.xlsx",1)
str(db_mh)
db_mh1<-db_mh[which(db_mh$treatment=="C"),c(2:4,7:9,11:12,15,19:20,29)]
db_mh2<-db_mh1[-which(db_mh1$done=="N"),]
str(db_mh2)
db_mh2$site<-"AD"

unique(db_mh2$class1)

db_mh2$class1<-as.character(db_mh2$class)
db_mh2$class1[which(db_mh2$subclass=="Sedentaria")]<-"Polychaeta_sedentaria"
db_mh2$class1[which(db_mh2$subclass=="Errantia")]<-"Polychaeta_errantia"
db_mh2$class1[which(db_mh2$subclass=="NID")]<-"Polychaeta_NID"
db_mh2$class1[which(db_mh2$class=="Enteropneusta")]<-"Other"


unique(db_mh2$class1)
db_mh2[db_mh2$class=="Enteropneusta",]

db_ac<-read.xlsx("D:/Work/FCUL/Doutoramento/Capitulos/Spatial_temporal_variability_invertebrates/Data_Ana/Cores_V22.xlsx",1)
str(db_ac)
db_ac1<-db_ac[,c(2,4,5:7,11,15:16)]
db_ac1$numb<-1
db_ac1$coreID<-paste(db_ac1$site,db_ac1$month,db_ac1$core,sep="")
db_ac1$class1<-ifelse(db_ac1$class=="BIV","Bivalvia",
                      ifelse(db_ac1$class=="SED","Polychaeta_sedentaria",
                             ifelse(db_ac1$class=="ERR","Polychaeta_errantia",
                                    ifelse(db_ac1$class=="GAS","Gastropoda",
                                           ifelse(db_ac1$class=="CRU","Malacostraca",
                                                  ifelse(db_ac1$class=="OTH","Other","empty"))))))
unique(db_ac1$class1)
unique(db_mh2$class1)


### Uniformizar
db_mh3<-db_mh2[,c(13,2:5,7,14,9:11)]
db_ac2<-db_ac1[,c(3,1,2,10,5,6,11,8,7,9)]
names(db_ac2)<-names(db_mh3)

DB<-rbind(db_mh3,db_ac2)



colnames(db_ac1)[6]<-"coreID"
