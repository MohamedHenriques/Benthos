setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("ggplot2","viridis","effects","RColorBrewer","xlsx","psych")
lapply(packs,require,character.only=T)

##LOad database
db<-read.xlsx("D:/Work/FCUL/Doutoramento/Capitulos/Exclosure_experiments/Data/Cores_ID/DBMH/Cores_DB_all_20201105.xlsx",1)
#str(db)

###subset only for Control cores with existing tubes
db1<-db[which(db$treatment=="C"&db$exists=="Y"&db$SIA!="Y"&db$done=="Y"),]

###subset for ll cores in control treatment that have been already identified
DB1<-db[which(db$treatment=="C"&db$done=="Y"),]

### Create base db for cores with existing animals
db2<-aggregate(db1$numb,by=list(year=db1$year,month=db1$month,yearmonth=db1$yearmonth,coreID=db1$coreID,depth=db1$depth,empty=db1$empty,class=db1$class,family=db1$family,low_taxa=db1$low_taxa),FUN=sum)
colnames(db2)[10]<-"numb"

###Create base db for all identified cores
DB2<-aggregate(DB1$numb,by=list(year=DB1$year,month=DB1$month,yearmonth=DB1$yearmonth,coreID=DB1$coreID,depth=DB1$depth,empty=DB1$empty,class=DB1$class,family=DB1$family,low_taxa=DB1$low_taxa),FUN=sum)
colnames(DB2)[10]<-"numb"

levels(db2$low_taxa)
levels(DB2$low_taxa)

## Add a unique identifier for each core
DB2$code<-paste(DB2$coreID,DB2$yearmonth,sep="_")

### Freq absoluta for tubes with existing animals
db3<-aggregate(db1$numb,by=list(family=db1$family,low_taxa=db1$low_taxa),FUN=sum)
colnames(db3)[3]<-"numb"

#write.table(db3,"./data_out/db/Total_sum_inv.csv",sep=";",row.names=F)

###Freq absoluta para todo os cores
DB3<-aggregate(DB1$numb,by=list(family=DB1$family,low_taxa=DB1$low_taxa),FUN=sum)
colnames(DB3)[3]<-"FA"

##Freq numérica para todos os cores
DB3$FN<-round(DB3$FA/sum(DB3$FA)*100,digits=1)

#write.table(DB3,"./data_out/DB/Descriptive_MH.csv",sep=";",row.names=F)

###Calculate sample sizes for each month
Ncore<-length(unique(DB2$code))
Ncore_1810<-length(unique(DB2$code[which(DB2$yearmonth==201810)]))
Ncore_1903<-length(unique(DB2$code[which(DB2$yearmonth==201903)]))
Ncore_1910<-length(unique(DB2$code[which(DB2$yearmonth==201910)]))


