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


db2<-aggregate(db1$numb,by=list(year=db1$year,month=db1$month,yearmonth=db1$yearmonth,coreID=db1$coreID,depth=db1$depth,empty=db1$empty,class=db1$class,family=db1$family,low_taxa=db1$low_taxa),FUN=sum)
colnames(db2)[10]<-"numb"

levels(db2$low_taxa)

db3<-aggregate(db1$numb,by=list(family=db1$family,low_taxa=db1$low_taxa),FUN=sum)
colnames(db3)[3]<-"numb"

#write.table(db3,"./data_out/db/Total_sum_inv.csv",sep=";",row.names=F)

a<-describeBy(db2$numb,group=list(db2$low_taxa,db2$yearmonth))
