setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("vegan","ggplot2","viridis","RColorBrewer","psych","reshape2","beepr","data.table")
lapply(packs,require,character.only=T)


DB66<-fread("data_out/db/Final_DB_lowtaxa_density_polyexcl_20201208.csv") ### created in script called Database_cleanup_joining
str(DB66)

##Remove non-target benthos
unique(DB66$low_taxa[which(DB66$class1=="Other")])

db<-DB66[-which(DB66$class1=="Other"),]
unique(db$low_taxa)

###remove data from 2020
db1<-db[!year==2020]
db1[year==2020]

###aggregate and reshape database for analysis
db2<-db1[,lapply(.SD,sum,na.rm=T),.SDcols="numb",by=c("site","coreID","low_taxa")]
DB<-dcast.data.table(db2,coreID+site~low_taxa,value.var="numb")
setkey(DB,coreID,site)

###Calculate densities
dens1<-function(x){x/0.00866}
dens2<-function(x){x/0.00817}
DB1<-DB[,lapply(.SD,ifelse(site=="AD",dens1,dens2)),by=c("coreID","site")]

###NMDS with Bray curtis
##Prepare data and remove rows with zeros in all columns

data<-DB1[apply(DB1[,!(1:2)],1,sum)!=0]

data1<-data[,!(1:2)]
table(is.na(data1))

data2<-data[,1:2]

set.seed(1)
NMDS<-metaMDS(data1,distance = "bray",k=3,trymax=100)
beep()

plot(NMDS)
