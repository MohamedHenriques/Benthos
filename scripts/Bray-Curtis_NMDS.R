setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("vegan","ggplot2","viridis","RColorBrewer","psych","reshape2","beepr","data.table")
lapply(packs,require,character.only=T)


DB66<-fread("data_out/db/Final_DB_lowtaxa_density_polyexcl_20201208.csv") ### created in script called Database_cleanup_joining
str(DB66)

DB67<-fread("data_out/db/Final_DB_family_density_polyexcl_20201208.csv") ### created in script called Database_cleanup_joining
str(DB67)
names(DB67)[8]<-"numb"


##Remove non-target benthos
unique(DB66$low_taxa[which(DB66$class1=="Other")])

db<-DB66[-which(DB66$class1=="Other"),]
unique(db$low_taxa)


##Remove non-target benthos
unique(DB67$family[which(DB67$class1=="Other")])

dbf<-DB67[-which(DB67$class1=="Other"),]
unique(dbf$family)

###remove data from 2020
db1<-db[!year==2020]
db1[year==2020]

###remove data from 2020 fam
db1f<-dbf[!year==2020]
db1f[year==2020]

###aggregate and reshape database for analysis
db2<-db1[,lapply(.SD,sum,na.rm=T),.SDcols="numb",by=c("site","month","coreID","low_taxa")]
DB<-dcast.data.table(db2,coreID+site+month~low_taxa,value.var="numb")
setkey(DB,coreID,site,month)

###aggregate and reshape database for analysis
db2f<-db1f[,lapply(.SD,sum,na.rm=T),.SDcols="numb",by=c("site","month","coreID","family")]
DBf<-dcast.data.table(db2f,coreID+site+month~family,value.var="numb")
setkey(DBf,coreID,site,month)

###Calculate densities
dens1<-function(x){x/0.00866}
dens2<-function(x){x/0.00817}
DB1<-DB[,lapply(.SD,ifelse(site=="AD",dens1,dens2)),by=c("coreID","site","month")]
setkey(DB1,coreID,month,site)
DB1f<-DBf[,lapply(.SD,ifelse(site=="AD",dens1,dens2)),by=c("coreID","site","month")]


###Limit database to the 75% most abundant species
perc<-function(x){x/sum(x)}
test<-DB1[,!c("coreID","month")][,lapply(.SD,mean),by=site]
test1<-sapply(test,perc)
#nam<-names(test)[2:ncol(test)]
test[,lapply(.SD,perc),by=site]
rows
###NMDS with Bray curtis
##Prepare data and remove rows with zeros in all columns

data<-DB1[apply(DB1[,!(1:3)],1,sum)!=0]
dataf<-DB1f[apply(DB1f[,!(1:3)],1,sum)!=0]

data1<-data[,!(1:3)]
table(is.na(data1))
data1f<-dataf[,!(1:3)]
table(is.na(data1f))

data2<-data[,1:3]
head(data2)
data2f<-dataf[,1:3]
head(data2f)

set.seed(2)
NMDS<-metaMDS(data1,distance = "bray",k=2,trymax=100)
beep()

plot(NMDS)
stressplot(NMDS)

set.seed(31)
NMDSf<-metaMDS(data1f,distance = "bray",k=3,trymax=1000)
beep()

plot(NMDSf)
stressplot(NMDSf)




#extract NMDS scores (x and y coordinates)
data.scores = as.data.table(scores(NMDS))

#add columns to data frame 
data.scores$coreID = data2$coreID
data.scores$site = data2$site
data.scores$month = data2$month

head(data.scores)

ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(colour=factor(month),shape=site))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Month", y = "NMDS2",shape="Site")+
  scale_color_brewer(palette="Spectral")



ordiplot(NMDS,type="n")
orditorp(NMDS,display="species",col="red",air=0.01)
orditorp(NMDS,display="sites",cex=0.5,air=0.01)


ordiplot(NMDS,type="n")
ordihull(NMDS,groups=data2$site,draw="polygon",col="grey90",label=F)
orditorp(NMDS,display="species",col="red",air=0.01)
orditorp(NMDS,display="sites",col=c(rep("green",5),rep("blue",5)),
         air=0.01,cex=1)
