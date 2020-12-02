setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("gridExtra","vegan","ggplot2","viridis","effects","RColorBrewer","xlsx","psych","reshape2","tidyr")
lapply(packs,require,character.only=T)

db_mh<-read.xlsx("D:/Work/FCUL/Doutoramento/Capitulos/Exclosure_experiments/Data/Cores_ID/DBMH/Cores_DB_all_20201105.xlsx",1)
str(db_mh)

unique(db_mh$low_taxa)

db_mh1<-db_mh[which(db_mh$treatment=="C"),c(2:5,7:9,11:12,15,19:20,29)]
db_mh2<-db_mh1[-which(db_mh1$done=="N"),]
str(db_mh2)
db_mh2$site<-"AD"
db_mh2$coreID_ok<-paste(db_mh2$coreID,db_mh2$round,sep="")
unique(db_mh2$coreID_ok)

unique(db_mh2$class)

db_mh2$class1<-as.character(db_mh2$class)
db_mh2$class1[which(db_mh2$subclass=="Sedentaria")]<-"Polychaeta_sedentaria"
db_mh2$class1[which(db_mh2$subclass=="Errantia")]<-"Polychaeta_errantia"
#db_mh2$class1[which(db_mh2$subclass=="NID")]<-"Polychaeta_NID" ## ja nao existem casos destes, mudei tudo paa sedentaria na base de dados original
db_mh2$class1[which(db_mh2$class1=="Enteropneusta")]<-"Other"
db_mh2$class1[which(db_mh2$class1=="NID")]<-"Other"


unique(db_mh2$class1)
db_mh2[db_mh2$class1=="Enteropneusta",]

db_ac<-read.xlsx("D:/Work/FCUL/Doutoramento/Capitulos/Spatial_temporal_variability_invertebrates/Data_Ana/Cores_V22.xlsx",1)
str(db_ac)
db_ac1<-db_ac[,c(2,4,5:7,9,11,15:16)]
db_ac1$numb<-ifelse(db_ac1$final_name=="empty",0,1)
db_ac1$coreID<-substr(db_ac1$point,1,nchar(as.character(db_ac1$point))-1)
db_ac1$class1<-ifelse(db_ac1$class=="BIV","Bivalvia",
                      ifelse(db_ac1$class=="SED","Polychaeta_sedentaria",
                             ifelse(db_ac1$class=="ERR","Polychaeta_errantia",
                                    ifelse(db_ac1$class=="GAS","Gastropoda",
                                           ifelse(db_ac1$class=="CRU","Malacostraca",
                                                  ifelse(db_ac1$class=="OTH","Other","empty"))))))
length(unique(db_ac1$coreID))
unique(db_mh2$class1)


### Uniformizar
db_mh3<-db_mh2[,c(14,2:3,15,6,8,16,10:12)]
colnames(db_mh3)[4]<-"coreID"
db_ac2<-db_ac1[,c(3,1,2,11,5,7,12,9,8,10)]
names(db_ac2)<-names(db_mh3)

### aggregate to have one line per taxa per core per all other stuff
db_mh4<-aggregate(db_mh3$numb,by=list(site=db_mh3$site,year=db_mh3$year,month=db_mh3$month,coreID=db_mh3$coreID,depth=db_mh3$depth,
                               class1=db_mh3$class1,family=db_mh3$family,low_taxa=db_mh3$low_taxa),FUN=sum)
str(db_mh4)
#db_mh4$coreID<-as.character(db_mh4$coreID)
db_mh4$depth<-as.character(db_mh4$depth)
db_mh4$family<-as.character(db_mh4$family)
db_mh4$low_taxa<-as.character(db_mh4$low_taxa)
unique(db_mh4$low_taxa)



db_ac3<-aggregate(db_ac2$numb,by=list(site=db_ac2$site,year=db_ac2$year,month=db_ac2$month,coreID=db_ac2$coreID,depth=db_ac2$depth,
                                      class1=db_ac2$class1,family=db_ac2$family,low_taxa=db_ac2$low_taxa),FUN=sum)
unique(db_ac3$low_taxa)

names(db_mh3)
names(db_ac2)

DB<-rbind(db_mh3,db_ac2)
str(DB)
unique(DB$site)
unique(DB$year)
unique(DB$month)
length(unique(DB$coreID))
unique(DB$depth)
unique(as.character(DB$class))
unique(DB$class1)
unique(DB$family)


DB1<-aggregate(DB$numb,by=list(site=DB$site,year=DB$year,month=DB$month,coreID=DB$coreID,depth=DB$depth,
                               class1=DB$class1,family=DB$family,low_taxa=DB$low_taxa),FUN=sum)
str(DB1)
colnames(DB1)[9]<-"numb"
unique(as.character(DB1$low_taxa))
length(unique(DB1$coreID))

#write.table(DB1,"data_out/db/FinalDB_aggregated.csv",sep=";",row.names = F)

### Calculate density of individuals

##### Complete database with zeros

DB12<-DB1[,c("low_taxa","year","month","coreID","depth","numb","site")]

DB13<-DB1[,c("low_taxa","class1","family")]
DB131<-aggregate(DB13$class1,by=list(low_taxa=DB13$low_taxa,class1=DB13$class1),FUN=length)
DB131<-DB131[,-3]

DB132<-aggregate(DB13$family,by=list(low_taxa=DB13$low_taxa,family=DB13$family),FUN=length)
DB132<-DB132[,-3]


DB12$low_taxa<-as.character(DB12$low_taxa)
length(unique(DB12$coreID))

DB2<-complete(DB12,low_taxa,nesting(year,month,site,coreID),fill = list(numb=0))
#write.table(DB2,"data_out/db/Total_complete_db.csv",row.names=F,sep=";")




length(unique(DB2$coreID))
length(unique(DB2$coreID[DB2$low_taxa=="Afruca_tangeri"]))
length(unique(DB2$coreID[DB2$low_taxa=="Austromacoma_nymphalis"]))
str(DB2)
unique(DB2$low_taxa[DB2$coreID=="C14T2"])
unique(DB2$low_taxa[DB2$coreID=="C7A7"])

unique(DB1$low_taxa)
unique(DB2$low_taxa)
  
DB3<-merge(DB2,DB132,by="low_taxa", all.x=T)
DB4<-merge(DB3,DB131,by="low_taxa", all.x=T)
DB5<-DB4[-which(DB4$low_taxa=="empty"),] # too remove taxa lines of empty for each core, but also repeated for T and B 
#write.table(DB5,"data_out/db/Final_DB_complete_20201123.csv",row.names=F,sep=";")

DB6<-aggregate(DB5$numb,by=list(year=DB5$year,month=DB5$month,site=DB5$site,coreID=DB5$coreID,class1=DB5$class1,
                                family=DB5$family,low_taxa=DB5$low_taxa),FUN=sum)
colnames(DB6)[8]<-"numb"


str(DB6)
length(unique(DB6$low_taxa))
length(unique(DB6$coreID))
length(unique(DB6$year))
length(unique(DB6$month))
length(unique(DB6$site))
length(unique(DB6$family))
length(unique(DB6$class1))
table(is.na(DB6))

#### Calcular densidades

DB6$dens<-ifelse(DB6$site=="AD",DB6$numb/0.0113,DB6$numb/0.00785)

#write.table(DB6,"data_out/db/DB6.csv",row.names=F,sep=";")

### include column with cut for most abundant low taxa
DB7<-aggregate(DB6$dens,by=list(family=DB6$family),FUN=mean)
DB7$cut<-ifelse(DB7$x>=10,"Y","N")
DB8<-merge(DB6,DB7,by="family",all.x=T)

##Order site by island
DB8$site1<-factor(DB8$site,levels=c("AD","BI","BR","E","A","AB"))


####plots de densidade total (todos os sitios todos os meses)

ggplot(DB8[DB8$cut=="Y",],aes(x=dens,y=reorder(low_taxa, dens),col=class1)) +
  #geom_point(size=3)+
  stat_summary(size=1)+
  theme_bw() +
  labs(x="Mean density (ind.m-2)",y="Lower tax.level (taxa above 10 ind.m-2)")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey60", linetype="dashed"),
        axis.title = element_text(size=16,face="bold"))+
  scale_colour_brewer(palette="Dark2")

ggplot(DB8[DB8$cut=="Y",],aes(x=dens,y=reorder(family, dens),col=class1)) +
  #geom_point(size=3)+
  stat_summary(size=1)+
  theme_bw() +
  labs(x="Mean density (ind.m-2)",y="Family (families above 10 ind.m-2)")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey60", linetype="dashed"),
        axis.title = element_text(size=16,face="bold"))+
  scale_colour_brewer(palette="Dark2")

ggplot(DB8,aes(x=dens,y=reorder(class1, dens))) +
  #geom_point(size=3)+
  stat_summary(size=1)+
  theme_bw() +
  labs(x="Mean density (ind.m-2)",y="(Sub)class")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey60", linetype="dashed"),
        axis.title = element_text(size=16,face="bold"))

#### Plots por sítio

###Barras percentuais
ggplot(DB6, aes(x=dens,y=reorder(low_taxa, dens),fill=site)) +
  geom_bar(stat="identity",position="fill")+
  scale_fill_brewer(palette="Set1")+
  theme_bw()+
  labs(x="Density (ind.m-2)",y="Lower taxa (ordered by total density)")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=8),
        axis.title = element_text(size=16))
#stat_summary()


ggplot(DB6, aes(x=dens,y=reorder(family, dens),fill=site)) +
  geom_bar(stat="identity",position="fill")+
  scale_fill_brewer(palette="Set1")+
  theme_bw()+
  labs(x="Density (ind.m-2)",y="Family (ordered by total density)")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=16))
  #stat_summary()

ggplot(DB6, aes(x=dens,y=reorder(class1, dens),fill=site)) +
  geom_bar(stat="identity",position="fill")+
  scale_fill_brewer(palette="Set1")+
  theme_bw()+
  labs(x="Density (ind.m-2)",y="Class (ordered by total density)")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=16))
#stat_summary()

### Barras normais

ggplot(DB8[DB8$cut=="Y",], aes(x=dens,y=reorder(low_taxa, dens),fill=site1)) +
  stat_summary(geom="bar",size=1)+
  scale_fill_manual(values=c("red","steelblue2","royalblue3","darkblue","limegreen","darkgreen"))+
  theme_bw() +
  labs(x="Mean density (ind.m-2)",y="Lower tax.level (families above 10 ind.m-2)")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey60", linetype="dashed"),
        axis.title = element_text(size=16,face="bold"))

XX<-aggregate(DB8$dens,by=list(family=DB8$family,site1=DB8$site1),FUN=mean)
colnames(XX)[3]<-"dens"
XX1<-merge(XX,DB7,by="family", all.x=T)

XXX<-dcast(data=XX,family~site1)

ggplot(XX1[XX1$cut=="Y",],aes(x=dens,y=reorder(family, dens),fill=site1)) +
  #stat_summary(geom="bar",size=1)+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("red","steelblue2","royalblue3","darkblue","limegreen","darkgreen"))+
  theme_bw() +
  labs(x="Mean density (ind.m-2)",y="Family (families above 10 ind.m-2)")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey60", linetype="dashed"),
        axis.title = element_text(size=16,face="bold"))


ggplot(DB8[DB8$cut=="Y",], aes(x=dens,y=reorder(family, dens),col=site1)) +
  stat_summary(geom="pointrange",size=1)+
  scale_fill_manual(values=c("red","steelblue2","royalblue3","darkblue","limegreen","darkgreen"))+
  theme_bw() +
  labs(x="Mean density (ind.m-2)",y="Family (families above 10 ind.m-2)")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey60", linetype="dashed"),
        axis.title = element_text(size=16,face="bold"))

ggplot(DB8,aes(y=dens,x=reorder(class1, -dens),fill=site1)) +
  stat_summary(geom="bar",size=1,position="dodge")+
  scale_fill_manual(values=c("red","steelblue2","royalblue3","darkblue","limegreen","darkgreen"))+
  theme_bw() +
  labs(x="Mean density (ind.m-2)",y="(Sub)class")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey60", linetype="dashed"),
        axis.title = element_text(size=16,face="bold"))


### 


ggplot(DB6, aes(x=dens,y=reorder(family, dens))) +
  #geom_point(size=3)+
  stat_summary(size=1)+
  facet_grid(. ~ site,scales="free")+
  theme_bw() +
  labs(x="Density (ind.m-2)",y="Family (ordered by total density)")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey60", linetype="dashed"))

ggplot(DB6, aes(x=dens,y=reorder(class1, dens)))+
  stat_summary(size=1)+
  facet_grid(. ~ site,scales="free")+
  theme_bw() +
  labs(x="Density (ind.m-2)",y="Class (ordered by total density)")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey60", linetype="dashed"))

#### Temporal



DB8$time<-paste(DB8$year,DB8$month,sep="")
DB8$time1<-factor(DB8$time, levels=c("20181","20182","20183","201810","201811","201812","20191","20193","20194","201910","20201","20202","20203"))
DB8$time2<-ifelse(DB8$time1=="20181","01/01/2018",ifelse(DB8$time1=="20182","01/02/2018",ifelse(DB8$time1=="20183","01/03/2018",
                  ifelse(DB8$time1=="201810","01/10/2018",ifelse(DB8$time1=="201811","01/11/2018",ifelse(DB8$time1=="201812","01/12/2018",
                  ifelse(DB8$time1=="20191","01/01/2019",ifelse(DB8$time1=="20193","01/03/2019",ifelse(DB8$time1=="20194","01/04/2019",
                  ifelse(DB8$time1=="201910","01/10/2019",ifelse(DB8$time1=="20201","01/01/2020",ifelse(DB8$time1=="20202","01/02/2020",
                  ifelse(DB8$time1=="20203","01/03/2020",NA)))))))))))))

DB8$time3<-as.POSIXct(DB8$time2,format="%d/%m/%Y",tz="GMT")

DB8$month1<-factor(DB8$month, levels=c("10","11","12","1","2","3","4"))

class(DB8$time3)
unique(DB8$time3)

DATES<-c("JAN18","FEB18","MAR18","OCT18","NOV18","DEC18","JAN19","MAR19","APR19","OCT19","JAN20","FEB20","MAR20")
names(DATES)<-levels(factor(unique(DB8$time3)))
SITES<-c("Adonga","Bijante","Bruce","Escadinhas","Anrumai","Abu")
names(SITES)<-levels(factor(unique(DB8$site1)))
#coll<-c("red","steelblue2","royalblue3","darkblue","limegreen","darkgreen")

#All dates
ggplot(DB8,aes(x=dens,y=reorder(class1, dens),col=class1)) +
  #geom_point(size=3)+
  stat_summary(size=1)+
  facet_grid(time3 ~ site1,scales="free",labeller = labeller(time3=DATES,site1=SITES))+
  theme_bw() +
  labs(x="Density (ind.m-2)",y="Class (ordered by total density)")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey60", linetype="dashed"))+
  scale_colour_brewer(palette="Dark2")


#Merging months

DATES1<-c("OCT","NOV","DEZ","JAN","FEB","MAR","APR")
names(DATES1)<-levels(factor(unique(DB8$month1)))

ggplot(DB8, aes(x=dens,y=reorder(class1, dens),col=class1)) +
  #geom_point(size=3)+
  stat_summary(size=1)+
  facet_grid(month1 ~ site1,scales="free",labeller = labeller(month1=DATES1,site1=SITES))+
  theme_bw() +
  labs(x="Density (ind.m-2)",y="Class (ordered by total density)")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey60", linetype="dashed"))+
  scale_colour_brewer(palette="Dark2")



#### More detailed temporal graphs
### per site vs class


ggplot(DB8,aes(x=month1,y=dens,col=class1,group=class1))+
  stat_summary(geom="pointrange",position=position_dodge(width=0.3),size=0.7)+
  stat_summary(geom="line",lwd=1.5,alpha=0.50)+
  facet_wrap(.~site1,scales="free",labeller = labeller(site1=SITES))+
  theme_bw()+
  labs(x="Month (aggregation of all years)",y="Mean density (ind.m-2)")+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=12),axis.title = element_text(size=15,face="bold"))+
  theme(strip.text = element_text(face="bold", size=rel(1.5)))




ggplot(DB6,aes(x=time3,y=dens,col=site,group=site))+
  stat_summary(geom="pointrange",position=position_dodge(width=0.6),size=0.7)+
  stat_summary(geom="line")+
  facet_wrap(.~class1,scales="free")+
  theme_bw()+
  labs(x="Time (months of sampling)",y="mean of density (ind.m-2)")+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  scale_x_datetime(breaks=c(unique(DB6$time3)))


###Aggregates to check values of mean densities per site 

DB77<-aggregate(DB6$dens,by=list(low_taxa=DB6$low_taxa,site=DB6$site),FUN=mean)
DB777<-dcast(data=DB77,low_taxa~site)

DB7777<-aggregate(DB6$dens,by=list(low_taxa=DB6$low_taxa),FUN=mean)

DB_exp<-merge(DB777,DB7777,by="low_taxa")

#write.table(DB_exp,"data_out/db/DB_exp.csv",sep=";",row.names=F)





ggplot(DB8[which(DB8$cut=="Y"),],aes(x=month1,y=dens,col=site,group=site))+
  stat_summary(geom="pointrange",position=position_dodge(width=0.6),size=0.7)+
  stat_summary(geom="line")+
  facet_wrap(.~family,scales="free")+
  theme_bw()+
  labs(x="month (aggregation of all years",y="mean of density (ind.m-2)")

ggplot(DB8[which(DB8$cut=="Y"),],aes(x=time3,y=dens,col=site,group=site))+
  stat_summary(geom="pointrange",position=position_dodge(width=0.6),size=0.7)+
  stat_summary(geom="line")+
  facet_wrap(.~family,scales="free")+
  theme_bw()+
  labs(x="Time (months of sampling)",y="mean of density (ind.m-2)")+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  scale_x_datetime(breaks=c(unique(DB6$time3)),labels=DATES)



####






