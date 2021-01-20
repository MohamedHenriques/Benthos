setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("gridExtra","vegan","ggplot2","viridis","effects","RColorBrewer","xlsx","psych","reshape2","tidyr")
lapply(packs,require,character.only=T)

DB66<-read.table("data_out/db/Final_DB_lowtaxa_density_polyexcl_20201208.csv",header=T,sep=";") ### created in script called Database_cleanup_joining
DB67<-read.table("data_out/db/Final_DB_family_density_polyexcl_20201208.csv",header=T,sep=";") ### created in script called Database_cleanup_joining
DB68<-read.table("data_out/db/Final_DB_class1_density_polyexcl_20201208.csv",header=T,sep=";") ### created in script called Database_cleanup_joining


ggplot(DB68,aes(x=class1,y=dens,colour=site))+
  stat_summary(geom = "pointrange",position=position_dodge(width=0.4),size=1)+
  theme_bw()
  
ggplot(DB68,aes(x=class1,y=dens1,colour=site))+
  stat_summary(geom = "pointrange",position=position_dodge(width=0.4),size=1)+
  theme_bw()

D<-aggregate(DB68$dens,by=list(site=DB68$site,class1=DB68$class1),FUN=mean)
D1<-dcast(class1~site,data=D)


### include column with cut for most abundant families
DB7<-aggregate(DB67$dens,by=list(family=DB67$family),FUN=mean)
colnames(DB7)[2]<-"dens"

DB7$cut<-ifelse(DB7$dens>=10,"Y","N")

colnames(DB67)[8]<-"numb"
DB8<-merge(DB67,DB7,by="family",all.x=T)
colnames(DB8)
DB8$dens.y<-NULL
colnames(DB8)[9]<-"dens"

### include column with cut for most abundant low_taxa
DB77<-aggregate(DB66$dens,by=list(low_taxa=DB66$low_taxa),FUN=mean)
colnames(DB77)[2]<-"dens"

DB77$cut<-ifelse(DB77$dens>=10,"Y","N")

DB88<-merge(DB66,DB77,by="low_taxa",all.x=T)
colnames(DB88)
DB88$dens.y<-NULL
colnames(DB88)[10]<-"dens"

##Order site by island
DB88$site1<-factor(DB88$site,levels=c("AD","BI","BR","E","A","AB"))


####plots de densidade total (todos os sitios todos os meses)

ggplot(DB88[DB88$cut=="Y",],aes(x=dens,y=reorder(low_taxa, dens),col=class1)) +
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

ggplot(DB68,aes(x=dens,y=reorder(class1, dens))) +
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
ggplot(DB66, aes(x=dens,y=reorder(low_taxa, dens),fill=site)) +
  geom_bar(stat="identity",position="fill")+
  scale_fill_brewer(palette="Set1")+
  theme_bw()+
  labs(x="Density (ind.m-2)",y="Lower taxa (ordered by total density)")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=8),
        axis.title = element_text(size=16))
#stat_summary()


ggplot(DB67, aes(x=dens,y=reorder(family, dens),fill=site)) +
  geom_bar(stat="identity",position="fill")+
  scale_fill_brewer(palette="Set1")+
  theme_bw()+
  labs(x="Density (ind.m-2)",y="Family (ordered by total density)")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=16))
#stat_summary()

ggplot(DB68, aes(x=dens,y=reorder(class1, dens),fill=site)) +
  geom_bar(stat="identity",position="fill")+
  scale_fill_brewer(palette="Set1")+
  theme_bw()+
  labs(x="Density (ind.m-2)",y="Class (ordered by total density)")+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=16))
#stat_summary()

### Barras normais

#ggplot(DB88[DB88$cut=="Y",], aes(x=dens,y=reorder(low_taxa, dens),fill=site1)) +
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

XX<-aggregate(DB67$dens,by=list(family=DB67$family,site1=DB67$site),FUN=mean)
colnames(XX)[3]<-"dens"
XX1<-merge(XX,DB7,by="family", all.x=T)

XXX<-dcast(data=XX,family~site1)

YY<-aggregate(DB66$dens,by=list(low_taxa=DB66$low_taxa,site1=DB66$site),FUN=mean)
colnames(YY)[3]<-"dens"

#DB77<-aggregate(DB61$dens,by=list(low_taxa=DB61$low_taxa),FUN=mean)
YY1<-merge(YY,DB66,by="low_taxa", all.x=T)

YYY<-dcast(data=YY,low_taxa~site1)

ggplot(XX1[XX1$cut=="Y",],aes(x=dens.x,y=reorder(family, dens.x),fill=site1)) +
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

ggplot(DB88,aes(y=dens,x=reorder(class1, -dens),fill=site1)) +
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


ggplot(DB8, aes(x=dens,y=reorder(family, dens))) +
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

ggplot(DB88, aes(x=dens,y=reorder(class1, dens)))+
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



DB68$time<-paste(DB68$year,DB68$month,sep="")
DB68$time1<-factor(DB68$time, levels=c("20181","20182","20183","201810","201811","201812","20191","20192","20193","20194","201910","20201","20202","20203"))
DB68$time2<-ifelse(DB68$time1=="20181","01/01/2018",ifelse(DB68$time1=="20182","01/02/2018",ifelse(DB68$time1=="20183","01/03/2018",
                ifelse(DB68$time1=="201810","01/10/2018",ifelse(DB68$time1=="201811","01/11/2018",ifelse(DB68$time1=="201812","01/12/2018",
                ifelse(DB68$time1=="20191","01/01/2019",ifelse(DB68$time1=="20192","01/02/2019",ifelse(DB68$time1=="20193","01/03/2019",ifelse(DB68$time1=="20194","01/04/2019",
                ifelse(DB68$time1=="201910","01/10/2019",ifelse(DB68$time1=="20201","01/01/2020",ifelse(DB68$time1=="20202","01/02/2020",
                ifelse(DB68$time1=="20203","01/03/2020",NA))))))))))))))

DB68$time3<-as.POSIXct(DB68$time2,format="%d/%m/%Y",tz="GMT")

DB68$month1<-factor(DB68$month, levels=c("10","11","12","1","2","3","4"))


DB68$site1<-factor(DB68$site,levels=c("AD","BI","BR","E","A","AB"))
class(DB68$time3)
unique(DB68$time3)


DATES<-c("JAN18","FEB18","MAR18","OCT18","NOV18","DEC18","JAN19","FEB19","MAR19","APR19","OCT19","JAN20","FEB20","MAR20")
names(DATES)<-levels(factor(unique(DB68$time1)))
SITES<-c("Adonga","Bijante","Bruce","Escadinhas","Anrumai","Abu")
names(SITES)<-levels(factor(unique(DB68$site1)))
#coll<-c("red","steelblue2","royalblue3","darkblue","limegreen","darkgreen")

#All dates
ggplot(DB68,aes(x=dens,y=reorder(class1, dens),col=class1)) +
  #geom_point(size=3)+
  stat_summary(size=1)+
  facet_grid(time1 ~ site1,scales="free",labeller = labeller(time1=DATES,site1=SITES))+
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
names(DATES1)<-levels(factor(unique(DB68$month1)))

ggplot(DB68, aes(x=dens,y=reorder(class1, dens),col=class1)) +
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


ggplot(DB68,aes(x=month1,y=dens,col=class1,group=class1))+
  stat_summary(geom="pointrange",position=position_dodge(width=0.3),size=0.7)+
  stat_summary(geom="line",lwd=1.5,alpha=0.50)+
  facet_wrap(.~site1,labeller = labeller(site1=SITES))+
  theme_bw()+
  labs(x="Month (aggregation of all years)",y="Mean density (ind.m-2)")+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=12),axis.title = element_text(size=15,face="bold"))+
  theme(strip.text = element_text(face="bold", size=rel(1.5)))


DB88$month1<-factor(DB88$month, levels=c("10","11","12","1","2","3","4"))

ggplot(DB88[DB88$cut=="Y",],aes(x=month1,y=dens,col=low_taxa,group=low_taxa))+
  #stat_summary(geom="pointrange",position=position_dodge(width=0.6),size=0.7)+
  geom_smooth(geom="smooth",method="auto",se=F,size=1.5)+
  facet_wrap(.~site1,scales="free")+
  theme_bw()+
  labs(x="Time (months of sampling)",y="mean of density (ind.m-2)")+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
  #scale_colour_manual(values=c("red","steelblue2","royalblue3","darkblue","limegreen","darkgreen"))
  #scale_x_datetime(breaks=c(unique(DB68$time3)))


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
