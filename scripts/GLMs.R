setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("MASS","ggplot2","viridis","RColorBrewer","psych","reshape2","beepr","data.table","car")
lapply(packs,require,character.only=T)


DB68<-fread("data_out/db/Final_DB_class1_density_polyexcl_20201208.csv") ### created in script called Database_cleanup_joining
str(DB68)

time<-c("10","11","12","1","2","3","4")
DB68[,month1:=factor(month,levels=time)]
DB68[,month2:=ifelse(month==10,1,ifelse(month==11,2,ifelse(month==12,3,ifelse(month==1,4,ifelse(month==2,5,ifelse(month==3,6,7))))))]
DB68[,unique(month2)]

DB68[class1=="Other",class1]

DB<-DB68[!class1=="Other"]

##Take a look at data
###site
ggplot(DB,aes(x=class1,y=dens,colour=site))+
  stat_summary(geom = "pointrange",position=position_dodge(width=0.4),size=1)+
  theme_bw()

###month as factor
ggplot(DB,aes(x=month1,y=dens,colour=site))+
  stat_summary(geom = "pointrange",position=position_dodge(width=0.4),size=1)+
  theme_bw()

###month as numeric
ggplot(DB,aes(x=month2,y=dens,colour=site))+
  stat_summary(geom = "pointrange",position=position_dodge(width=0.4),size=1)+
  scale_x_continuous(breaks=1:7,labels=c(10,11,12,1,2,3,4))+
  stat_smooth(method="loess",se=F)+
  theme_bw()

D<-aggregate(DB$dens,by=list(site=DB$site,class1=DB$class1),FUN=mean)
D1<-dcast(class1~site,data=D)


#### GLM

###### Include area of cores in DB

DB[,areacore:=ifelse(site=="AD",0.0113,0.00785)]
DB[,areacore]

#######Overall

m1<-glm.nb(x~site*month2+offset(log(areacore)),data=DB)
summary(m1)
anova(m1)

m1a<-update(m1,~.-site:month2)
drop1(m1a, test="Chi")
anova(m1,m1a)
AIC(m1,m1a)
summary(m1a)

######## Bivalves
dbb<-DB[class1=="Bivalvia"]

m2<-glm.nb(x~site*month2+offset(log(areacore)),data=dbb)
summary(m2)
anova(m2)
m2a<-update(m2,~.-site:month2)
drop1(m2a,test="Chi")
anova(m2,m2a)
AIC(m2,m2a)
summary(m2a)

residualPlots(m2a)

est<-cbind(coef(m2a),confint(m2a))
exp(est)

######## errant polychaetes
dbpe<-DB68[class1=="Polychaeta_errantia"]

m3<-glm.nb(x~site*month2+offset(log(areacore)),data=dbpe)
summary(m3)
anova(m3)
m3a<-update(m3,~.-site:month2)
drop1(m3a,test="Chi")
anova(m3,m3a)
AIC(m3,m3a)
summary(m3)

residualPlots(m3)

### etc for the others



