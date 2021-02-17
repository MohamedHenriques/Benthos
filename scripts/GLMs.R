setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("MASS","ggplot2","viridis","RColorBrewer","psych","reshape2","beepr","data.table","car")
lapply(packs,require,character.only=T)


DB69<-fread("data_out/db/Final_DB_total_otherexc_density_polyexcl_20210202.csv") ### created in script called Database_cleanup_joining
str(DB69)

time<-c("10","11","12","1","2","3","4")
DB69[,month1:=factor(month,levels=time)]
DB69[,month2:=ifelse(month==10,1,ifelse(month==11,2,ifelse(month==12,3,ifelse(month==1,4,ifelse(month==2,5,ifelse(month==3,6,7))))))]
DB69[,unique(month2)]

DB<-DB69

DB68<-fread("data_out/db/Final_DB_class1_density_polyexcl_20210202.csv")

DB68[,month1:=factor(month,levels=time)]
DB68[,month2:=ifelse(month==10,1,ifelse(month==11,2,ifelse(month==12,3,ifelse(month==1,4,ifelse(month==2,5,ifelse(month==3,6,7))))))]
DB68[,unique(month2)]

DB1<-DB68
str(DB)

sum(DB$dens)

tot<-mean(DB1[class1=="Bivalvia",dens])+
mean(DB1[class1=="Polychaeta_sedentaria",dens])+
mean(DB1[class1=="Polychaeta_errantia",dens])+
mean(DB1[class1=="Malacostraca",dens])+
mean(DB1[class1=="Gastropoda",dens])+
mean(DB1[class1=="Other",dens])

DB1[class1=="Bivalvia",dens==0]

b<-mean(DB1[class1=="Bivalvia",dens])/tot*100
ps<-mean(DB1[class1=="Polychaeta_sedentaria",dens])/tot*100
pe<-mean(DB1[class1=="Polychaeta_errantia",dens])/tot*100
m<-mean(DB1[class1=="Malacostraca",dens])/tot*100
g<-mean(DB1[class1=="Gastropoda",dens])/tot*100
o<-mean(DB1[class1=="Other",dens])/tot*100
sum(b+ps+pe+m+g+o)

tot1<-sum(DB1[class1=="Bivalvia",dens])+
  sum(DB1[class1=="Polychaeta_sedentaria",dens])+
  sum(DB1[class1=="Polychaeta_errantia",dens])+
  sum(DB1[class1=="Malacostraca",dens])+
  sum(DB1[class1=="Gastropoda",dens])+
  sum(DB1[class1=="Other",dens])

b1<-sum(DB1[class1=="Bivalvia",dens])/tot1*100
ps1<-sum(DB1[class1=="Polychaeta_sedentaria",dens])/tot1*100
pe1<-sum(DB1[class1=="Polychaeta_errantia",dens])/tot1*100
m1<-sum(DB1[class1=="Malacostraca",dens])/tot1*100
g1<-sum(DB1[class1=="Gastropoda",dens])/tot1*100
o1<-sum(DB1[class1=="Other",dens])/tot1*100
sum(b1+ps1+pe1+m1+g1+o1)


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


#### GLM

DB[site=="E",mean(dens)]

###### Include area of cores in DB

DB[,areacore:=ifelse(site=="AD",0.0113,0.00785)]
DB[,areacore]

###### Include area of cores in DB1

DB1[,areacore:=ifelse(site=="AD",0.0113,0.00785)]
DB1[,areacore]


#######Overall

DB$site1<-factor(DB$site,levels=c("AD","A","AB","BI","BR","E"))

m1<-glm.nb(x~site1*month2+offset(log(areacore)),data=DB)
summary(m1)
anova(m1,test="Chisq")
anova.glm(m1)
unique(DB$month2)



m1a<-update(m1,~.-site:month2)
drop1(m1a, test="Chi")
anova(m1,m1a)
AIC(m1,m1a)
summary(m1a)

######## Bivalves
dbb<-DB1[class1=="Bivalvia"]

m2<-glm.nb(x~site*month2+offset(log(areacore)),data=dbb)
summary(m2)
anova(m2)
m2a<-update(m2,~.-site:month2)
drop1(m2a,test="Chi")
anova(m2,m2a)
AIC(m2,m2a)
summary(m2a)

#####Checking residuals
residualPlots(m2a)

E2<-resid(m2a,type='pearson')
F2<-predict(m2a,type="link")
gamma.dispersion(m2a,modeltype="nb")

ggplot()+
  geom_boxplot(aes(x=factor(dbb$site),y=E2))+
  geom_hline(yintercept=0, linetype='dashed', col='blue')+
  theme_bw(16)+ylab("Residuals")+xlab("sites")


##laticce plot
xyplot(E2 ~ month2 | factor(site), 
       data = dbb,
       xlab = list(label = "month (1=oct,7=Apr)", cex = 1.5),
       ylab = list(label = "Pearson residuals", cex = 1.5),
       panel = function(x,y)
       {
         panel.points(x,y, col = 1, pch = 16, cex = 0.7)
         panel.loess(x,y, col = 1, lwd = 2)
         panel.abline(h=0)
       }
)

est<-cbind(coef(m2a),confint(m2a))
exp(est)

###Predicted values
####create database ot plot predicted values over observed data points
newdata<-expand.grid(month2=seq(from=range(dbb$month2)[1], to=range(dbb$month2)[2],length=length(unique(dbb$month2))),
                     site=factor(c(1:6),labels=levels(factor(dbb$site))))
newdata$areacore<-ifelse(newdata$site=="AD",log(unique(dbb$areacore[dbb$site=="AD"])+1),log(unique(dbb$areacore[dbb$site!="AD"])+1))
newdata$x<-predict(m2a,newdata,type="response")

ggplot()+
  geom_point(data=dbb,aes(x=month2, y=x, col=site),  size=2, alpha=0.7,position=position_dodge(width=0.4))+
  geom_path(aes(x=month2, y=x, col=site), data=newdata, size=1)+
  theme_bw()+
  scale_x_continuous(breaks=1:7,labels=c(10,11,12,1,2,3,4))+
  stat_summary(data=dbb,aes(x=month2, y=x, col=site),geom="pointrange",shape=34,size=3,position=position_dodge(width=0.4))+
  ylab("Total Abundance")+xlab("Month")




######## errant polychaetes
dbpe<-DB1[class1=="Polychaeta_errantia"]

m3<-glm.nb(x~site*month2+offset(log(areacore)),data=dbpe)
summary(m3)
anova(m3)
m3a<-update(m3,~.-site:month2)
drop1(m3a,test="Chi")
anova(m3,m3a)
AIC(m3,m3a)
summary(m3)

residualPlots(m3)


est<-cbind(coef(m3),confint(m3))
exp(est)

###Predicted values
####create database ot plot predicted values over observed data points
newdata<-expand.grid(month2=seq(from=range(dbpe$month2)[1], to=range(dbpe$month2)[2],length=length(unique(dbpe$month2))),
                     site=factor(c(1:6),labels=levels(factor(dbpe$site))))
newdata$areacore<-ifelse(newdata$site=="AD",log(unique(dbpe$areacore[dbpe$site=="AD"])+1),log(unique(dbpe$areacore[dbpe$site!="AD"])+1))
newdata$x<-predict(m3,newdata,type="response")

ggplot()+
  geom_point(data=dbpe,aes(x=month2, y=x, col=site),  size=2, alpha=0.7,position=position_dodge(width=0.4))+
  geom_path(aes(x=month2, y=x, col=site), data=newdata, size=1)+
  theme_bw()+
  scale_x_continuous(breaks=1:7,labels=c(10,11,12,1,2,3,4))+
  stat_summary(data=dbpe,aes(x=month2, y=x, col=site),geom="pointrange",shape=34,size=3,position=position_dodge(width=0.4))+
  ylab("Abundance (ind per core)")+xlab("Month")+
  ggtitle("Poly errantia")


### sed poly

dbps<-DB1[class1=="Polychaeta_sedentaria"]

m4<-glm.nb(x~site*month2+offset(log(areacore)),data=dbps)
summary(m4)
anova(m4)
m4a<-update(m4,~.-site:month2)
drop1(m4a,test="Chi")
anova(m4,m4a)
AIC(m4,m4a)
summary(m4)

residualPlots(m4)


est<-cbind(coef(m4),confint(m4))
exp(est)

###Predicted values
####create database ot plot predicted values over observed data points
newdata<-expand.grid(month2=seq(from=range(dbps$month2)[1], to=range(dbps$month2)[2],length=length(unique(dbps$month2))),
                     site=factor(c(1:6),labels=levels(factor(dbps$site))))
newdata$areacore<-ifelse(newdata$site=="AD",log(unique(dbps$areacore[dbps$site=="AD"])+1),log(unique(dbps$areacore[dbps$site!="AD"])+1))
newdata$x<-predict(m4,newdata,type="response")

ggplot()+
  geom_point(data=dbps,aes(x=month2, y=x, col=site),  size=2, alpha=0.7,position=position_dodge(width=0.4))+
  geom_path(aes(x=month2, y=x, col=site), data=newdata, size=1)+
  theme_bw()+
  scale_x_continuous(breaks=1:7,labels=c(10,11,12,1,2,3,4))+
  stat_summary(data=dbps,aes(x=month2, y=x, col=site),geom="pointrange",shape=34,size=3,position=position_dodge(width=0.4))+
  ylab("Abundance (ind per core)")+xlab("Month")+
  ggtitle("Poly sedentaria")

### Gastropoda

######## errant polychaetes
dbg<-DB1[class1=="Gastropoda"]

m5<-glm.nb(x~site*month2+offset(log(areacore)),data=dbg)
summary(m5)
anova(m5)
m5a<-update(m5,~.-site:month2)
drop1(m5a,test="Chi")
anova(m5,m5a)
AIC(m5,m5a)
summary(m5)

residualPlots(m5)


est<-cbind(coef(m5),confint(m5))
exp(est)

###Predicted values
####create database ot plot predicted values over observed data points
newdata<-expand.grid(month2=seq(from=range(dbg$month2)[1], to=range(dbg$month2)[2],length=length(unique(dbg$month2))),
                     site=factor(c(1:6),labels=levels(factor(dbg$site))))
newdata$areacore<-ifelse(newdata$site=="AD",log(unique(dbg$areacore[dbg$site=="AD"])+1),log(unique(dbg$areacore[dbg$site!="AD"])+1))
newdata$x<-predict(m5,newdata,type="response")

ggplot()+
  geom_point(data=dbg,aes(x=month2, y=x, col=site),  size=2, alpha=0.7,position=position_dodge(width=0.4))+
  geom_path(aes(x=month2, y=x, col=site), data=newdata, size=1)+
  theme_bw()+
  scale_x_continuous(breaks=1:7,labels=c(10,11,12,1,2,3,4))+
  stat_summary(data=dbg,aes(x=month2, y=x, col=site),geom="pointrange",shape=34,size=3,position=position_dodge(width=0.4))+
  ylab("Abundance (ind per core)")+xlab("Month")+
  ggtitle("Gastropoda")


#### Malacostraca

dbm<-DB1[class1=="Malacostraca"]

m6<-glm.nb(x~site*month2+offset(log(areacore)),data=dbm)
summary(m6)
anova(m6)
m6a<-update(m6,~.-site:month2)
drop1(m6a,test="Chi")
anova(m6,m6a,test="Chisq")
AIC(m6,m6a)
summary(m6)

residualPlots(m6)


est<-cbind(coef(m6),confint(m6))
exp(est)

###Predicted values
####create database ot plot predicted values over observed data points
newdata<-expand.grid(month2=seq(from=range(dbm$month2)[1], to=range(dbm$month2)[2],length=length(unique(dbm$month2))),
                     site=factor(c(1:6),labels=levels(factor(dbm$site))))
newdata$areacore<-ifelse(newdata$site=="AD",log(unique(dbm$areacore[dbm$site=="AD"])+1),log(unique(dbm$areacore[dbm$site!="AD"])+1))
newdata$x<-predict(m6,newdata,type="response")

ggplot()+
  geom_point(data=dbm,aes(x=month2, y=x, col=site),  size=2, alpha=0.7,position=position_dodge(width=0.4))+
  geom_path(aes(x=month2, y=x, col=site), data=newdata, size=1)+
  theme_bw()+
  scale_x_continuous(breaks=1:7,labels=c(10,11,12,1,2,3,4))+
  stat_summary(data=dbm,aes(x=month2, y=x, col=site),geom="pointrange",shape=34,size=3,position=position_dodge(width=0.4))+
  ylab("Abundance (ind per core)")+xlab("Month")+
  ggtitle("Malacostraca")





