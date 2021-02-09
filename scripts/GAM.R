setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("MASS","ggplot2","viridis","RColorBrewer","psych","reshape2","beepr","data.table","car","mgcv")
lapply(packs,require,character.only=T)


DB68<-fread("data_out/db/Final_DB_class1_density_polyexcl_20210202.csv") ### created in script called Database_cleanup_joining
str(DB68)

time<-c("10","11","12","1","2","3","4")
DB68[,month1:=factor(month,levels=time)]
DB68[,month2:=ifelse(month==10,1,ifelse(month==11,2,ifelse(month==12,3,ifelse(month==1,4,ifelse(month==2,5,ifelse(month==3,6,7))))))]
DB68[,unique(month2)]

DB68<-DB68[-which(DB68$class1=="Other"),]

DB68$site1<-factor(DB68$site,levels=c("A","AB","BI","E","BR","AD"))
DB68$month1<-factor(DB68$month,levels=c("10","11","12","1","2","3","4"))

DB68$month2<-ifelse(DB68$month1==10,1,ifelse(DB68$month1==11,2,ifelse(DB68$month1==12,3,ifelse(DB68$month1==1,4,
                    ifelse(DB68$month1==2,5,ifelse(DB68$month1==3,6,ifelse(DB68$month1==4,7,NA)))))))

DB68[,areacore:=ifelse(site=="AD",0.0113,0.00785)]


DBAD<-DB68[site1=="AD"]
DBAD_PS<-DBAD[class1=="Polychaeta_sedentaria"]
DBAD_PE<-DBAD[class1=="Polychaeta_errantia"]
DBAD_B<-DBAD[class1=="Bivalvia"]
DBAD_G<-DBAD[class1=="Gastropoda"]
DBAD_M<-DBAD[class1=="Malacostraca"]

DBA<-DB68[site1=="A"]
DBA_PS<-DBA[class1=="Polychaeta_sedentaria"]
DBA_PE<-DBA[class1=="Polychaeta_errantia"]
DBA_B<-DBA[class1=="Bivalvia"]
DBA_G<-DBA[class1=="Gastropoda"]
DBA_M<-DBA[class1=="Malacostraca"]

DBAB<-DB68[site1=="AB"]
DBAB_PS<-DBAB[class1=="Polychaeta_sedentaria"]
DBAB_PE<-DBAB[class1=="Polychaeta_errantia"]
DBAB_B<-DBAB[class1=="Bivalvia"]
DBAB_G<-DBAB[class1=="Gastropoda"]
DBAB_M<-DBAB[class1=="Malacostraca"]

DBBI<-DB68[site1=="BI"]
DBBI_PS<-DBBI[class1=="Polychaeta_sedentaria"]
DBBI_PE<-DBBI[class1=="Polychaeta_errantia"]
DBBI_B<-DBBI[class1=="Bivalvia"]
DBBI_G<-DBBI[class1=="Gastropoda"]
DBBI_M<-DBBI[class1=="Malacostraca"]

DBBR<-DB68[site1=="BR"]
DBBR_PS<-DBBR[class1=="Polychaeta_sedentaria"]
DBBR_PE<-DBBR[class1=="Polychaeta_errantia"]
DBBR_B<-DBBR[class1=="Bivalvia"]
DBBR_G<-DBBR[class1=="Gastropoda"]
DBBR_M<-DBBR[class1=="Malacostraca"]

DBE<-DB68[site1=="E"]
DBE_PS<-DBE[class1=="Polychaeta_sedentaria"]
DBE_PE<-DBE[class1=="Polychaeta_errantia"]
DBE_B<-DBE[class1=="Bivalvia"]
DBE_G<-DBE[class1=="Gastropoda"]
DBE_M<-DBE[class1=="Malacostraca"]


g0<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DB68,gamma=1.4)
plot.gam(g0,residuals=F,pch=16,all.terms=T)

ggplot(DB68,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  #geom_line(aes(y=fitted(g0)),color = "green", size=1.2) +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = F, color="blue") +
  theme_bw()

summary(g0)

#Anrumai

##PS
g1<-gam(x~s(month2,k=6)+offset(log(areacore)),family="nb",data=DBA_PS,method="REML")
g1a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBA_PS,method="REML")

summary(g1)
anova(g1)

summary(g1a)
anova(g1,g1a,test="Chisq")

AIC(g1)
AIC(g1a)

summary(g1)$sp.criterion
summary(g1a)$sp.criterion

summary(g1)$r.sq
summary(g1a)$r.sq

gam.check(g1, k.rep=1000)


plot.gam(g1,residuals=T,pch=16,all.terms=T)

ggplot(DBA_PS,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  #geom_line(aes(y=fitted(g1)),color = "green", size=1.2) +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()

## PE
g2<-gam(x~s(month2,k=6)+offset(log(areacore)),family="nb",data=DBA_PE) ##better model
g2a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBA_PE)

summary(g2)
summary(g2a)

anova(g2,g2a,test="Chisq")

AIC(g2,g2a)

summary(g2)$sp.criterion
summary(g2a)$sp.criterion

gam.check(g2, k.rep=1000)

plot.gam(g2,residuals=T,pch=16,all.terms=T)

ggplot(DBA_PE,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()

## B
g3<-gam(x~s(month2,k=6)+offset(log(areacore)),family="nb",data=DBA_B) ##better model
g3a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBA_B)

summary(g3)
summary(g3a)

anova(g3,g3a,test="Chisq")

AIC(g3,g3a)

summary(g3)$sp.criterion
summary(g3a)$sp.criterion

gam.check(g3, k.rep=1000)

plot.gam(g3,residuals=T,pch=16,all.terms=T)

ggplot(DBA_B,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()


## G
g4<-gam(x~s(month2,k=6)+offset(log(areacore)),family="nb",data=DBA_G) ##better model
g4a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBA_G)

summary(g4)
summary(g4a)

anova(g4,g4a,test="Chisq")

AIC(g4,g4a)

summary(g4)$sp.criterion
summary(g4a)$sp.criterion

gam.check(g4, k.rep=1000)

plot.gam(g4,residuals=T,pch=16,all.terms=T)

ggplot(DBA_G,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()


## M
g5<-gam(x~s(month2,k=6)+offset(log(areacore)),family="nb",data=DBA_M) ##better model
g5a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBA_M)

summary(g5)
summary(g5a)

anova(g5,g5a,test="Chisq")

AIC(g5,g5a)

summary(g5)$sp.criterion
summary(g5a)$sp.criterion

gam.check(g5, k.rep=1000)

plot.gam(g5,residuals=T,pch=16,all.terms=T)

ggplot(DBA_M,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()


#################ABU

##PS
g1<-gam(x~s(month2,k=6)+offset(log(areacore)),family="nb",data=DBAB_PS,method="REML")
g1a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBAB_PS,method="REML")

summary(g1)
summary(g1a)

anova(g1,g1a,test="Chisq")

AIC(g1)
AIC(g1a)

summary(g1)$sp.criterion
summary(g1a)$sp.criterion

gam.check(g1, k.rep=1000)

plot.gam(g1,residuals=T,pch=16,all.terms=T)

ggplot(DBAB_PS,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_line(aes(y=fitted(g1)),color = "red", size=1.2) +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()

## PE
g2<-gam(x~s(month2,k=6)+offset(log(areacore)),family="nb",data=DBAB_PE) ##better model
g2a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBAB_PE)

summary(g2)
summary(g2a)

anova(g2,g2a,test="Chisq")

AIC(g2,g2a)

summary(g2)$sp.criterion
summary(g2a)$sp.criterion

gam.check(g2, k.rep=1000)

plot.gam(g2,residuals=T,pch=16,all.terms=T)

ggplot(DBAB_PE,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()

## B
g3<-gam(x~s(month2,k=6)+offset(log(areacore)),family="nb",data=DBAB_B,method="REML") ##better model
g3a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBAB_B,method="REML")

summary(g3)
summary(g3a)

anova(g3,g3a,test="Chisq")

AIC(g3,g3a)

summary(g3)$sp.criterion
summary(g3a)$sp.criterion

gam.check(g3, k.rep=1000)

plot.gam(g3,residuals=T,pch=16,all.terms=T)

ggplot(DBAB_B,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  geom_line(aes(y=fitted(g3)),color = "red", size=1.2) +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()


## G
g4<-gam(x~s(month2,k=6)+offset(log(areacore)),family="nb",data=DBAB_G) ##better model
g4a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBAB_G)

summary(g4)
summary(g4a)

anova(g4,g4a,test="Chisq")

AIC(g4,g4a)

summary(g4)$sp.criterion
summary(g4a)$sp.criterion

gam.check(g4, k.rep=1000)

plot.gam(g4,residuals=T,pch=16,all.terms=T)

ggplot(DBAB_G,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()


## M
g5<-gam(x~s(month2,k=6)+offset(log(areacore)),family="nb",data=DBAB_M) ##better model
g5a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBAB_M)

summary(g5)
summary(g5a)

anova(g5,g5a,test="Chisq")

AIC(g5,g5a)

summary(g5)$sp.criterion
summary(g5a)$sp.criterion

gam.check(g5, k.rep=1000)

plot.gam(g5,residuals=T,pch=16,all.terms=T)

ggplot(DBAB_M,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()


##Bijante

##PS
g1<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DBBI_PS)
g1a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBBI_PS)

summary(g1)
summary(g1a)

anova(g1,g1a,test="Chisq")

AIC(g1)
AIC(g1a)

summary(g1)$sp.criterion
summary(g1a)$sp.criterion

gam.check(g1, k.rep=1000)

plot.gam(g1,residuals=T,pch=16,all.terms=T)

ggplot(DBBI_PS,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  #geom_line(aes(y=fitted(g1)),color = "green", size=1.2) +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()

## PE
g2<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DBBI_PE) ##better model
g2a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBBI_PE)

summary(g2)
summary(g2a)

anova(g2,g2a,test="Chisq")

AIC(g2,g2a)

summary(g2)$sp.criterion
summary(g2a)$sp.criterion

gam.check(g2, k.rep=1000)

plot.gam(g2,residuals=T,pch=16,all.terms=T)

ggplot(DBBI_PE,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()

## B
g3<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DBBI_B) ##better model
g3a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBBI_B)

summary(g3)
summary(g3a)

anova(g3,g3a,test="Chisq")

AIC(g3,g3a)

summary(g3)$sp.criterion
summary(g3a)$sp.criterion

gam.check(g3, k.rep=1000)

plot.gam(g3,residuals=T,pch=16,all.terms=T)

ggplot(DBBI_B,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()


## G
g4<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DBBI_G) ##better model
g4a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBBI_G)

summary(g4)
summary(g4a)

anova(g4,g4a,test="Chisq")

AIC(g4,g4a)

summary(g4)$sp.criterion
summary(g4a)$sp.criterion

gam.check(g4, k.rep=1000)

plot.gam(g4,residuals=T,pch=16,all.terms=T)

ggplot(DBBI_G,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()


## M
g5<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DBBI_M) ##better model
g5a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBBI_M)

summary(g5)
summary(g5a)

anova(g5,g5a,test="Chisq")

AIC(g5,g5a)

summary(g5)$sp.criterion
summary(g5a)$sp.criterion

gam.check(g5, k.rep=1000)

plot.gam(g5,residuals=T,pch=16,all.terms=T)

ggplot(DBBI_M,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()


###Bruce
##PS
g1<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DBBR_PS)
g1a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBBR_PS)

summary(g1)
summary(g1a)

anova(g1,g1a,test="Chisq")

AIC(g1)
AIC(g1a)

summary(g1)$sp.criterion
summary(g1a)$sp.criterion

gam.check(g1, k.rep=1000)

plot.gam(g1,residuals=T,pch=16,all.terms=T)

ggplot(DBBR_PS,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  #geom_line(aes(y=fitted(g1)),color = "green", size=1.2) +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()

## PE
g2<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DBBR_PE) ##better model
g2a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBBR_PE)

summary(g2)
summary(g2a)

anova(g2,g2a,test="Chisq")

AIC(g2,g2a)

summary(g2)$sp.criterion
summary(g2a)$sp.criterion

gam.check(g2, k.rep=1000)

plot.gam(g2,residuals=T,pch=16,all.terms=T)

ggplot(DBBR_PE,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()

## B
g3<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DBBR_B) ##better model
g3a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBBR_B)

summary(g3)
summary(g3a)

anova(g3,g3a,test="Chisq")

AIC(g3,g3a)

summary(g3)$sp.criterion
summary(g3a)$sp.criterion

gam.check(g3, k.rep=1000)

plot.gam(g3,residuals=T,pch=16,all.terms=T)

ggplot(DBBR_B,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()


## G
g4<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DBBR_G) ##better model
g4a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBBR_G)

summary(g4)
summary(g4a)

anova(g4,g4a,test="Chisq")

AIC(g4,g4a)

summary(g4)$sp.criterion
summary(g4a)$sp.criterion

gam.check(g4, k.rep=1000)

plot.gam(g4,residuals=T,pch=16,all.terms=T)

ggplot(DBBR_G,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()


## M
g5<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DBBR_M) ##better model
g5a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBBR_M)

summary(g5)
summary(g5a)

anova(g5,g5a,test="Chisq")

AIC(g5,g5a)

summary(g5)$sp.criterion
summary(g5a)$sp.criterion

gam.check(g5, k.rep=1000)

plot.gam(g5,residuals=T,pch=16,all.terms=T)

ggplot(DBBR_M,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()



####Escadinhas
##PS
g1<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DBE_PS)
g1a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBE_PS)

summary(g1)
summary(g1a)

anova(g1,g1a,test="Chisq")

AIC(g1)
AIC(g1a)

summary(g1)$sp.criterion
summary(g1a)$sp.criterion

gam.check(g1, k.rep=1000)

plot.gam(g1,residuals=T,pch=16,all.terms=T)

ggplot(DBE_PS,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  #geom_line(aes(y=fitted(g1)),color = "green", size=1.2) +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()

## PE
g2<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DBE_PE) ##better model
g2a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBE_PE)

summary(g2)
summary(g2a)

anova(g2,g2a,test="Chisq")

AIC(g2,g2a)

summary(g2)$sp.criterion
summary(g2a)$sp.criterion

gam.check(g2, k.rep=1000)

plot.gam(g2,residuals=T,pch=16,all.terms=T)

ggplot(DBE_PE,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()

## B
g3<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DBE_B) ##better model
g3a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBE_B)

summary(g3)
summary(g3a)

anova(g3,g3a,test="Chisq")

AIC(g3,g3a)

summary(g3)$sp.criterion
summary(g3a)$sp.criterion

gam.check(g3, k.rep=1000)

plot.gam(g3,residuals=T,pch=16,all.terms=T)

ggplot(DBE_B,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=7), se = T, color="blue") +
  geom_line(data=predict(g3),se = T, color="red") +
  theme_bw()


## G
g4<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DBE_G) ##better model
g4a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBE_G)

summary(g4)
summary(g4a)

anova(g4,g4a,test="Chisq")

AIC(g4,g4a)

summary(g4)$sp.criterion
summary(g4a)$sp.criterion

gam.check(g4, k.rep=1000)

plot.gam(g4,residuals=T,pch=16,all.terms=T)

ggplot(DBE_G,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()


## M
g5<-gam(x~s(month2,k=7)+offset(log(areacore)),family="nb",data=DBE_M) ##better model
g5a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBE_M)

summary(g5)
summary(g5a)

anova(g5,g5a,test="Chisq")

AIC(g5,g5a)

summary(g5)$sp.criterion
summary(g5a)$sp.criterion

gam.check(g5, k.rep=1000)

plot.gam(g5,residuals=T,pch=16,all.terms=T)

ggplot(DBE_M,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=7), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()



########Adonga
##PS
g1<-gam(x~s(month2,k=3)+offset(log(areacore)),family="nb",data=DBAD_PS)
g1a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBAD_PS)

summary(g1)
summary(g1a)

anova(g1,g1a,test="Chisq")

AIC(g1)
AIC(g1a)

summary(g1)$sp.criterion
summary(g1a)$sp.criterion

gam.check(g1, k.rep=1000)

plot.gam(g1,residuals=T,pch=16,all.terms=T)

ggplot(DBAD_PS,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  #geom_line(aes(y=fitted(g1)),color = "green", size=1.2) +
  geom_smooth(method="gam", formula=y~s(x,k=6), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()

## PE
g2<-gam(x~s(month2,k=3)+offset(log(areacore)),family="nb",data=DBAD_PE) ##better model
g2a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBAD_PE)

summary(g2)
summary(g2a)

anova(g2,g2a,test="Chisq")

AIC(g2,g2a)

summary(g2)$sp.criterion
summary(g2a)$sp.criterion

gam.check(g2, k.rep=1000)

plot.gam(g2,residuals=T,pch=16,all.terms=T)

ggplot(DBAD_PE,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=3), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()

## B
g3<-gam(x~s(month2,k=3)+offset(log(areacore)),family="nb",data=DBAD_B) ##better model
g3a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBAD_B)

summary(g3)
summary(g3a)

anova(g3,g3a,test="Chisq")

AIC(g3,g3a)

summary(g3)$sp.criterion
summary(g3a)$sp.criterion

gam.check(g3, k.rep=1000)

plot.gam(g3,residuals=T,pch=16,all.terms=T)

ggplot(DBAD_B,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=3), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()


## G
g4<-gam(x~s(month2,k=3)+offset(log(areacore)),family="nb",data=DBAD_G) ##better model
g4a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBAD_G)

summary(g4)
summary(g4a)

anova(g4,g4a,test="Chisq")

AIC(g4,g4a)

summary(g4)$sp.criterion
summary(g4a)$sp.criterion

gam.check(g4, k.rep=1000)

plot.gam(g4,residuals=T,pch=16,all.terms=T)

ggplot(DBAD_G,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=3), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()


## M
g5<-gam(x~s(month2,k=3)+offset(log(areacore)),family="nb",data=DBAD_M) ##better model
g5a<-gam(x~month2+offset(log(areacore)),family="nb",data=DBAD_M)

summary(g5)
summary(g5a)

anova(g5,g5a,test="Chisq")

AIC(g5,g5a)

summary(g5)$sp.criterion
summary(g5a)$sp.criterion

gam.check(g5, k.rep=1000)

plot.gam(g5,residuals=T,pch=16,all.terms=T)

ggplot(DBAD_M,aes(x=month2,y=x)) +
  #geom_point() +
  stat_summary()+
  labs(x = "Month (1=Oct,7=Apr)", y = "Density (ind/m2)") +
  geom_smooth(method="gam", formula=y~s(x,k=3), se = T, color="blue") +
  #geom_smooth(method=lm, se = T, color="red") +
  theme_bw()

