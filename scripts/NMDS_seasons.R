setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("vegan","ggplot2","viridis","RColorBrewer","psych","reshape2","beepr","data.table","BiodiversityR","cluster","pairwiseAdonis")
lapply(packs,require,character.only=T)

p<-c("#4DAF4A","darkgreen","#377EB8","#984EA3","#FF7F00","#E41A1C")

DB1<-fread("Data_out/db/DB_multianal_20210312.csv") ##Data created in script Bray-Curtis_NMDS

##Prepare data and remove rows with zeros in all columns

data<-DB1[apply(DB1[,!(1:3)],1,sum)!=0] ### Bray curtis does not work with empty cores. So we have to remove all of them
table(apply(data[,!c(1:3)], 2, sum)==0)

data[,season:="season"][month==10|month==11,season:="begining"][month==12|month==1|month==2,season:="mid"][month==3|month==4,season:="end"]
data[,season:=factor(season, levels=c("begining","mid","end"))]
data[,.(season,month)]

data1<-data[,!c(1:3,86)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables


data1[,dummy:=88.49558]

data2<-data[,c(1:3,86)] ###store aggregating variables to use latter

############Cutting db
m<-data1[,apply(.SD,2,mean)]
str(m)
names(m)

m1<-data.table(sp=names(m),m)

totil<-sum(m1$m)

m1[,freq:=m/totil*100]
m1[,sum(freq)]

m2<-m1[order(-freq)]

sum(m2[1:22,"freq"]) ###chosing taxa totaling 95% of overall mean density

###Creating cuts for db to reduce less abundant sps
m2[,p95:="N"][1:35,p95:="Y"]
m2[,p75:="N"][1:14,p75:="Y"]
m2[,p1p:="N"][1:22,p1p:="Y"]


####### 75%
nam<-m2[p75=="N",sp]
str(nam)
data75<-data[,!..nam]
table(is.na(data75))

##Remove cores with zero species after selection of p75
data75_z<-data75[apply(data75[,!c(1:3,17)],1,sum)!=0]
data75_1<-data75_z[,!c(1:3,17)]
data75_1[,dummy:=88.49558]

data75_2<-data75_z[,c(1:3,17)]

###log transform
data1log<-log(data1+1)
data75_1log<-log(data75_1+1)

#### run NMDS
set.seed(0)
NMDS<-metaMDS(data75_1log,distance="bray",k=3,trymax=1000,autotransform = F)
beep()

stressplot(NMDS)

#extract NMDS scores (x and y coordinates)
data.scores=as.data.table(scores(NMDS))

#add columns to data frame 
data.scores$coreID = data75_2$coreID
data.scores$site = data75_2$site
data.scores$season = data75_2$season

head(data.scores)

#p<-c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","black")

PP<-ggplot(data.scores,aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=season))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=season, colour=season))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS, with log, with dummy")+
  theme_bw()

###fit sps data in Adonga
fit<-envfit(NMDS,data75_1,permutations=999)
arrow<-data.frame(fit$vectors$arrows,R = fit$vectors$r, P = fit$vectors$pvals)
arrow$FG <- rownames(arrow)
arrow.p<-arrow[arrow$P<=0.05,]

PP+
  geom_segment(data=arrow.p,aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrow.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


###Permanova

dist<-vegdist(data75_1log,method="bray")
per1<-adonis(dist~site*season,data=data75_2,permutations=999)


dispersion<-betadisper(dist,group=data75_2$site)
permutest(dispersion)

plot(dispersion,hull=F,ellipse=T,lwd=1)

set.seed(200)
pairwise.adonis(data75_1log,factors=data75_2$season,p.adjust.m="bonferroni", perm=100000)


################OVERALL SPATIAL AND TEMPORAL######### CANONICAL 

### CAP Total
###
data2[,site:=factor(site,levels=c("A","AB","BI","E","BR","AD"))]

set.seed(118)
system.time(
  OM1Tot<-CAPdiscrim(data1log~site*season,data=data2,dist="bray",axes=3,m=0,add=F, permutations = 1000)
)
beep(3)
summary(OM1Tot)
OM1Tot$manova
OM1Tot$m

##PERMANOVA
perTot<-adonis(data1log~site*season,data=data2,permutations = 1000)
perTot

distTot<-vegdist(data1log,method="bray")

set.seed(1000)
dispTot<-betadisper(distTot,group=data2$site)
permutest(dispTot)
boxplot(dispTot)

TukTot<-TukeyHSD(dispTot)
plot(TukTot)



#extract CAP scores (x and y coordinates)
data.scoresCAP_Tot = as.data.table(scores(OM1Tot))

#add columns to data frame 
data.scoresCAP_Tot$site = data2$site
data.scoresCAP_Tot$season = data2$season
head(data.scoresCAP_Tot)


###fit sps data in Tot
set.seed(119)
fitCAP_Tot<-envfit(OM1Tot,data1log,permutations=100000)
arrowCAP_Tot<-data.frame(fitCAP_Tot$vectors$arrows,R=fitCAP_Tot$vectors$r,P=fitCAP_Tot$vectors$pvals)
arrowCAP_Tot$FG <- rownames(arrowCAP_Tot)
arrowarrowCAP_Tot.p<-arrowCAP_Tot[arrowCAP_Tot$P<=0.05&arrowCAP_Tot$R>=0.15,]

##Site
PP_CAP_site<-ggplot(data.scoresCAP_Tot,aes(x=LD1,y=LD2))+
  geom_point(size=2,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x="LD1",y="LD2",colour="Site",fill="Site")+
  stat_ellipse(size=2,type="t",aes(group=site,colour=site,fill=site),level=.6,geom="polygon",alpha=.2)+
  scale_colour_manual(values=p)+
  scale_fill_manual(values=p)+
  ggtitle("Canonical Discriminant per site")+
  theme_bw()
PP_CAP_site+
  geom_segment(data=arrowarrowCAP_Tot.p,aes(x=0,y=0,xend=LD1*R*7,yend=LD2*R*7),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowarrowCAP_Tot.p,aes(x=LD1*R*7,y=LD2*R*7,label=FG),cex=5,direction="both",segment.size=0.25)

##season
PP_CAP_S<-ggplot(data.scoresCAP_Tot,aes(x=LD1,y=LD2))+
  geom_point(size=4,aes(colour=season))+
  #geom_text(aes(label=month),size=7)+
  labs(x="LD1",y="LD2",colour="Season",fill="Season")+
  stat_ellipse(size=2,type="t",aes(group=season,colour=season,fill=season),level=.6,geom="polygon",alpha=.2)+
  #scale_colour_manual(values=p)+
  #scale_fill_manual(values=p)+
  ggtitle("Canonical Discriminant per season")+
  theme_bw()
PP_CAP_S+
  geom_segment(data=arrowarrowCAP_Tot.p,aes(x=0,y=0,xend=LD1*R*7,yend=LD2*R*7),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowarrowCAP_Tot.p,aes(x=LD1*R*7,y=LD2*R*7,label=FG),cex=5,direction="both",segment.size=0.25)







######Anrumai

dataA<-data[site=="A"]
dataA1<-dataA[,!c(1:3,89)]
dataA1[,dummy:=127.3885]

dataA2<-dataA[,c(1:3,89)]

dataA1log<-log(dataA1+1)


#### run NMDS
set.seed(1)
NMDSA<-metaMDS(dataA1log,distance="bray",k=3,trymax=1000,autotransform = F)
beep()

#stressplot(NMDSA)

#extract NMDS scores (x and y coordinates)
data.scoresA=as.data.table(scores(NMDSA))

#add columns to data frame 
data.scoresA$coreID = dataA2$coreID
data.scoresA$site = dataA2$site
data.scoresA$season = factor(dataA2$season)

head(data.scoresA)


###fit sps data in Adonga
fitA<-envfit(NMDSA,dataA1log,permutations=999)
arrowA<-data.frame(fitA$vectors$arrows,R = fitA$vectors$r, P = fitA$vectors$pvals)
arrowA$FG <- rownames(arrowA)
arrowA.p<-arrowA[arrowA$P<=0.05,]

PPA<-ggplot(data.scoresA,aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=season))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Season",fill="Season")+
  stat_ellipse(size=2, type="t",aes(group=season, colour=season,fill=season),level=.6,geom="polygon",alpha=.2)+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Anrumai, with log, with dummy")+
  theme_bw()
PPA+
  geom_segment(data=arrowA.p,aes(x=0,y=0,xend=NMDS1*R*2.5,yend=NMDS2*R*2.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowA.p,aes(x=NMDS1*R*2.5,y=NMDS2*R*2.5,label=FG),cex=5,direction="both",segment.size=0.25)

###Permanova

distA<-vegdist(dataA1log,method="bray")
perA<-adonis(distA~season,data=dataA2,permutations=999)


dispersionA<-betadisper(distA,group=dataA2$season)
permutest(dispersionA)

plot(dispersionA,hull=F,ellipse=T,lwd=1)

set.seed(200)
pairwise.adonis(distA,factors=dataA2$season,p.adjust.m="bonferroni", perm=100000)

### CAP Anrumai
###

system.time(
  OM1<-CAPdiscrim(dataA1log~season,data=dataA2,dist="bray",axes=3,m=0,add=F, permutations = 2000)
)
beep(3)
summary(OM1)
OM1$manova

#extract CAP scores (x and y coordinates)
data.scoresCAPA = as.data.table(scores(OM1))

#add columns to data frame 
data.scoresCAPA$season = dataA2$season
head(data.scoresCAPA)


###fit sps data in Adonga
fitCAPA<-envfit(OM1,dataA1log,permutations=999)
arrowCAPA<-data.frame(fitCAPA$vectors$arrows,R=fitCAPA$vectors$r,P=fitCAPA$vectors$pvals)
arrowCAPA$FG <- rownames(arrowCAPA)
arrowCAPA.p<-arrowCAPA[arrowCAPA$P<=0.05,]

PP_CAPA<-ggplot(data.scoresCAPA,aes(x = LD1, y = LD2)) + 
  geom_point(size = 4,aes(colour=season))+
  #geom_text(aes(label=month),size=7)+
  labs(x="LD1",y="LD2",colour="Season",fill="Season")+
  stat_ellipse(size=2, type="t",aes(group=season, colour=season,fill=season),level=.6,geom="polygon",alpha=.2)+
  #scale_colour_manual(values=p)+
  ggtitle("LD in Anrumai, with log, with dummy")+
  theme_bw()
PP_CAPA+
  geom_segment(data=arrowCAPA.p,aes(x=0,y=0,xend=LD1*R*6,yend=LD2*R*6),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowCAPA.p,aes(x=LD1*R*6,y=LD2*R*6,label=FG),cex=5,direction="both",segment.size=0.25)



######Abu

dataAB<-data[site=="AB"]
dataAB1<-dataAB[,!c(1:3,89)]
dataAB1[,dummy:=88.49558]

dataAB2<-dataAB[,c(1:3,89)]

dataAB1log<-log(dataAB1+1)


#### run NMDS
set.seed(1)
NMDSAB<-metaMDS(dataAB1log,distance="bray",k=3,trymax=1000,autotransform = F)
beep(3)

#stressplot(NMDSAB)

#extract NMDS scores (x and y coordinates)
data.scoresAB=as.data.table(scores(NMDSAB))

#add columns to data frame 
data.scoresAB$coreID = dataAB2$coreID
data.scoresAB$site = dataAB2$site
data.scoresAB$season = factor(dataAB2$season)

head(data.scoresAB)

###fit sps data in Abu
fitAB<-envfit(NMDSAB,dataAB1log,permutations=999)
arrowAB<-data.frame(fitAB$vectors$arrows,R = fitAB$vectors$r, P = fitAB$vectors$pvals)
arrowAB$FG <- rownames(arrowAB)
arrowAB.p<-arrowAB[arrowAB$P<=0.05,]

PPAB<-ggplot(data.scoresAB,aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=season))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Season",fill="Season")+
  stat_ellipse(size=2, type="t",aes(group=season, colour=season,fill=season),level=.6,geom="polygon",alpha=.15)+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Abu, with log, with dummy")+
  theme_bw()
PPAB+
  geom_segment(data=arrowAB.p,aes(x=0,y=0,xend=NMDS1*R*2.5,yend=NMDS2*R*2.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowAB.p,aes(x=NMDS1*R*2.5,y=NMDS2*R*2.5,label=FG),cex=5,direction="both",segment.size=0.25)

###Permanova

distAB<-vegdist(dataAB1log,method="bray")
perAB<-adonis(distAB~season,data=dataAB2,permutations=999)

set.seed(150)
dispersionAB<-betadisper(distAB,group=dataAB2$season)
permutest(dispersionAB)

plot(dispersionAB,hull=F,ellipse=T,lwd=1)

set.seed(201)
pairwise.adonis(distAB,factors=dataAB2$season,p.adjust.m="bonferroni", perm=100000)


### CAP Abu
###
system.time(
  OM1AB<-CAPdiscrim(dataAB1log~season,data=dataAB2,dist="bray",axes=3,m=0,add=F, permutations = 1000)
)
beep(3)
summary(OM1AB)
OM1AB$manova

#extract CAP scores (x and y coordinates)
data.scoresCAP_AB = as.data.table(scores(OM1AB))

#add columns to data frame 
data.scoresCAP_AB$season = dataAB2$season
head(data.scoresCAP_AB)


###fit sps data in Adonga
fitCAP_AB<-envfit(OM1AB,dataAB1log,permutations=999)
arrowCAP_AB<-data.frame(fitCAP_AB$vectors$arrows,R=fitCAP_AB$vectors$r,P=fitCAP_AB$vectors$pvals)
arrowCAP_AB$FG <- rownames(arrowCAP_AB)
arrowarrowCAP_AB.p<-arrowCAP_AB[arrowCAP_AB$P<=0.05&arrowCAP_AB$R>=0.15,]

PP_CAP_AB<-ggplot(data.scoresCAP_AB,aes(x=LD1,y=LD2))+
  geom_point(size=4,aes(colour=season))+
  #geom_text(aes(label=month),size=7)+
  labs(x="LD1",y="LD2",colour="Season",fill="Season")+
  stat_ellipse(size=2,type="t",aes(group=season,colour=season,fill=season),level=.6,geom="polygon",alpha=.2)+
  #scale_colour_manual(values=p)+
  ggtitle("LD in Abu, with log, with dummy")+
  theme_bw()
PP_CAP_AB+
  geom_segment(data=arrowarrowCAP_AB.p,aes(x=0,y=0,xend=LD1*R*6,yend=LD2*R*6),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowarrowCAP_AB.p,aes(x=LD1*R*6,y=LD2*R*6,label=FG),cex=5,direction="both",segment.size=0.25)






######Bijante

dataBI<-data[site=="BI"]
dataBI1<-dataBI[,!c(1:3,89)]
dataBI1[,dummy:=127.3885]

dataBI2<-dataBI[,c(1:3,89)]

dataBI1log<-log(dataBI1+1)


#### run NMDS
set.seed(1)
NMDSBI<-metaMDS(dataBI1log,distance="bray",k=3,trymax=1000,autotransform = F)
beep()

#stressplot(NMDSBI)

#extract NMDS scores (x and y coordinates)
data.scoresBI=as.data.table(scores(NMDSBI))

#add columns to data frame 
data.scoresBI$coreID = dataBI2$coreID
data.scoresBI$site = dataBI2$site
data.scoresBI$season = factor(dataBI2$season)

head(data.scoresBI)

###fit sps data in Bijante
fitBI<-envfit(NMDSBI,dataBI1log,permutations=999)
arrowBI<-data.frame(fitBI$vectors$arrows,R = fitBI$vectors$r, P = fitBI$vectors$pvals)
arrowBI$FG <- rownames(arrowBI)
arrowBI.p<-arrowBI[arrowBI$P<=0.05&arrowBI$R>=.2,]

PPBI<-ggplot(data.scoresBI,aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=season))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Season",fill="Season")+
  stat_ellipse(size=2, type="t",aes(group=season, colour=season,fill=season),level=.6,geom="polygon",alpha=.15)+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Bijante, with log, with dummy")+
  theme_bw()
PPBI+
  geom_segment(data=arrowBI.p,aes(x=0,y=0,xend=NMDS1*R*2.5,yend=NMDS2*R*2.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowBI.p,aes(x=NMDS1*R*2.5,y=NMDS2*R*2.5,label=FG),cex=5,direction="both",segment.size=0.25)

###Permanova

distBI<-vegdist(dataBI1log,method="bray")
perBI<-adonis(distBI~season,data=dataBI2,permutations=999)
perBI

set.seed(401)
dispersionBI<-betadisper(distBI,group=dataBI2$season)
permutest(dispersionBI)
TukeyHSD(dispersionBI)

plot(dispersionBI,hull=T,ellipse=F,lwd=1)

set.seed(400)
pairwise.adonis(distBI,factors=dataBI2$season,p.adjust.m="bonferroni", perm=100000)


### CAP Bijante
###
system.time(
  OM1BI<-CAPdiscrim(dataBI1log~season,data=dataBI2,dist="bray",axes=3,m=0,add=F, permutations = 1000)
)
beep()
summary(OM1BI)
OM1BI$manova

#extract CAP scores (x and y coordinates)
data.scoresCAP_BI = as.data.table(scores(OM1BI))

#add columns to data frame 
data.scoresCAP_BI$season = dataBI2$season
head(data.scoresCAP_BI)


###fit sps data in Adonga
fitCAP_BI<-envfit(OM1BI,dataBI1log,permutations=999)
arrowCAP_BI<-data.frame(fitCAP_BI$vectors$arrows,R=fitCAP_BI$vectors$r,P=fitCAP_BI$vectors$pvals)
arrowCAP_BI$FG <- rownames(arrowCAP_BI)
arrowarrowCAP_BI.p<-arrowCAP_BI[arrowCAP_BI$P<=0.05&arrowCAP_BI$R>=0.15,]

PP_CAP_BI<-ggplot(data.scoresCAP_BI,aes(x=LD1,y=LD2))+
  geom_point(size=4,aes(colour=season))+
  #geom_text(aes(label=month),size=7)+
  labs(x="LD1",y="LD2",colour="Season",fill="Season")+
  stat_ellipse(size=2,type="t",aes(group=season,colour=season,fill=season),level=.6,geom="polygon",alpha=.2)+
  #scale_colour_manual(values=p)+
  ggtitle("LD in Bijante, with log, with dummy")+
  theme_bw()
PP_CAP_BI+
  geom_segment(data=arrowarrowCAP_BI.p,aes(x=0,y=0,xend=LD1*R*6,yend=LD2*R*6),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowarrowCAP_BI.p,aes(x=LD1*R*6,y=LD2*R*6,label=FG),cex=5,direction="both",segment.size=0.25)


######Bruce

dataBR<-data[site=="BR"]
dataBR1<-dataBR[,!c(1:3,89)]
dataBR1[,dummy:=127.3885]

dataBR2<-dataBR[,c(1:3,89)]

dataBR1log<-log(dataBR1+1)


#### run NMDS
set.seed(1)
NMDSBR<-metaMDS(dataBR1log,distance="bray",k=3,trymax=1000,autotransform = F)
beep()

#stressplot(NMDSBR)

#extract NMDS scores (x and y coordinates)
data.scoresBR=as.data.table(scores(NMDSBR))

#add columns to data frame 
data.scoresBR$coreID = dataBR2$coreID
data.scoresBR$site = dataBR2$site
data.scoresBR$season = factor(dataBR2$season)

head(data.scoresBR)

###fit sps data in Bruce
fitBR<-envfit(NMDSBR,dataBR1log,permutations=999)
arrowBR<-data.frame(fitBR$vectors$arrows,R = fitBR$vectors$r, P = fitBR$vectors$pvals)
arrowBR$FG <- rownames(arrowBR)
arrowBR.p<-arrowBR[arrowBR$P<=0.05&arrowBR$R>=.2,]

PPBR<-ggplot(data.scoresBR,aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=season))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Season",fill="Season")+
  stat_ellipse(size=2, type="t",aes(group=season, colour=season,fill=season),level=.6,geom="polygon",alpha=.15)+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Bruce, with log, with dummy")+
  theme_bw()
PPBR+
  geom_segment(data=arrowBR.p,aes(x=0,y=0,xend=NMDS1*R*2.5,yend=NMDS2*R*2.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowBR.p,aes(x=NMDS1*R*2.5,y=NMDS2*R*2.5,label=FG),cex=5,direction="both",segment.size=0.25)

###Permanova
set.seed(501)
distBR<-vegdist(dataBR1log,method="bray")
perBR<-adonis(distBR~season,data=dataBR2,permutations=999)
perBR

set.seed(502)
dispersionBR<-betadisper(distBR,group=dataBR2$season)
permutest(dispersionBR)
TukeyHSD(dispersionBR)

plot(dispersionBR,hull=F,ellipse=T,lwd=1)


set.seed(500)
pairwise.adonis(distBR,factors=dataBR2$season,p.adjust.m="bonferroni", perm=100000)


### CAP Bruce
###
system.time(
  OM1BR<-CAPdiscrim(dataBR1log~season,data=dataBR2,dist="bray",axes=3,m=0,add=F, permutations = 1000)
)
beep()
summary(OM1BR)
OM1BR$manova

#extract CAP scores (x and y coordinates)
data.scoresCAP_BR = as.data.table(scores(OM1BR))

#add columns to data frame 
data.scoresCAP_BR$season = dataBR2$season
head(data.scoresCAP_BR)


###fit sps data in Adonga
fitCAP_BR<-envfit(OM1BR,dataBR1log,permutations=999)
arrowCAP_BR<-data.frame(fitCAP_BR$vectors$arrows,R=fitCAP_BR$vectors$r,P=fitCAP_BR$vectors$pvals)
arrowCAP_BR$FG <- rownames(arrowCAP_BR)
arrowarrowCAP_BR.p<-arrowCAP_BR[arrowCAP_BR$P<=0.05&arrowCAP_BR$R>=0.15,]

PP_CAP_BR<-ggplot(data.scoresCAP_BR,aes(x=LD1,y=LD2))+
  geom_point(size=4,aes(colour=season))+
  #geom_text(aes(label=month),size=7)+
  labs(x="LD1",y="LD2",colour="Season",fill="Season")+
  stat_ellipse(size=2,type="t",aes(group=season,colour=season,fill=season),level=.6,geom="polygon",alpha=.2)+
  #scale_colour_manual(values=p)+
  ggtitle("LD in Bruce, with log, with dummy")+
  theme_bw()
PP_CAP_BR+
  geom_segment(data=arrowarrowCAP_BR.p,aes(x=0,y=0,xend=LD1*R*6,yend=LD2*R*6),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowarrowCAP_BR.p,aes(x=LD1*R*6,y=LD2*R*6,label=FG),cex=5,direction="both",segment.size=0.25)




######Escadinhas

dataE<-data[site=="E"]
dataE1<-dataE[,!c(1:3,89)]
dataE1[,dummy:=127.3885]

dataE2<-dataE[,c(1:3,89)]

dataE1log<-log(dataE1+1)


#### run NMDS
set.seed(1)
NMDSE<-metaMDS(dataE1log,distance="bray",k=3,trymax=1000,autotransform = F)
beep()

#stressplot(NMDSE)

#extract NMDS scores (x and y coordinates)
data.scoresE=as.data.table(scores(NMDSE))

#add columns to data frame 
data.scoresE$coreID = dataE2$coreID
data.scoresE$site = dataE2$site
data.scoresE$season = factor(dataE2$season)

head(data.scoresE)

###fit sps data in Escadinhas
fitE<-envfit(NMDSE,dataE1log,permutations=999)
arrowE<-data.frame(fitE$vectors$arrows,R = fitE$vectors$r, P = fitE$vectors$pvals)
arrowE$FG <- rownames(arrowE)
arrowE.p<-arrowE[arrowE$P<=0.05&arrowE$R>=.2,]

PPE<-ggplot(data.scoresE,aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=season))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Season",fill="Season")+
  stat_ellipse(size=2, type="t",aes(group=season, colour=season,fill=season),level=.6,geom="polygon",alpha=.15)+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Escadinhas, with log, with dummy")+
  theme_bw()
PPE+
  geom_segment(data=arrowE.p,aes(x=0,y=0,xend=NMDS1*R*2.5,yend=NMDS2*R*2.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowE.p,aes(x=NMDS1*R*2.5,y=NMDS2*R*2.5,label=FG),cex=5,direction="both",segment.size=0.25)

###Permanova
set.seed(601)
distE<-vegdist(dataE1log,method="bray")
perE<-adonis(distE~season,data=dataE2,permutations=999)
perE

set.seed(602)
dispersionE<-betadisper(distE,group=dataE2$season)
permutest(dispersionE)

TukeyHSD(dispersionE)

plot(dispersionE,hull=F,ellipse=T,lwd=1)

set.seed(600)
pairwise.adonis(distE,factors=dataE2$season,p.adjust.m="bonferroni", perm=100000)


### CAP Escadinhas
###
system.time(
  OM1E<-CAPdiscrim(dataE1log~season,data=dataE2,dist="bray",axes=3,m=0,add=F, permutations = 1000)
)
beep()
summary(OM1E)
OM1E$manova

#extract CAP scores (x and y coordinates)
data.scoresCAP_E = as.data.table(scores(OM1E))

#add columns to data frame 
data.scoresCAP_E$season = dataE2$season
head(data.scoresCAP_E)


###fit sps data in Adonga
fitCAP_E<-envfit(OM1E,dataE1log,permutations=999)
arrowCAP_E<-data.frame(fitCAP_E$vectors$arrows,R=fitCAP_E$vectors$r,P=fitCAP_E$vectors$pvals)
arrowCAP_E$FG <- rownames(arrowCAP_E)
arrowarrowCAP_E.p<-arrowCAP_E[arrowCAP_E$P<=0.05&arrowCAP_E$R>=0.15,]

PP_CAP_E<-ggplot(data.scoresCAP_E,aes(x=LD1,y=LD2))+
  geom_point(size=4,aes(colour=season))+
  #geom_text(aes(label=month),size=7)+
  labs(x="LD1",y="LD2",colour="Season",fill="Season")+
  stat_ellipse(size=2,type="t",aes(group=season,colour=season,fill=season),level=.6,geom="polygon",alpha=.2)+
  #scale_colour_manual(values=p)+
  ggtitle("LD in Escadinhas, with log, with dummy")+
  theme_bw()
PP_CAP_E+
  geom_segment(data=arrowarrowCAP_E.p,aes(x=0,y=0,xend=LD1*R*6,yend=LD2*R*6),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowarrowCAP_E.p,aes(x=LD1*R*6,y=LD2*R*6,label=FG),cex=5,direction="both",segment.size=0.25)




######Adonga

dataAD<-data[site=="AD"]
dataAD1<-dataAD[,!c(1:3,89)]
dataAD1[,dummy:=88.49558]

dataAD2<-dataAD[,c(1:3,89)]

dataAD1log<-log(dataAD1+1)


#### run NMDS
set.seed(1)
NMDSAD<-metaMDS(dataAD1log,distance="bray",k=3,trymax=1000,autotransform = F)
beep()

#stressplot(NMDSAD)

#extract NMDS scores (x and y coordinates)
data.scoresAD=as.data.table(scores(NMDSAD))

#add columns to data frame 
data.scoresAD$coreID = dataAD2$coreID
data.scoresAD$site = dataAD2$site
data.scoresAD$season = factor(dataAD2$season)

head(data.scoresAD)

###fit sps data in Adonga
fitAD<-envfit(NMDSAD,dataAD1log,permutations=999)
arrowAD<-data.frame(fitAD$vectors$arrows,R = fitAD$vectors$r, P = fitAD$vectors$pvals)
arrowAD$FG <- rownames(arrowAD)
arrowAD.p<-arrowAD[arrowAD$P<=0.05&arrowAD$R>=.2,]

PPAD<-ggplot(data.scoresAD,aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=season))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Season",fill="Season")+
  stat_ellipse(size=2, type="t",aes(group=season, colour=season,fill=season),level=.6,geom="polygon",alpha=.15)+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Adonga, with log, with dummy")+
  theme_bw()
PPAD+
  geom_segment(data=arrowAD.p,aes(x=0,y=0,xend=NMDS1*R*2.5,yend=NMDS2*R*2.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowAD.p,aes(x=NMDS1*R*2.5,y=NMDS2*R*2.5,label=FG),cex=5,direction="both",segment.size=0.25)

###Permanova

distAD<-vegdist(dataAD1log,method="bray")
perAD<-adonis(distAD~season,data=dataAD2,permutations=999)
perAD

dispersionAD<-betadisper(distAD,group=dataAD2$season)
permutest(dispersionAD)

plot(dispersionAD,hull=F,ellipse=T,lwd=1)

set.seed(700)
pairwise.adonis(distAD,factors=dataAD2$season,p.adjust.m="bonferroni", perm=100000)


### CAP Adonga
###
system.time(
  OM1AD<-CAPdiscrim(dataAD1log~season,data=dataAD2,dist="bray",axes=3,m=0,add=F, permutations = 1000)
)
beep()
summary(OM1AD)
OM1AD$manova

#extract CAP scores (x and y coordinates)
data.scoresCAP_AD = as.data.table(scores(OM1AD))

#add columns to data frame 
data.scoresCAP_AD$season = dataAD2$season
head(data.scoresCAP_AD)


###fit sps data in Adonga
fitCAP_AD<-envfit(OM1AD,dataAD1log,permutations=999)
arrowCAP_AD<-data.frame(fitCAP_AD$vectors$arrows,R=fitCAP_AD$vectors$r,P=fitCAP_AD$vectors$pvals)
arrowCAP_AD$FG <- rownames(arrowCAP_AD)
arrowarrowCAP_AD.p<-arrowCAP_AD[arrowCAP_AD$P<=0.05&arrowCAP_AD$R>=0.15,]

PP_CAP_AD<-ggplot(data.scoresCAP_AD,aes(x=LD1,y=LD2))+
  geom_point(size=4,aes(colour=season))+
  #geom_text(aes(label=month),size=7)+
  labs(x="LD1",y="LD2",colour="Season",fill="Season")+
  stat_ellipse(size=2,type="t",aes(group=season,colour=season,fill=season),level=.6,geom="polygon",alpha=.2)+
  #scale_colour_manual(values=p)+
  ggtitle("LD in Adonga, with log, with dummy")+
  theme_bw()
PP_CAP_AD+
  geom_segment(data=arrowarrowCAP_AD.p,aes(x=0,y=0,xend=LD1*R*6,yend=LD2*R*6),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowarrowCAP_AD.p,aes(x=LD1*R*6,y=LD2*R*6,label=FG),cex=5,direction="both",segment.size=0.25)

############################################
###############################################################
##############################################################################
###############################################################################################

##### Now by season
p<-c("#4DAF4A","darkgreen","#377EB8","#984EA3","#FF7F00","#E41A1C")
####################
######Begining (Oct+Nov)

dataBeg<-data[season=="begining"]
dataBeg[,site:=factor(site,levels=c("A","AB","BI","E","BR","AD"))]
dataBeg1<-dataBeg[,!c(1:3,86)]
dataBeg1[,dummy:=88.49558]

dataBeg2<-dataBeg[,c(1:3,86)]

dataBeg1log<-log(dataBeg1+1)


#### run NMDS
set.seed(100)
NMDSBeg<-metaMDS(dataBeg1log,distance="bray",k=3,trymax=1000,autotransform = F)
beep()

#stressplot(NMDSBeg)

#extract NMDS scores (x and y coordinates)
data.scoresBeg=as.data.table(scores(NMDSBeg))

#add columns to data frame 
data.scoresBeg$coreID = dataBeg2$coreID
data.scoresBeg$site = dataBeg2$site
data.scoresBeg$season = factor(dataBeg2$season)

head(data.scoresBeg)

###fit sps data in Begining
set.seed(101)
fitBeg<-envfit(NMDSBeg,dataBeg1log,permutations=999)
arrowBeg<-data.frame(fitBeg$vectors$arrows,R = fitBeg$vectors$r, P = fitBeg$vectors$pvals)
arrowBeg$FG <- rownames(arrowBeg)
arrowBeg.p<-arrowBeg[arrowBeg$P<=0.05,]

PPBeg<-ggplot(data.scoresBeg,aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Site",fill="Site")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site,fill=site),level=.6,geom="polygon",alpha=.15)+
  scale_colour_manual(values=p)+
  scale_fill_manual(values=p)+
  ggtitle("NMDS in Begining, with log, with dummy")+
  theme_bw()
PPBeg+
  geom_segment(data=arrowBeg.p,aes(x=0,y=0,xend=NMDS1*R*2.5,yend=NMDS2*R*2.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowBeg.p,aes(x=NMDS1*R*2.5,y=NMDS2*R*2.5,label=FG),cex=5,direction="both",segment.size=0.25)

###Permanova
set.seed(103)
distBeg<-vegdist(dataBeg1log,method="bray")
perBeg<-adonis(distBeg~site,data=dataBeg2,permutations=999)

set.seed(104)
dispersionBeg<-betadisper(distBeg,group=dataBeg2$site)
permutest(dispersionBeg)

TukeyHSD(dispersionBeg)

plot(dispersionBeg,hull=F,ellipse=T,lwd=1)

##pairwise
set.seed(200)
pairwise.adonis(distBeg,factors=dataBeg2$site,p.adjust.m="bonferroni", perm=100000)

###test of pairwise validity
dataBIA<-dataBeg[site=="A"|site=="BI"]
table(dataBIA$site)
table(dataBIA$season)

distBIA<-vegdist(log(dataBIA[,-c(1:3,89)]+1),method="bray")
perBIA<-adonis(distBIA~site,data=dataBIA[,c(1:3,89)],permutations=999)
perBIA



### CAP Begining
###
set.seed(105)
system.time(
  OM1Beg<-CAPdiscrim(dataBeg1log~site,data=dataBeg2,dist="bray",axes=3,m=0,add=F, permutations = 1000)
)
beep(3)
summary(OM1Beg)
OM1Beg$manova
OM1Beg$m

#extract CAP scores (x and y coordinates)
data.scoresCAP_Beg = as.data.table(scores(OM1Beg))

#add columns to data frame 
data.scoresCAP_Beg$site = dataBeg2$site
head(data.scoresCAP_Beg)


###fit sps data in
set.seed(106)
fitCAP_Beg<-envfit(OM1Beg,dataBeg1log,permutations=10000)
arrowCAP_Beg<-data.frame(fitCAP_Beg$vectors$arrows,R=fitCAP_Beg$vectors$r,P=fitCAP_Beg$vectors$pvals)
arrowCAP_Beg$FG <- rownames(arrowCAP_Beg)
arrowCAP_Beg.p<-arrowCAP_Beg[arrowCAP_Beg$P<=0.05&arrowCAP_Beg$R>=0.15,]
#arrowCAP_Beg.p<-arrowCAP_Beg[arrowCAP_Beg$P<=0.05,]

PP_CAP_Beg<-ggplot(data.scoresCAP_Beg,aes(x=LD1,y=LD2))+
  geom_point(size=4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x="LD1",y="LD2",colour="Site",fill="Site")+
  stat_ellipse(size=2,type="t",aes(group=site,colour=site,fill=site),level=.6,geom="polygon",alpha=.2)+
  scale_colour_manual(values=p)+
  scale_fill_manual(values=p)+
  ggtitle("LD in Begining, with log, with dummy")+
  theme_bw()
PP_CAP_Beg+
  geom_segment(data=arrowCAP_Beg.p,aes(x=0,y=0,xend=LD1*R*6,yend=LD2*R*6),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowCAP_Beg.p,aes(x=LD1*R*6,y=LD2*R*6,label=FG),cex=5,direction="both",segment.size=0.25)





######Mid(Dec+Jan+Feb)

dataMid<-data[season=="mid"]
dataMid[,site:=factor(site,levels=c("A","AB","BI","E","BR","AD"))]
dataMid1<-dataMid[,!c(1:3,86)]
dataMid1[,dummy:=88.49558]

dataMid2<-dataMid[,c(1:3,86)]

dataMid1log<-log(dataMid1+1)


#### run NMDS
set.seed(107)
NMDSMid<-metaMDS(dataMid1log,distance="bray",k=3,trymax=1000,autotransform = F)
beep()

#stressplot(NMDSMid)

#extract NMDS scores (x and y coordinates)
data.scoresMid=as.data.table(scores(NMDSMid))

#add columns to data frame 
data.scoresMid$coreID = dataMid2$coreID
data.scoresMid$site = dataMid2$site
data.scoresMid$season = factor(dataMid2$season)

head(data.scoresMid)

###fit sps data in Mid
set.seed(108)
fitMid<-envfit(NMDSMid,dataMid1log,permutations=999)
arrowMid<-data.frame(fitMid$vectors$arrows,R = fitMid$vectors$r, P = fitMid$vectors$pvals)
arrowMid$FG <- rownames(arrowMid)
arrowMid.p<-arrowMid[arrowMid$P<=0.05&arrowMid$R>=0.15,]

PPMid<-ggplot(data.scoresMid,aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Site",fill="Site")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site,fill=site),level=.6,geom="polygon",alpha=.15)+
  scale_colour_manual(values=p)+
  scale_fill_manual(values=p)+
  ggtitle("NMDS in Mid(Dec+Jan+Feb), with log, with dummy")+
  theme_bw()
PPMid+
  geom_segment(data=arrowMid.p,aes(x=0,y=0,xend=NMDS1*R*2.5,yend=NMDS2*R*2.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowMid.p,aes(x=NMDS1*R*2.5,y=NMDS2*R*2.5,label=FG),cex=5,direction="both",segment.size=0.25)

###Permanova
set.seed(109)
distMid<-vegdist(dataMid1log,method="bray")
perMid<-adonis(distMid~site,data=dataMid2,permutations=999)

set.seed(110)
dispersionMid<-betadisper(distMid,group=dataMid2$site)
permutest(dispersionMid)

TukeyHSD(dispersionMid)

plot(dispersionMid,hull=F,ellipse=T,lwd=1)

##pairwise
set.seed(300)
pairwise.adonis(distMid,factors=dataMid2$site,p.adjust.m="bonferroni", perm=100000)


### CAP Mid
###
set.seed(11111111)

OM1Mid<-CAPdiscrim(dataMid1log~site,data=dataMid2,dist="bray",axes=3,m=0,add=F, permutations = 1000)
beep(3)

summary(OM1Mid)
OM1Mid$manova
OM1Mid$m

#extract CAP scores (x and y coordinates)
data.scoresCAP_Mid = as.data.table(scores(OM1Mid))

#add columns to data frame 
data.scoresCAP_Mid$site = dataMid2$site
head(data.scoresCAP_Mid)


###fit sps data in Mid
set.seed(11211111)
fitCAP_Mid<-envfit(OM1Mid,dataMid1log,permutations=10000)
arrowCAP_Mid<-data.frame(fitCAP_Mid$vectors$arrows,R=fitCAP_Mid$vectors$r,P=fitCAP_Mid$vectors$pvals)
arrowCAP_Mid$FG <- rownames(arrowCAP_Mid)
arrowarrowCAP_Mid.p<-arrowCAP_Mid[arrowCAP_Mid$P<=0.05&arrowCAP_Mid$R>=0.15,]
#arrowarrowCAP_Mid.p<-arrowCAP_Mid[arrowCAP_Mid$P<=0.05,]

PP_CAP_Mid<-ggplot(data.scoresCAP_Mid,aes(x=LD1,y=LD2))+
  geom_point(size=4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x="LD1",y="LD2",colour="Site",fill="Site")+
  stat_ellipse(size=2,type="t",aes(group=site,colour=site,fill=site),level=.6,geom="polygon",alpha=.2)+
  scale_colour_manual(values=p)+
  scale_fill_manual(values=p)+
  ggtitle("LD in Mid(Dec+Jan+Feb), with log, with dummy")+
  theme_bw()
PP_CAP_Mid+
  geom_segment(data=arrowarrowCAP_Mid.p,aes(x=0,y=0,xend=LD1*R*6,yend=LD2*R*6),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowarrowCAP_Mid.p,aes(x=LD1*R*6,y=LD2*R*6,label=FG),cex=5,direction="both",segment.size=0.25)





######End (Mar+Apr)

dataEnd<-data[season=="end"]
dataEnd[,site:=factor(site,levels=c("A","AB","BI","E","BR","AD"))]
dataEnd1<-dataEnd[,!c(1:3,86)]
dataEnd1[,dummy:=88.49558]

dataEnd2<-dataEnd[,c(1:3,86)]

dataEnd1log<-log(dataEnd1+1)


#### run NMDS
set.seed(113)
NMDSEnd<-metaMDS(dataEnd1log,distance="bray",k=3,trymax=1000,autotransform = F)
beep()

#stressplot(NMDSEnd)

#extract NMDS scores (x and y coordinates)
data.scoresEnd=as.data.table(scores(NMDSEnd))

#add columns to data frame 
data.scoresEnd$site = dataEnd2$site
head(data.scoresEnd)

###fit sps data in End
set.seed(114)
fitEnd<-envfit(NMDSEnd,dataEnd1log,permutations=999)
arrowEnd<-data.frame(fitEnd$vectors$arrows,R = fitEnd$vectors$r, P = fitEnd$vectors$pvals)
arrowEnd$FG <- rownames(arrowEnd)
arrowEnd.p<-arrowEnd[arrowEnd$P<=0.05,]

PPEnd<-ggplot(data.scoresEnd,aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Site",fill="Site")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site,fill=site),level=.6,geom="polygon",alpha=.15)+
  scale_colour_manual(values=p)+
  scale_fill_manual(values=p)+
  ggtitle("NMDS in End(Dec+Jan+Feb), with log, with dummy")+
  theme_bw()
PPEnd+
  geom_segment(data=arrowEnd.p,aes(x=0,y=0,xend=NMDS1*R*2.5,yend=NMDS2*R*2.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowEnd.p,aes(x=NMDS1*R*2.5,y=NMDS2*R*2.5,label=FG),cex=5,direction="both",segment.size=0.25)

###Permanova
set.seed(115)
distEnd<-vegdist(dataEnd1log,method="bray")
perEnd<-adonis(distEnd~site,data=dataEnd2,permutations=999)

set.seed(116)
dispersionEnd<-betadisper(distEnd,group=dataEnd2$site)
permutest(dispersionEnd)

TukeyHSD(dispersionEnd)

plot(dispersionEnd,hull=F,ellipse=T,lwd=1)


##pairwise
set.seed(300)
pairwise.adonis(distEnd,factors=dataEnd2$site,p.adjust.m="bonferroni", perm=100000)


### CAP End
###

test1<-dataEnd1log[-which(dataEnd2$site=="AD")]
test2<-dataEnd2[-which(dataEnd2$site=="AD")]
set.seed(117)
system.time(
  OM1End<-CAPdiscrim(dataEnd1log~site,data=dataEnd2,dist="bray",axes=3,m=0,add=F, permutations = 1000)
)
beep(4)
summary(OM1End)
OM1End$manova
OM1End$m

#extract CAP scores (x and y coordinates)
data.scoresCAP_End = as.data.table(scores(OM1End))

#add columns to data frame 
data.scoresCAP_End$site = dataEnd2$site
head(data.scoresCAP_End)


###fit sps data in End
set.seed(118)
fitCAP_End<-envfit(OM1End,dataEnd1log,permutations=10000)
arrowCAP_End<-data.frame(fitCAP_End$vectors$arrows,R=fitCAP_End$vectors$r,P=fitCAP_End$vectors$pvals)
arrowCAP_End$FG <- rownames(arrowCAP_End)
arrowarrowCAP_End.p<-arrowCAP_End[arrowCAP_End$P<=0.05&arrowCAP_End$R>=0.15,]
#arrowarrowCAP_End.p<-arrowCAP_End[arrowCAP_End$P<=0.05,]

PP_CAP_End<-ggplot(data.scoresCAP_End,aes(x=LD1,y=LD2))+
  geom_point(size=4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x="LD1",y="LD2",colour="Site",fill="Site")+
  stat_ellipse(size=2,type="t",aes(group=site,colour=site,fill=site),level=.6,geom="polygon",alpha=.2)+
  scale_colour_manual(values=p)+
  scale_fill_manual(values=p)+
  ggtitle("LD in End, with log, with dummy")+
  theme_bw()
PP_CAP_End+
  geom_segment(data=arrowarrowCAP_End.p,aes(x=0,y=0,xend=LD1*R*6,yend=LD2*R*6),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowarrowCAP_End.p,aes(x=LD1*R*6,y=LD2*R*6,label=FG),cex=5,direction="both",segment.size=0.25)






