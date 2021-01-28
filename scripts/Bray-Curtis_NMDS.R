setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("vegan","ggplot2","viridis","RColorBrewer","psych","reshape2","beepr","data.table")
lapply(packs,require,character.only=T)

## Load database 
DB66<-fread("data_out/db/Final_DB_lowtaxa_density_polyexcl_20210121.csv") ### created in script called Database_cleanup_joining
str(DB66)


##Remove non-target benthos
unique(DB66$low_taxa[which(DB66$class1=="Other")])
db<-DB66[-which(DB66$class1=="Other"),]
unique(db$low_taxa)

##Check data
sum(db[low_taxa=="Polychaeta_errantia",numb])/sum(db$numb)*100
sum(db[low_taxa=="Polychaeta_sedentaria",numb])/sum(db$numb)*100

xx<-db[low_taxa=="Capitella_sp"&numb!=0]
db[family=="Lasaeidae"&numb!=0]

## Reduce the dimension of taxa names 
db[,taxaf:=low_taxa]

# Juntar todos os Capitellidae num so (Capitellidae+Capitella_sp+Heteromastus_filiformis,Notomastus_fauveli)
db[taxaf=="Capitella_sp"|taxaf=="Heteromastus_filiformis"|taxaf=="Notomastus_fauveli",taxaf:="Capitellidae"]

# Juntar todos os Cirratulidae (Cirratulidae+Cirriformia_sp+Kirkegaardia_sp)
db[taxaf=="Cirratulidae"|taxaf=="Cirriformia_sp"|taxaf=="Kirkegaardia_sp",taxaf:="Cirratulidae"]

# Juntar marphysas de um lado (incluindo Marphysa_sp+Marphysa_sanguinea) e o resto dos eunicidae de outro (Eunice_sp+Eunicidae)
db[taxaf=="Marphysa_sp",taxaf:="Marphysa_sanguinea"]
db[taxaf=="Eunice_sp",taxaf:="Eunicidae"]

# Juntar MaldanidaeA e Petaloproctus_sp aos restantes Maldanidae
db[taxaf=="MaldanidaeA"|taxaf=="Petaloproctus_sp",taxaf:="Maldanidae"]

# Pilargidae to be changed to Sigambra_sp
db[taxaf=="Pilargidae",taxaf:="Sigambra_sp"]

# Polycirrus_sp+Streblosoma_sp to be joined to the rest of Terebellidae
db[taxaf=="Polycirrus_sp"|taxaf=="Streblosoma_sp",taxaf:="Terebellidae"]

# Nereis2 para ser joined a Nereididae
db[taxaf=="Nereis2",taxaf:="Nereididae"]

# Juntar aos Paraonidae: Aricidea_sp+Aricidea_spA+Aricidea_spB+ParaonidaeA
db[taxaf=="Aricidea_sp"|taxaf=="Aricidea_spA"|taxaf=="Aricidea_spB"|taxaf=="ParaonidaeA",taxaf:="Paraonidae"]

# Orbiniidae sera a juncao:OrbinidaeA+Leodamas_sp+Orbiniidae
db[taxaf=="OrbinidaeA"|taxaf=="Leodamas_sp",taxaf:="Orbiniidae"]

# Turbonilla_sp has to be eliminated from database, appears only with zeros
# Eliminar Megalopas da ase de dados para sp richness e bray curtis
# Remover Pachygraspus_gracilis de toda a base de dados (está tudo a zero))
db1<-db[!taxaf=="Turbonilla_sp"][!taxaf=="Megalopa"][!taxaf=="Pachygraspus_gracilis"]

db1[,unique(taxaf)]



###remove data from 2020
db2<-db1[!year==2020]
db2[year==2020]
db2[,unique(taxaf)]

###aggregate and reshape database for analysis
db3<-db2[,lapply(.SD,sum,na.rm=T),.SDcols="numb",by=c("site","month","coreID","taxaf")]
DB<-dcast.data.table(db3,coreID+site+month~taxaf,value.var="numb")
setkey(DB,coreID,site,month) ## isto define as variaveis core site e month como as variaveis de base para qualquer operação

###Calculate densities
dens1<-function(x){x/0.0113} #for Adonga
dens2<-function(x){x/0.00785} # for the rest of the sites
DB1<-DB[,lapply(.SD,ifelse(site=="AD",dens1,dens2)),by=c("coreID","site","month")]
setkey(DB1,coreID,month,site)

write.csv(DB1,"Data_out/db/DB_multianal_20210127.csv",row.names=F)
###NMDS with Bray curtis
##Prepare data and remove rows with zeros in all columns

data<-DB1[apply(DB1[,!(1:3)],1,sum)!=0] ### NMDS and bray curtis do not work with empty cores. So we have to remove all of them

data1<-data[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

data2<-data[,1:3] ###store aggregating variables to use latter
head(data2)

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



##############NMDS and bray curtis

table(is.na(data1)) ###check if there is any NAs

###subsetting for reduced sps numb

####### 75%
nam<-m2[p75=="N",sp]
str(nam)
dataS<-data[,!..nam]
table(is.na(dataS))


####### 95%
nam1<-m2[p95=="N",sp]
str(nam1)
dataS1<-data[,!..nam1]
table(is.na(dataS1))

####### sp totalling at least 1% of total abundance
nam2<-m2[p1p=="N",sp]
str(nam2)
dataS2<-data[,!..nam2]
table(is.na(dataS2))



##Remove cores with zero species after selection of p75
dataSS<-dataS[apply(dataS[,!1:3],1,sum)!=0]
dataSS1<-dataS1[apply(dataS1[,!1:3],1,sum)!=0]
dataSS2<-dataS2[apply(dataS2[,!1:3],1,sum)!=0]

dataSS$dummy<-88.49558
dataSS1$dummy<-88.49558
dataSS2$dummy<-88.49558

dataSS_1<-dataSS[,!1:3]
dataSS_2<-dataSS[,1:3]

dataSS1_1<-dataSS1[,!1:3]
dataSS1_2<-dataSS1[,1:3]

dataSS2_1<-dataSS2[,!1:3]
dataSS2_2<-dataSS2[,1:3]


ff<-function(x){log(x+1)}
dataSS_1log<-dataSS_1[,lapply(.SD,ff)]
dataSS1_1log<-dataSS1_1[,lapply(.SD,ff)]
dataSS2_1log<-dataSS2_1[,lapply(.SD,ff)]

####with dummy without log
data1$dummy<-88.49558
data1log<-log(data1+1)

table(data$month,data$site)

table(data$site)


set.seed(0)
NMDS<-metaMDS(data1log,distance="bray",k=2,trymax=1000,autotransform = F)
beep()

set.seed(1)
NMDS_d75<-metaMDS(log(dataSS_1+1),distance="bray",k=3,trymax=500,autotransform = F)
beep()

set.seed(2)
NMDS_d95<-metaMDS(dataSS1_1,distance="bray",k=2,trymax=1000,autotransform = F)
beep()

set.seed(3)
NMDS_dp1p<-metaMDS(dataSS2_1,distance="bray",k=2,trymax=1000,autotransform = F)
beep()

####with dummy with log
set.seed(4)
NMDS_d75log<-metaMDS(dataSS_1log,distance="bray",k=2,trymax=1000,autotransform = F)
beep()

set.seed(5)
NMDS_d95log<-metaMDS(dataSS1_1log,distance="bray",k=2,trymax=1000,autotransform = F)
beep()

set.seed(6)
NMDS_dp1plog<-metaMDS(dataSS2_1log,distance="bray",k=2,trymax=1000,autotransform = F)
beep()


###ploting
plot(NMDS)
stressplot(NMDS_d75)


#extract NMDS scores (x and y coordinates)
data.scores = as.data.table(scores(NMDS))

#add columns to data frame 
data.scores$coreID = data2$coreID
data.scores$site = data2$site
data.scores$month = data2$month

head(data.scores)

p<-c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","black")

ggplot(data.scores, aes(x = NMDS1, y = NMDS2, colour=site,group=site)) + 
  #geom_point(size = 4)+
  geom_text(aes(label=month),size=7)+
  #stat_summary(geom="pointrange",size = 1, aes(colour=factor(month),shape=site))+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", y = "NMDS2",colour="Site")+
  stat_ellipse(size=2, type="t")+
  scale_colour_manual(values=p)+
  ggtitle("NMDS all sp, with log, with dummy")+
  theme_bw()


###fit sps data
fit<-envfit(NMDS,log(data1+1), permutations=999)
arrow<-data.frame(fit$vectors$arrows,R = fit$vectors$r, P = fit$vectors$pvals)
arrow$FG <- rownames(arrow)
arrow.p<-arrow[arrow$P<=0.05,]


ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Site")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site))+
  scale_colour_manual(values=p)+
  ggtitle("NMDS all sp, with log, with dummy")+
  geom_segment(data=arrow.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  geom_text(data=arrow.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),size=5)+
  theme_bw()

###fit site and month data
dataSS1_2[,month1:=factor(month)]

fit1<-envfit(NMDS_d95,dataSS1_2[,-c(1,3)], permutations=999)
cent<-as.data.frame(scores(fit1, display = "factors"))
cent$var<-row.names(cent)
cent$R<-c(rep(fit1$factors$r[1],6),rep(fit1$factors$r[2],7))
centmonth<-cent[-c(1:6),]

ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 2,aes(colour=site))+
  #geom_text(aes(label=month,colour=site),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Site")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site))+
  scale_colour_manual(values=p)+
  ggtitle("NMDS d95, without log, with dummy")+
  #geom_segment(data=arrow.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  #geom_text(data=arrow.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),size=5)+
  #ggrepel::geom_text_repel(data=arrow.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)+
  geom_text(data=centmonth,aes(x=NMDS1,y=NMDS2,label=var),size=5.5)+
  theme_bw()


################separating sites #########################################################################
#####Adonga

###NMDS with Bray curtis
##Prepare data and remove rows with zeros in all columns

dataAD<-data[data$site=="AD",] ###filtrar para adonga

dataAD1<-dataAD[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

dataAD2<-dataAD[,1:3] ###store aggregating variables to use latter
head(dataAD2)


##############NMDS and bray curtis

table(is.na(dataAD1)) ###check if there is any NAs

###subsetting for reduced sps numb

####### 75%
#nam<-m2[p75=="N",sp]
#str(nam)
dataADS<-dataAD[,!..nam]
table(is.na(dataADS))


####### 95%
#nam1<-m2[p95=="N",sp]
#str(nam1)
dataADS1<-dataAD[,!..nam1]
table(is.na(dataADS1))

####### sp totalling at least 1% of total abundance
#nam2<-m2[p1p=="N",sp]
#str(nam2)
dataADS2<-dataAD[,!..nam2]
table(is.na(dataADS2))



##Remove cores with zero species after selection of p75, P95, P1P
dataADSS<-dataADS[apply(dataADS[,!1:3],1,sum)!=0]
dataADSS1<-dataADS1[apply(dataADS1[,!1:3],1,sum)!=0]
dataADSS2<-dataADS2[apply(dataADS2[,!1:3],1,sum)!=0]

dataAD1$dummy<-88.49558
dataADSS$dummy<-88.49558
dataADSS1$dummy<-88.49558
dataADSS2$dummy<-88.49558

dataADSS_1<-dataADSS[,!1:3]
dataADSS_2<-dataADSS[,1:3]

dataADSS1_1<-dataADSS1[,!1:3]
dataADSS1_2<-dataADSS1[,1:3]

dataADSS2_1<-dataADSS2[,!1:3]
dataADSS2_2<-dataADSS2[,1:3]


####with dummy without log

set.seed(0)
NMDSAD<-metaMDS(log(dataAD1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSAD)
stressplot(NMDSAD)

#extract NMDS scores (x and y coordinates)
data.scoresAD = as.data.table(scores(NMDSAD))

#add columns to data frame 
data.scoresAD$coreID = dataAD2$coreID
data.scoresAD$site = dataAD2$site
data.scoresAD$month = factor(dataAD2$month)

head(data.scoresAD)

#p<-c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","black")

PP_AD<-ggplot(data.scoresAD, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Adonga, without log, with dummy")+
  theme_bw()

###fit sps data in Adonga
fitAD<-envfit(NMDSAD,log(dataAD1+1),permutations=999)
arrowAD<-data.frame(fitAD$vectors$arrows,R = fitAD$vectors$r, P = fitAD$vectors$pvals)
arrowAD$FG <- rownames(arrowAD)
arrowAD.p<-arrowAD[arrowAD$P<=0.05,]

PP_AD+
  geom_segment(data=arrowAD.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowAD.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


set.seed(10)
NMDSAD_d75<-metaMDS(log(dataADSS_1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSAD_d75)
stressplot(NMDSAD_d75)

#extract NMDS scores (x and y coordinates)
data.scoresAD75 = as.data.table(scores(NMDSAD_d75))

#add columns to data frame 
data.scoresAD75$coreID = dataADSS_2$coreID
data.scoresAD75$site = dataADSS_2$site
data.scoresAD75$month = factor(dataADSS_2$month)

head(data.scoresAD75)

PP_AD75<-ggplot(data.scoresAD75, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p75 in Adonga, with log, with dummy")+
  theme_bw()

###fit sps data in Adonga
fitAD75<-envfit(NMDSAD_d75,log(dataADSS_1+1),permutations=999)
arrowAD75<-data.frame(fitAD75$vectors$arrows,R = fitAD75$vectors$r, P = fitAD75$vectors$pvals)
arrowAD75$FG <- rownames(arrowAD75)
arrowAD75.p<-arrowAD75[arrowAD75$P<=0.05,]

PP+geom_segment(data=arrowAD75.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowAD75.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)




set.seed(11)
NMDSAD_d95<-metaMDS(dataADSS1_1,distance="bray",k=2,trymax=1000, autotransform = F)
beep()

plot(NMDSAD_d95)
stressplot(NMDSAD_d95)

#extract NMDS scores (x and y coordinates)
data.scoresAD95 = as.data.table(scores(NMDSAD_d95))

#add columns to data frame 
data.scoresAD95$coreID = dataADSS1_2$coreID
data.scoresAD95$site = dataADSS1_2$site
data.scoresAD95$month = factor(dataADSS1_2$month)

head(data.scoresAD95)

PP_AD95<-ggplot(data.scoresAD95, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p95 in Adonga, without log, with dummy")+
  theme_bw()

###fit sps data in Adonga
fitAD95<-envfit(NMDSAD_d95,dataADSS1_1,permutations=999)
arrowAD95<-data.frame(fitAD95$vectors$arrows,R = fitAD95$vectors$r, P = fitAD95$vectors$pvals)
arrowAD95$FG <- rownames(arrowAD95)
arrowAD95.p<-arrowAD95[arrowAD95$P<=0.05,]

PP_AD95+geom_segment(data=arrowAD95.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowAD95.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


set.seed(12)
NMDSAD_dp1p<-metaMDS(dataADSS2_1,distance="bray",k=2,trymax=1000,autotransform = F)
beep()


plot(NMDSAD_dp1p)
stressplot(NMDSAD_dp1p)

#extract NMDS scores (x and y coordinates)
data.scoresdp1p = as.data.table(scores(NMDSAD_dp1p))

#add columns to data frame 
data.scoresdp1p$coreID = dataADSS2_2$coreID
data.scoresdp1p$site = dataADSS2_2$site
data.scoresdp1p$month = factor(dataADSS2_2$month)

head(data.scoresdp1p)

PP_ADdp1p<-ggplot(data.scoresdp1p, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p1p in Adonga, without log, with dummy")+
  theme_bw()

###fit sps data in Adonga
fitdp1p<-envfit(NMDSAD_dp1p,dataADSS2_1,permutations=999)
arrowdp1p<-data.frame(fitdp1p$vectors$arrows,R = fitdp1p$vectors$r, P = fitdp1p$vectors$pvals)
arrowdp1p$FG <- rownames(arrowdp1p)
arrowdp1p.p<-arrowdp1p[arrowdp1p$P<=0.05,]

PP_ADdp1p+geom_segment(data=arrowdp1p.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowdp1p.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


#####Anrumai

###NMDS with Bray curtis
##Prepare data and remove rows with zeros in all columns

dataA<-data[data$site=="A",] ###filtrar para Anrumai

dataA1<-dataA[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

dataA2<-dataA[,1:3] ###store aggregating variables to use latter
head(dataA2)


##############NMDS and bray curtis

table(is.na(dataA1)) ###check if there is any NAs

###subsetting for reduced sps numb

####### 75%
#nam<-m2[p75=="N",sp]
#str(nam)
dataAS<-dataA[,!..nam]
table(is.na(dataAS))


####### 95%
#nam1<-m2[p95=="N",sp]
#str(nam1)
dataAS1<-dataA[,!..nam1]
table(is.na(dataAS1))

####### sp totalling at least 1% of total Anrumaindance
#nam2<-m2[p1p=="N",sp]
#str(nam2)
dataAS2<-dataA[,!..nam2]
table(is.na(dataAS2))



##Remove cores with zero species after selection of p75, P95, P1P
dataASS<-dataAS[apply(dataAS[,!1:3],1,sum)!=0]
dataASS1<-dataAS1[apply(dataAS1[,!1:3],1,sum)!=0]
dataASS2<-dataAS2[apply(dataAS2[,!1:3],1,sum)!=0]

dataA1$dummy<-88.49558
dataASS$dummy<-88.49558
dataASS1$dummy<-88.49558
dataASS2$dummy<-88.49558

dataASS_1<-dataASS[,!1:3]
dataASS_2<-dataASS[,1:3]

dataASS1_1<-dataASS1[,!1:3]
dataASS1_2<-dataASS1[,1:3]

dataASS2_1<-dataASS2[,!1:3]
dataASS2_2<-dataASS2[,1:3]

####with dummy without log

set.seed(0)
NMDSA<-metaMDS(log(dataA1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSA)
stressplot(NMDSA)

#extract NMDS scores (x and y coordinates)
data.scoresA = as.data.table(scores(NMDSA))

#add columns to data frame 
data.scoresA$coreID = dataA2$coreID
data.scoresA$site = dataA2$site
data.scoresA$month = factor(dataA2$month)

head(data.scoresA)

#p<-c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","black")

PP_A<-ggplot(data.scoresA, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Anrumai, with log, with dummy")+
  theme_bw()

###fit sps data in Anrumai
fitA<-envfit(NMDSA,log(dataA1+1),permutations=999)
arrowA<-data.frame(fitA$vectors$arrows,R = fitA$vectors$r, P = fitA$vectors$pvals)
arrowA$FG <- rownames(arrowA)
arrowA.p<-arrowA[arrowA$P<=0.05,]

PP_A+
  geom_segment(data=arrowA.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowA.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


set.seed(10)
NMDSA_d75<-metaMDS(log(dataASS_1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSA_d75)
stressplot(NMDSA_d75)

#extract NMDS scores (x and y coordinates)
data.scoresA75 = as.data.table(scores(NMDSA_d75))

#add columns to data frame 
data.scoresA75$coreID = dataASS_2$coreID
data.scoresA75$site = dataASS_2$site
data.scoresA75$month = factor(dataASS_2$month)

head(data.scoresA75)

PP_A75<-ggplot(data.scoresA75, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p75 in Anrumai, with log, with dummy")+
  theme_bw()

###fit sps data in Anrumai
fitA75<-envfit(NMDSA_d75,log(dataASS_1+1),permutations=999)
arrowA75<-data.frame(fitA75$vectors$arrows,R = fitA75$vectors$r, P = fitA75$vectors$pvals)
arrowA75$FG <- rownames(arrowA75)
arrowA75.p<-arrowA75[arrowA75$P<=0.05,]

PP_A75+geom_segment(data=arrowA75.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowA75.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)




set.seed(11)
NMDSA_d95<-metaMDS(dataASS1_1,distance="bray",k=2,trymax=1000, autotransform = F)
beep()

plot(NMDSA_d95)
stressplot(NMDSA_d95)

#extract NMDS scores (x and y coordinates)
data.scoresA95 = as.data.table(scores(NMDSA_d95))

#add columns to data frame 
data.scoresA95$coreID = dataASS1_2$coreID
data.scoresA95$site = dataASS1_2$site
data.scoresA95$month = factor(dataASS1_2$month)

head(data.scoresA95)

PP_A95<-ggplot(data.scoresA95, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p95 in Anrumai, without log, with dummy")+
  theme_bw()

###fit sps data in Anrumai
fitA95<-envfit(NMDSA_d95,dataASS1_1,permutations=999)
arrowA95<-data.frame(fitA95$vectors$arrows,R = fitA95$vectors$r, P = fitA95$vectors$pvals)
arrowA95$FG <- rownames(arrowA95)
arrowA95.p<-arrowA95[arrowA95$P<=0.05,]

PP_A95+geom_segment(data=arrowA95.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowA95.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


set.seed(12)
NMDSA_dp1p<-metaMDS(dataASS2_1,distance="bray",k=2,trymax=1000,autotransform = F)
beep()


plot(NMDSA_dp1p)
stressplot(NMDSA_dp1p)

#extract NMDS scores (x and y coordinates)
data.scoresdp1p = as.data.table(scores(NMDSA_dp1p))

#add columns to data frame 
data.scoresdp1p$coreID = dataASS2_2$coreID
data.scoresdp1p$site = dataASS2_2$site
data.scoresdp1p$month = factor(dataASS2_2$month)

head(data.scoresdp1p)

PP_dp1p<-ggplot(data.scoresdp1p, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p1p in Anrumai, without log, with dummy")+
  theme_bw()

###fit sps data in Anrumai
fitdp1p<-envfit(NMDSA_dp1p,dataASS2_1,permutations=999)
arrowdp1p<-data.frame(fitdp1p$vectors$arrows,R = fitdp1p$vectors$r, P = fitdp1p$vectors$pvals)
arrowdp1p$FG <- rownames(arrowdp1p)
arrowdp1p.p<-arrowdp1p[arrowdp1p$P<=0.05,]

PP_dp1p+geom_segment(data=arrowdp1p.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowdp1p.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)

#####Abu

###NMDS with Bray curtis
##Prepare data and remove rows with zeros in all columns

dataAB<-data[data$site=="AB",] ###filtrar para Abu

dataAB1<-dataAB[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

dataAB2<-dataAB[,1:3] ###store aggregating variables to use latter
head(dataAB2)


##############NMDS and bray curtis

table(is.na(dataAB1)) ###check if there is any NABs

###subsetting for reduced sps numb

####### 75%
#nam<-m2[p75=="N",sp]
#str(nam)
dataABS<-dataAB[,!..nam]
table(is.na(dataABS))


####### 95%
#nam1<-m2[p95=="N",sp]
#str(nam1)
dataABS1<-dataAB[,!..nam1]
table(is.na(dataABS1))

####### sp totalling at least 1% of total Abundance
#nam2<-m2[p1p=="N",sp]
#str(nam2)
dataABS2<-dataAB[,!..nam2]
table(is.na(dataABS2))



##Remove cores with zero species after selection of p75, P95, P1P
dataABSS<-dataABS[apply(dataABS[,!1:3],1,sum)!=0]
dataABSS1<-dataABS1[apply(dataABS1[,!1:3],1,sum)!=0]
dataABSS2<-dataABS2[apply(dataABS2[,!1:3],1,sum)!=0]

dataAB1$dummy<-88.49558
dataABSS$dummy<-88.49558
dataABSS1$dummy<-88.49558
dataABSS2$dummy<-88.49558

dataABSS_1<-dataABSS[,!1:3]
dataABSS_2<-dataABSS[,1:3]

dataABSS1_1<-dataABSS1[,!1:3]
dataABSS1_2<-dataABSS1[,1:3]

dataABSS2_1<-dataABSS2[,!1:3]
dataABSS2_2<-dataABSS2[,1:3]

####with dummy without log

set.seed(0)
NMDSAB<-metaMDS(log(dataAB1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSAB)
stressplot(NMDSAB)

#extract NMDS scores (x and y coordinates)
data.scoresAB = as.data.table(scores(NMDSAB))

#add columns to data frame 
data.scoresAB$coreID = dataAB2$coreID
data.scoresAB$site = dataAB2$site
data.scoresAB$month = factor(dataAB2$month)

head(data.scoresAB)

#p<-c("#E41AB1C","#377EB8","#4DABF4AB","#984EAB3","#FF7F00","black")

PP_AB<-ggplot(data.scoresAB, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Abu, with log, with dummy")+
  theme_bw()

###fit sps data in Abu
fitAB<-envfit(NMDSAB,log(dataAB1+1),permutations=999)
arrowAB<-data.frame(fitAB$vectors$arrows,R = fitAB$vectors$r, P = fitAB$vectors$pvals)
arrowAB$FG <- rownames(arrowAB)
arrowAB.p<-arrowAB[arrowAB$P<=0.05,]

PP_AB+
  geom_segment(data=arrowAB.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowAB.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


set.seed(10)
NMDSAB_d75<-metaMDS(log(dataABSS_1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSAB_d75)
stressplot(NMDSAB_d75)

#extract NMDS scores (x and y coordinates)
data.scoresAB75 = as.data.table(scores(NMDSAB_d75))

#add columns to data frame 
data.scoresAB75$coreID = dataABSS_2$coreID
data.scoresAB75$site = dataABSS_2$site
data.scoresAB75$month = factor(dataABSS_2$month)

head(data.scoresAB75)

PP_AB75<-ggplot(data.scoresAB75, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p75 in Abu, with log, with dummy")+
  theme_bw()

###fit sps data in Abu
fitAB75<-envfit(NMDSAB_d75,log(dataABSS_1+1),permutations=999)
arrowAB75<-data.frame(fitAB75$vectors$arrows,R = fitAB75$vectors$r, P = fitAB75$vectors$pvals)
arrowAB75$FG <- rownames(arrowAB75)
arrowAB75.p<-arrowAB75[arrowAB75$P<=0.05,]

PP_AB75+geom_segment(data=arrowAB75.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowAB75.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)




set.seed(11)
NMDSAB_d95<-metaMDS(dataABSS1_1,distance="bray",k=2,trymax=1000, autotransform = F)
beep()

plot(NMDSAB_d95)
stressplot(NMDSAB_d95)

#extract NMDS scores (x and y coordinates)
data.scoresAB95 = as.data.table(scores(NMDSAB_d95))

#add columns to data frame 
data.scoresAB95$coreID = dataABSS1_2$coreID
data.scoresAB95$site = dataABSS1_2$site
data.scoresAB95$month = factor(dataABSS1_2$month)

head(data.scoresAB95)

PP_AB95<-ggplot(data.scoresAB95, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p95 in Abu, without log, with dummy")+
  theme_bw()

###fit sps data in Abu
fitAB95<-envfit(NMDSAB_d95,dataABSS1_1,permutations=999)
arrowAB95<-data.frame(fitAB95$vectors$arrows,R = fitAB95$vectors$r, P = fitAB95$vectors$pvals)
arrowAB95$FG <- rownames(arrowAB95)
arrowAB95.p<-arrowAB95[arrowAB95$P<=0.05,]

PP_AB95+geom_segment(data=arrowAB95.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowAB95.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


set.seed(12)
NMDSAB_dp1p<-metaMDS(dataABSS2_1,distance="bray",k=2,trymax=1000,autotransform = F)
beep()


plot(NMDSAB_dp1p)
stressplot(NMDSAB_dp1p)

#extract NMDS scores (x and y coordinates)
data.scoresdp1p = as.data.table(scores(NMDSAB_dp1p))

#add columns to data frame 
data.scoresdp1p$coreID = dataABSS2_2$coreID
data.scoresdp1p$site = dataABSS2_2$site
data.scoresdp1p$month = factor(dataABSS2_2$month)

head(data.scoresdp1p)

PP_dp1p<-ggplot(data.scoresdp1p, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p1p in Abu, without log, with dummy")+
  theme_bw()

###fit sps data in Abu
fitdp1p<-envfit(NMDSAB_dp1p,dataABSS2_1,permutations=999)
arrowdp1p<-data.frame(fitdp1p$vectors$arrows,R = fitdp1p$vectors$r, P = fitdp1p$vectors$pvals)
arrowdp1p$FG <- rownames(arrowdp1p)
arrowdp1p.p<-arrowdp1p[arrowdp1p$P<=0.05,]

PP_dp1p+geom_segment(data=arrowdp1p.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowdp1p.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


#####Escadinhas

###NMDS with Bray curtis
##Prepare data and remove rows with zeros in all columns

dataE<-data[data$site=="E",] ###filtrar para Escadinhas

dataE1<-dataE[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

dataE2<-dataE[,1:3] ###store aggregating variables to use latter
head(dataE2)


##############NMDS and bray curtis

table(is.na(dataE1)) ###check if there is any NEs

###subsetting for reduced sps numb

####### 75%
#nam<-m2[p75=="N",sp]
#str(nam)
dataES<-dataE[,!..nam]
table(is.na(dataES))


####### 95%
#nam1<-m2[p95=="N",sp]
#str(nam1)
dataES1<-dataE[,!..nam1]
table(is.na(dataES1))

####### sp totalling at least 1% of total Escadinhasndance
#nam2<-m2[p1p=="N",sp]
#str(nam2)
dataES2<-dataE[,!..nam2]
table(is.na(dataES2))



##Remove cores with zero species after selection of p75, P95, P1P
dataESS<-dataES[apply(dataES[,!1:3],1,sum)!=0]
dataESS1<-dataES1[apply(dataES1[,!1:3],1,sum)!=0]
dataESS2<-dataES2[apply(dataES2[,!1:3],1,sum)!=0]

dataE1$dummy<-88.49558
dataESS$dummy<-88.49558
dataESS1$dummy<-88.49558
dataESS2$dummy<-88.49558

dataESS_1<-dataESS[,!1:3]
dataESS_2<-dataESS[,1:3]

dataESS1_1<-dataESS1[,!1:3]
dataESS1_2<-dataESS1[,1:3]

dataESS2_1<-dataESS2[,!1:3]
dataESS2_2<-dataESS2[,1:3]

####with dummy without log

set.seed(0)
NMDSE<-metaMDS(log(dataE1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSE)
stressplot(NMDSE)

#extract NMDS scores (x and y coordinates)
data.scoresE = as.data.table(scores(NMDSE))

#add columns to data frame 
data.scoresE$coreID = dataE2$coreID
data.scoresE$site = dataE2$site
data.scoresE$month = factor(dataE2$month)

head(data.scoresE)

#p<-c("#E41E1C","#377EB8","#4DEF4E","#984EE3","#FF7F00","black")

PP_E<-ggplot(data.scoresE, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Escadinhas, with log, with dummy")+
  theme_bw()

###fit sps data in Escadinhas
fitE<-envfit(NMDSE,log(dataE1+1),permutations=999)
arrowE<-data.frame(fitE$vectors$arrows,R = fitE$vectors$r, P = fitE$vectors$pvals)
arrowE$FG <- rownames(arrowE)
arrowE.p<-arrowE[arrowE$P<=0.05,]

PP_E+
  geom_segment(data=arrowE.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowE.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


set.seed(10)
NMDSE_d75<-metaMDS(log(dataESS_1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSE_d75)
stressplot(NMDSE_d75)

#extract NMDS scores (x and y coordinates)
data.scoresE75 = as.data.table(scores(NMDSE_d75))

#add columns to data frame 
data.scoresE75$coreID = dataESS_2$coreID
data.scoresE75$site = dataESS_2$site
data.scoresE75$month = factor(dataESS_2$month)

head(data.scoresE75)

PP_E75<-ggplot(data.scoresE75, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p75 in Escadinhas, with log, with dummy")+
  theme_bw()

###fit sps data in Escadinhas
fitE75<-envfit(NMDSE_d75,log(dataESS_1+1),permutations=999)
arrowE75<-data.frame(fitE75$vectors$arrows,R = fitE75$vectors$r, P = fitE75$vectors$pvals)
arrowE75$FG <- rownames(arrowE75)
arrowE75.p<-arrowE75[arrowE75$P<=0.05,]

PP_E75+geom_segment(data=arrowE75.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowE75.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)



#####this one does not converge
set.seed(11)
NMDSE_d95<-metaMDS(dataESS1_1,distance="bray",k=2,trymax=1000, autotransform = F)
beep()

plot(NMDSE_d95)
stressplot(NMDSE_d95)

#extract NMDS scores (x and y coordinates)
data.scoresE95 = as.data.table(scores(NMDSE_d95))

#add columns to data frame 
data.scoresE95$coreID = dataESS1_2$coreID
data.scoresE95$site = dataESS1_2$site
data.scoresE95$month = factor(dataESS1_2$month)

head(data.scoresE95)

PP_E95<-ggplot(data.scoresE95, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p95 in Escadinhas, without log, with dummy")+
  theme_bw()

###fit sps data in Escadinhas
fitE95<-envfit(NMDSE_d95,dataESS1_1,permutations=999)
arrowE95<-data.frame(fitE95$vectors$arrows,R = fitE95$vectors$r, P = fitE95$vectors$pvals)
arrowE95$FG <- rownames(arrowE95)
arrowE95.p<-arrowE95[arrowE95$P<=0.05,]

PP_E95+geom_segment(data=arrowE95.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowE95.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


set.seed(12)
NMDSE_dp1p<-metaMDS(dataESS2_1,distance="bray",k=2,trymax=1000,autotransform = F)
beep()


plot(NMDSE_dp1p)
stressplot(NMDSE_dp1p)

#extract NMDS scores (x and y coordinates)
data.scoresdp1p = as.data.table(scores(NMDSE_dp1p))

#add columns to data frame 
data.scoresdp1p$coreID = dataESS2_2$coreID
data.scoresdp1p$site = dataESS2_2$site
data.scoresdp1p$month = factor(dataESS2_2$month)

head(data.scoresdp1p)

PP_dp1p<-ggplot(data.scoresdp1p, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p1p in Escadinhas, without log, with dummy")+
  theme_bw()

###fit sps data in Escadinhas
fitdp1p<-envfit(NMDSE_dp1p,dataESS2_1,permutations=999)
arrowdp1p<-data.frame(fitdp1p$vectors$arrows,R = fitdp1p$vectors$r, P = fitdp1p$vectors$pvals)
arrowdp1p$FG <- rownames(arrowdp1p)
arrowdp1p.p<-arrowdp1p[arrowdp1p$P<=0.05,]

PP_dp1p+geom_segment(data=arrowdp1p.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowdp1p.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)



#####Bijante

###NMDS with Bray curtis
##Prepare data and remove rows with zeros in all columns

dataBI<-data[data$site=="BI",] ###filtrar para Bijante

dataBI1<-dataBI[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

dataBI2<-dataBI[,1:3] ###store aggregating variables to use latter
head(dataBI2)


##############NMDS and bray curtis

table(is.na(dataBI1)) ###check if there is any NBIs

###subsetting for reduced sps numb

####### 75%
#nam<-m2[p75=="N",sp]
#str(nam)
dataBIS<-dataBI[,!..nam]
table(is.na(dataBIS))


####### 95%
#nam1<-m2[p95=="N",sp]
#str(nam1)
dataBIS1<-dataBI[,!..nam1]
table(is.na(dataBIS1))

####### sp totalling at least 1% of total Bijantendance
#nam2<-m2[p1p=="N",sp]
#str(nam2)
dataBIS2<-dataBI[,!..nam2]
table(is.na(dataBIS2))



##Remove cores with zero species after selection of p75, P95, P1P
dataBISS<-dataBIS[apply(dataBIS[,!1:3],1,sum)!=0]
dataBISS1<-dataBIS1[apply(dataBIS1[,!1:3],1,sum)!=0]
dataBISS2<-dataBIS2[apply(dataBIS2[,!1:3],1,sum)!=0]

dataBI1$dummy<-88.49558
dataBISS$dummy<-88.49558
dataBISS1$dummy<-88.49558
dataBISS2$dummy<-88.49558

dataBISS_1<-dataBISS[,!1:3]
dataBISS_2<-dataBISS[,1:3]

dataBISS1_1<-dataBISS1[,!1:3]
dataBISS1_2<-dataBISS1[,1:3]

dataBISS2_1<-dataBISS2[,!1:3]
dataBISS2_2<-dataBISS2[,1:3]

####with dummy without log

set.seed(0)
NMDSBI<-metaMDS(log(dataBI1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSBI)
stressplot(NMDSBI)

#extract NMDS scores (x and y coordinates)
data.scoresBI = as.data.table(scores(NMDSBI))

#add columns to data frame 
data.scoresBI$coreID = dataBI2$coreID
data.scoresBI$site = dataBI2$site
data.scoresBI$month = factor(dataBI2$month)

head(data.scoresBI)

#p<-c("#E41BI1C","#377EB8","#4DBIF4BI","#984EBI3","#FF7F00","black")

PP_BI<-ggplot(data.scoresBI, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Bijante, with log, with dummy")+
  theme_bw()

###fit sps data in Bijante
fitBI<-envfit(NMDSBI,log(dataBI1+1),permutations=999)
arrowBI<-data.frame(fitBI$vectors$arrows,R = fitBI$vectors$r, P = fitBI$vectors$pvals)
arrowBI$FG <- rownames(arrowBI)
arrowBI.p<-arrowBI[arrowBI$P<=0.05,]

PP_BI+
  geom_segment(data=arrowBI.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowBI.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


set.seed(10)
NMDSBI_d75<-metaMDS(log(dataBISS_1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSBI_d75)
stressplot(NMDSBI_d75)

#extract NMDS scores (x and y coordinates)
data.scoresBI75 = as.data.table(scores(NMDSBI_d75))

#add columns to data frame 
data.scoresBI75$coreID = dataBISS_2$coreID
data.scoresBI75$site = dataBISS_2$site
data.scoresBI75$month = factor(dataBISS_2$month)

head(data.scoresBI75)

PP_BI75<-ggplot(data.scoresBI75, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p75 in Bijante, with log, with dummy")+
  theme_bw()

###fit sps data in Bijante
fitBI75<-envfit(NMDSBI_d75,log(dataBISS_1+1),permutations=999)
arrowBI75<-data.frame(fitBI75$vectors$arrows,R = fitBI75$vectors$r, P = fitBI75$vectors$pvals)
arrowBI75$FG <- rownames(arrowBI75)
arrowBI75.p<-arrowBI75[arrowBI75$P<=0.05,]

PP_BI75+geom_segment(data=arrowBI75.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowBI75.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)




set.seed(11)
NMDSBI_d95<-metaMDS(dataBISS1_1,distance="bray",k=2,trymax=1000, autotransform = F)
beep()

plot(NMDSBI_d95)
stressplot(NMDSBI_d95)

#extract NMDS scores (x and y coordinates)
data.scoresBI95 = as.data.table(scores(NMDSBI_d95))

#add columns to data frame 
data.scoresBI95$coreID = dataBISS1_2$coreID
data.scoresBI95$site = dataBISS1_2$site
data.scoresBI95$month = factor(dataBISS1_2$month)

head(data.scoresBI95)

PP_BI95<-ggplot(data.scoresBI95, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p95 in Bijante, without log, with dummy")+
  theme_bw()

###fit sps data in Bijante
fitBI95<-envfit(NMDSBI_d95,dataBISS1_1,permutations=999)
arrowBI95<-data.frame(fitBI95$vectors$arrows,R = fitBI95$vectors$r, P = fitBI95$vectors$pvals)
arrowBI95$FG <- rownames(arrowBI95)
arrowBI95.p<-arrowBI95[arrowBI95$P<=0.05,]

PP_BI95+geom_segment(data=arrowBI95.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowBI95.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


set.seed(12)
NMDSBI_dp1p<-metaMDS(dataBISS2_1,distance="bray",k=2,trymax=1000,autotransform = F)
beep()


plot(NMDSBI_dp1p)
stressplot(NMDSBI_dp1p)

#extract NMDS scores (x and y coordinates)
data.scoresdp1p = as.data.table(scores(NMDSBI_dp1p))

#add columns to data frame 
data.scoresdp1p$coreID = dataBISS2_2$coreID
data.scoresdp1p$site = dataBISS2_2$site
data.scoresdp1p$month = factor(dataBISS2_2$month)

head(data.scoresdp1p)

PP_dp1p<-ggplot(data.scoresdp1p, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p1p in Bijante, without log, with dummy")+
  theme_bw()

###fit sps data in Bijante
fitdp1p<-envfit(NMDSBI_dp1p,dataBISS2_1,permutations=999)
arrowdp1p<-data.frame(fitdp1p$vectors$arrows,R = fitdp1p$vectors$r, P = fitdp1p$vectors$pvals)
arrowdp1p$FG <- rownames(arrowdp1p)
arrowdp1p.p<-arrowdp1p[arrowdp1p$P<=0.05,]

PP_dp1p+geom_segment(data=arrowdp1p.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowdp1p.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)



#####Bruce

###NMDS with Bray curtis
##Prepare data and remove rows with zeros in all columns

dataBR<-data[data$site=="BR",] ###filtrar para Bruce

dataBR1<-dataBR[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

dataBR2<-dataBR[,1:3] ###store aggregating variables to use latter
head(dataBR2)


##############NMDS and bray curtis

table(is.na(dataBR1)) ###check if there is any NBRs

###subsetting for reduced sps numb

####### 75%
#nam<-m2[p75=="N",sp]
#str(nam)
dataBRS<-dataBR[,!..nam]
table(is.na(dataBRS))


####### 95%
#nam1<-m2[p95=="N",sp]
#str(nam1)
dataBRS1<-dataBR[,!..nam1]
table(is.na(dataBRS1))

####### sp totalling at least 1% of total Brucendance
#nam2<-m2[p1p=="N",sp]
#str(nam2)
dataBRS2<-dataBR[,!..nam2]
table(is.na(dataBRS2))



##Remove cores with zero species after selection of p75, P95, P1P
dataBRSS<-dataBRS[apply(dataBRS[,!1:3],1,sum)!=0]
dataBRSS1<-dataBRS1[apply(dataBRS1[,!1:3],1,sum)!=0]
dataBRSS2<-dataBRS2[apply(dataBRS2[,!1:3],1,sum)!=0]

dataBR1$dummy<-88.49558
dataBRSS$dummy<-88.49558
dataBRSS1$dummy<-88.49558
dataBRSS2$dummy<-88.49558

dataBRSS_1<-dataBRSS[,!1:3]
dataBRSS_2<-dataBRSS[,1:3]

dataBRSS1_1<-dataBRSS1[,!1:3]
dataBRSS1_2<-dataBRSS1[,1:3]

dataBRSS2_1<-dataBRSS2[,!1:3]
dataBRSS2_2<-dataBRSS2[,1:3]

####with dummy without log

set.seed(0)
NMDSBR<-metaMDS(log(dataBR1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSBR)
stressplot(NMDSBR)

#extract NMDS scores (x and y coordinates)
data.scoresBR = as.data.table(scores(NMDSBR))

#add columns to data frame 
data.scoresBR$coreID = dataBR2$coreID
data.scoresBR$site = dataBR2$site
data.scoresBR$month = factor(dataBR2$month)

head(data.scoresBR)

#p<-c("#E41BR1C","#377EB8","#4DBRF4BR","#984EBR3","#FF7F00","black")

PP_BR<-ggplot(data.scoresBR, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Bruce, with log, with dummy")+
  theme_bw()

###fit sps data in Bruce
fitBR<-envfit(NMDSBR,log(dataBR1+1),permutations=999)
arrowBR<-data.frame(fitBR$vectors$arrows,R = fitBR$vectors$r, P = fitBR$vectors$pvals)
arrowBR$FG <- rownames(arrowBR)
arrowBR.p<-arrowBR[arrowBR$P<=0.05,]

PP_BR+
  geom_segment(data=arrowBR.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowBR.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


set.seed(10)
NMDSBR_d75<-metaMDS(log(dataBRSS_1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSBR_d75)
stressplot(NMDSBR_d75)

#extract NMDS scores (x and y coordinates)
data.scoresBR75 = as.data.table(scores(NMDSBR_d75))

#add columns to data frame 
data.scoresBR75$coreID = dataBRSS_2$coreID
data.scoresBR75$site = dataBRSS_2$site
data.scoresBR75$month = factor(dataBRSS_2$month)

head(data.scoresBR75)

PP_BR75<-ggplot(data.scoresBR75, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p75 in Bruce, with log, with dummy")+
  theme_bw()

###fit sps data in Bruce
fitBR75<-envfit(NMDSBR_d75,log(dataBRSS_1+1),permutations=999)
arrowBR75<-data.frame(fitBR75$vectors$arrows,R = fitBR75$vectors$r, P = fitBR75$vectors$pvals)
arrowBR75$FG <- rownames(arrowBR75)
arrowBR75.p<-arrowBR75[arrowBR75$P<=0.05,]

PP_BR75+geom_segment(data=arrowBR75.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowBR75.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)




set.seed(11)
NMDSBR_d95<-metaMDS(dataBRSS1_1,distance="bray",k=2,trymax=1000, autotransform = F)
beep()

plot(NMDSBR_d95)
stressplot(NMDSBR_d95)

#extract NMDS scores (x and y coordinates)
data.scoresBR95 = as.data.table(scores(NMDSBR_d95))

#add columns to data frame 
data.scoresBR95$coreID = dataBRSS1_2$coreID
data.scoresBR95$site = dataBRSS1_2$site
data.scoresBR95$month = factor(dataBRSS1_2$month)

head(data.scoresBR95)

PP_BR95<-ggplot(data.scoresBR95, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p95 in Bruce, without log, with dummy")+
  theme_bw()

###fit sps data in Bruce
fitBR95<-envfit(NMDSBR_d95,dataBRSS1_1,permutations=999)
arrowBR95<-data.frame(fitBR95$vectors$arrows,R = fitBR95$vectors$r, P = fitBR95$vectors$pvals)
arrowBR95$FG <- rownames(arrowBR95)
arrowBR95.p<-arrowBR95[arrowBR95$P<=0.05,]

PP_BR95+geom_segment(data=arrowBR95.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowBR95.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


set.seed(12)
NMDSBR_dp1p<-metaMDS(dataBRSS2_1,distance="bray",k=2,trymax=1000,autotransform = F)
beep()


plot(NMDSBR_dp1p)
stressplot(NMDSBR_dp1p)

#extract NMDS scores (x and y coordinates)
data.scoresdp1p = as.data.table(scores(NMDSBR_dp1p))

#add columns to data frame 
data.scoresdp1p$coreID = dataBRSS2_2$coreID
data.scoresdp1p$site = dataBRSS2_2$site
data.scoresdp1p$month = factor(dataBRSS2_2$month)

head(data.scoresdp1p)

PP_dp1p<-ggplot(data.scoresdp1p, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=month))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=month, colour=month))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS p1p in Bruce, without log, with dummy")+
  theme_bw()

###fit sps data in Bruce
fitdp1p<-envfit(NMDSBR_dp1p,dataBRSS2_1,permutations=999)
arrowdp1p<-data.frame(fitdp1p$vectors$arrows,R = fitdp1p$vectors$r, P = fitdp1p$vectors$pvals)
arrowdp1p$FG <- rownames(arrowdp1p)
arrowdp1p.p<-arrowdp1p[arrowdp1p$P<=0.05,]

PP_dp1p+geom_segment(data=arrowdp1p.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowdp1p.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


################### NMDS per month ################################

###########################Oct##############################################################

dataOct<-data[data$month==10] ###filtrar para Outubro

dataOct1<-dataOct[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

dataOct2<-dataOct[,1:3] ###store aggregating variables to use latter

##############NMDS and bray curtis

table(is.na(dataBR1)) ###check if there is any NBRs

###subsetting for reduced sps numb

####### 75%
#nam<-m2[p75=="N",sp]
#str(nam)
dataOctS<-dataOct[,!..nam]
table(is.na(dataOctS))


##Remove cores with zero species after selection of p75, P95, P1P
dataOctSS<-dataOctS[apply(dataOctS[,!1:3],1,sum)!=0]

dataOct1$dummy<-88.49558
dataOctSS$dummy<-88.49558

dataOctSS_1<-dataOctSS[,!1:3]
dataOctSS_2<-dataOctSS[,1:3]

####with dummy with log

########All sp
set.seed(0)
NMDSOct<-metaMDS(log(dataOct1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSOct)
stressplot(NMDSOct)

#extract NMDS scores (x and y coordinates)
data.scoresOct = as.data.table(scores(NMDSOct))

#add columns to data frame 
data.scoresOct$coreID = dataOct2$coreID
data.scoresOct$site = dataOct2$site
data.scoresOct$month = factor(dataOct2$month)

head(data.scoresOct)

#p<-c("#E41Oct1C","#377EB8","#4DOctF4Oct","#984EOct3","#FF7F00","black")

PP_Oct<-ggplot(data.scoresOct, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in October, with log, with dummy")+
  theme_bw()

###fit sps data in Bruce
fitOct<-envfit(NMDSOct,log(dataOct1+1),permutations=999)
arrowOct<-data.frame(fitOct$vectors$arrows,R = fitOct$vectors$r, P = fitOct$vectors$pvals)
arrowOct$FG <- rownames(arrowOct)
arrowOct.p<-arrowOct[arrowOct$P<=0.05,]

PP_Oct+
  geom_segment(data=arrowOct.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowOct.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)



########All sp
set.seed(1)
NMDSOct75<-metaMDS(log(dataOctSS_1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSOct75)
stressplot(NMDSOct75)

#extract NMDS scores (x and y coordinates)
data.scoresOct75 = as.data.table(scores(NMDSOct75))

#add columns to data frame 
data.scoresOct75$coreID = dataOctSS_2$coreID
data.scoresOct75$site = dataOctSS_2$site
data.scoresOct75$month = factor(dataOctSS_2$month)

head(data.scoresOct75)

#p<-c("#E41Oct1C","#377EB8","#4DOctF4Oct","#984EOct3","#FF7F00","black")

PP_Oct75<-ggplot(data.scoresOct75, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in October p75, with log, with dummy")+
  theme_bw()

###fit sps data in Bruce
fitOct75<-envfit(NMDSOct75,log(dataOctSS_1+1),permutations=999)
arrowOct75<-data.frame(fitOct75$vectors$arrows,R = fitOct75$vectors$r, P = fitOct75$vectors$pvals)
arrowOct75$FG <- rownames(arrowOct75)
arrowOct75.p<-arrowOct75[arrowOct75$P<=0.05,]

PP_Oct75+
  geom_segment(data=arrowOct75.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowOct75.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)
















###########################Nov##############################################################

dataNov<-data[data$month==12] ###filtrar para Outubro

dataNov1<-dataNov[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

dataNov2<-dataNov[,1:3] ###store aggregating variables to use latter

##############NMDS and bray curtis

table(is.na(dataBR1)) ###check if there is any NBRs

###subsetting for reduced sps numb

####### 75%
#nam<-m2[p75=="N",sp]
#str(nam)
dataNovS<-dataNov[,!..nam]
table(is.na(dataNovS))


##Remove cores with zero species after selection of p75, P95, P1P
dataNovSS<-dataNovS[apply(dataNovS[,!1:3],1,sum)!=0]

dataNov1$dummy<-88.49558
dataNovSS$dummy<-88.49558

dataNovSS_1<-dataNovSS[,!1:3]
dataNovSS_2<-dataNovSS[,1:3]

####with dummy with log

########All sp
set.seed(0)
NMDSNov<-metaMDS(log(dataNov1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSNov)
stressplot(NMDSNov)

#extract NMDS scores (x and y coordinates)
data.scoresNov = as.data.table(scores(NMDSNov))

#add columns to data frame 
data.scoresNov$coreID = dataNov2$coreID
data.scoresNov$site = dataNov2$site
data.scoresNov$month = factor(dataNov2$month)

head(data.scoresNov)

#p<-c("#E41Nov1C","#377EB8","#4DNovF4Nov","#984ENov3","#FF7F00","black")

PP_Nov<-ggplot(data.scoresNov, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Novober, with log, with dummy")+
  theme_bw()

###fit sps data in Bruce
fitNov<-envfit(NMDSNov,log(dataNov1+1),permutations=999)
arrowNov<-data.frame(fitNov$vectors$arrows,R = fitNov$vectors$r, P = fitNov$vectors$pvals)
arrowNov$FG <- rownames(arrowNov)
arrowNov.p<-arrowNov[arrowNov$P<=0.05,]

PP_Nov+
  geom_segment(data=arrowNov.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowNov.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)



########75 sp
set.seed(1)
NMDSNov75<-metaMDS(log(dataNovSS_1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSNov75)
stressplot(NMDSNov75)

#extract NMDS scores (x and y coordinates)
data.scoresNov75 = as.data.table(scores(NMDSNov75))

#add columns to data frame 
data.scoresNov75$coreID = dataNovSS_2$coreID
data.scoresNov75$site = dataNovSS_2$site
data.scoresNov75$month = factor(dataNovSS_2$month)

head(data.scoresNov75)

#p<-c("#E41Nov1C","#377EB8","#4DNovF4Nov","#984ENov3","#FF7F00","black")

PP_Nov75<-ggplot(data.scoresNov75, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Novober p75, with log, with dummy")+
  theme_bw()

###fit sps data in Bruce
fitNov75<-envfit(NMDSNov75,log(dataNovSS_1+1),permutations=999)
arrowNov75<-data.frame(fitNov75$vectors$arrows,R = fitNov75$vectors$r, P = fitNov75$vectors$pvals)
arrowNov75$FG <- rownames(arrowNov75)
arrowNov75.p<-arrowNov75[arrowNov75$P<=0.05,]

PP_Nov75+
  geom_segment(data=arrowNov75.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowNov75.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)











###########################Dez##############################################################

dataDez<-data[data$month==11] ###filtrar para Outubro

dataDez1<-dataDez[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

dataDez2<-dataDez[,1:3] ###store aggregating variables to use latter

##############NMDS and bray curtis

table(is.na(dataBR1)) ###check if there is any NBRs

###subsetting for reduced sps numb

####### 75%
#nam<-m2[p75=="N",sp]
#str(nam)
dataDezS<-dataDez[,!..nam]
table(is.na(dataDezS))


##Remove cores with zero species after selection of p75, P95, P1P
dataDezSS<-dataDezS[apply(dataDezS[,!1:3],1,sum)!=0]

dataDez1$dummy<-88.49558
dataDezSS$dummy<-88.49558

dataDezSS_1<-dataDezSS[,!1:3]
dataDezSS_2<-dataDezSS[,1:3]

####with dummy with log

########All sp
set.seed(0)
NMDSDez<-metaMDS(log(dataDez1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSDez)
stressplot(NMDSDez)

#extract NMDS scores (x and y coordinates)
data.scoresDez = as.data.table(scores(NMDSDez))

#add columns to data frame 
data.scoresDez$coreID = dataDez2$coreID
data.scoresDez$site = dataDez2$site
data.scoresDez$month = factor(dataDez2$month)

head(data.scoresDez)

#p<-c("#E41Dez1C","#377EB8","#4DDezF4Dez","#984EDez3","#FF7F00","black")

PP_Dez<-ggplot(data.scoresDez, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Dezober, with log, with dummy")+
  theme_bw()

###fit sps data in Bruce
fitDez<-envfit(NMDSDez,log(dataDez1+1),permutations=999)
arrowDez<-data.frame(fitDez$vectors$arrows,R = fitDez$vectors$r, P = fitDez$vectors$pvals)
arrowDez$FG <- rownames(arrowDez)
arrowDez.p<-arrowDez[arrowDez$P<=0.05,]

PP_Dez+
  geom_segment(data=arrowDez.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowDez.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)



########All sp
set.seed(1)
NMDSDez75<-metaMDS(log(dataDezSS_1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSDez75)
stressplot(NMDSDez75)

#extract NMDS scores (x and y coordinates)
data.scoresDez75 = as.data.table(scores(NMDSDez75))

#add columns to data frame 
data.scoresDez75$coreID = dataDezSS_2$coreID
data.scoresDez75$site = dataDezSS_2$site
data.scoresDez75$month = factor(dataDezSS_2$month)

head(data.scoresDez75)

#p<-c("#E41Dez1C","#377EB8","#4DDezF4Dez","#984EDez3","#FF7F00","black")

PP_Dez75<-ggplot(data.scoresDez75, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Dezober p75, with log, with dummy")+
  theme_bw()

###fit sps data in Bruce
fitDez75<-envfit(NMDSDez75,log(dataDezSS_1+1),permutations=999)
arrowDez75<-data.frame(fitDez75$vectors$arrows,R = fitDez75$vectors$r, P = fitDez75$vectors$pvals)
arrowDez75$FG <- rownames(arrowDez75)
arrowDez75.p<-arrowDez75[arrowDez75$P<=0.05,]

PP_Dez75+
  geom_segment(data=arrowDez75.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowDez75.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)





###########################Jan##############################################################

dataJan<-data[data$month==1] ###filtrar para Outubro

dataJan1<-dataJan[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

dataJan2<-dataJan[,1:3] ###store aggregating variables to use latter

##############NMDS and bray curtis

table(is.na(dataBR1)) ###check if there is any NBRs

###subsetting for reduced sps numb

####### 75%
#nam<-m2[p75=="N",sp]
#str(nam)
dataJanS<-dataJan[,!..nam]
table(is.na(dataJanS))


##Remove cores with zero species after selection of p75, P95, P1P
dataJanSS<-dataJanS[apply(dataJanS[,!1:3],1,sum)!=0]

dataJan1$dummy<-88.49558
dataJanSS$dummy<-88.49558

dataJanSS_1<-dataJanSS[,!1:3]
dataJanSS_2<-dataJanSS[,1:3]

####with dummy with log

########All sp
set.seed(0)
NMDSJan<-metaMDS(log(dataJan1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSJan)
stressplot(NMDSJan)

#extract NMDS scores (x and y coordinates)
data.scoresJan = as.data.table(scores(NMDSJan))

#add columns to data frame 
data.scoresJan$coreID = dataJan2$coreID
data.scoresJan$site = dataJan2$site
data.scoresJan$month = factor(dataJan2$month)

head(data.scoresJan)

#p<-c("#E41Jan1C","#377EB8","#4DJanF4Jan","#984EJan3","#FF7F00","black")

PP_Jan<-ggplot(data.scoresJan, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Jan, with log, with dummy")+
  theme_bw()

###fit sps data in Bruce
fitJan<-envfit(NMDSJan,log(dataJan1+1),permutations=999)
arrowJan<-data.frame(fitJan$vectors$arrows,R = fitJan$vectors$r, P = fitJan$vectors$pvals)
arrowJan$FG <- rownames(arrowJan)
arrowJan.p<-arrowJan[arrowJan$P<=0.05,]

PP_Jan+
  geom_segment(data=arrowJan.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowJan.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)



########75 sp
set.seed(1)
NMDSJan75<-metaMDS(log(dataJanSS_1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSJan75)
stressplot(NMDSJan75)

#extract NMDS scores (x and y coordinates)
data.scoresJan75 = as.data.table(scores(NMDSJan75))

#add columns to data frame 
data.scoresJan75$coreID = dataJanSS_2$coreID
data.scoresJan75$site = dataJanSS_2$site
data.scoresJan75$month = factor(dataJanSS_2$month)

head(data.scoresJan75)

#p<-c("#E41Jan1C","#377EB8","#4DJanF4Jan","#984EJan3","#FF7F00","black")

PP_Jan75<-ggplot(data.scoresJan75, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Janober p75, with log, with dummy")+
  theme_bw()

###fit sps data in Bruce
fitJan75<-envfit(NMDSJan75,log(dataJanSS_1+1),permutations=999)
arrowJan75<-data.frame(fitJan75$vectors$arrows,R = fitJan75$vectors$r, P = fitJan75$vectors$pvals)
arrowJan75$FG <- rownames(arrowJan75)
arrowJan75.p<-arrowJan75[arrowJan75$P<=0.05,]

PP_Jan75+
  geom_segment(data=arrowJan75.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowJan75.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)
















###########################Feb##############################################################

dataFeb<-data[data$month==2] ###filtrar para Outubro

dataFeb1<-dataFeb[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

dataFeb2<-dataFeb[,1:3] ###store aggregating variables to use latter

##############NMDS and bray curtis

table(is.na(dataBR1)) ###check if there is any NBRs

###subsetting for reduced sps numb

####### 75%
#nam<-m2[p75=="N",sp]
#str(nam)
dataFebS<-dataFeb[,!..nam]
table(is.na(dataFebS))


##Remove cores with zero species after selection of p75, P95, P1P
dataFebSS<-dataFebS[apply(dataFebS[,!1:3],1,sum)!=0]

dataFeb1$dummy<-88.49558
dataFebSS$dummy<-88.49558

dataFebSS_1<-dataFebSS[,!1:3]
dataFebSS_2<-dataFebSS[,1:3]

####with dummy with log

########All sp
set.seed(0)
NMDSFeb<-metaMDS(log(dataFeb1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSFeb)
stressplot(NMDSFeb)

#extract NMDS scores (x and y coordinates)
data.scoresFeb = as.data.table(scores(NMDSFeb))

#add columns to data frame 
data.scoresFeb$coreID = dataFeb2$coreID
data.scoresFeb$site = dataFeb2$site
data.scoresFeb$month = factor(dataFeb2$month)

head(data.scoresFeb)

#p<-c("#E41Feb1C","#377EB8","#4DFebF4Feb","#984EFeb3","#FF7F00","black")

PP_Feb<-ggplot(data.scoresFeb, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Feb, with log, with dummy")+
  theme_bw()

###fit sps data in Bruce
fitFeb<-envfit(NMDSFeb,log(dataFeb1+1),permutations=999)
arrowFeb<-data.frame(fitFeb$vectors$arrows,R = fitFeb$vectors$r, P = fitFeb$vectors$pvals)
arrowFeb$FG <- rownames(arrowFeb)
arrowFeb.p<-arrowFeb[arrowFeb$P<=0.05,]

PP_Feb+
  geom_segment(data=arrowFeb.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowFeb.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)




###########################Mar##############################################################

dataMar<-data[data$month==3] ###filtrar para Outubro

dataMar1<-dataMar[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

dataMar2<-dataMar[,1:3] ###store aggregating variables to use latter

##############NMDS and bray curtis

table(is.na(dataBR1)) ###check if there is any NBRs

###subsetting for reduced sps numb

####### 75%
#nam<-m2[p75=="N",sp]
#str(nam)
dataMarS<-dataMar[,!..nam]
table(is.na(dataMarS))


##Remove cores with zero species after selection of p75, P95, P1P
dataMarSS<-dataMarS[apply(dataMarS[,!1:3],1,sum)!=0]

dataMar1$dummy<-88.49558
dataMarSS$dummy<-88.49558

dataMarSS_1<-dataMarSS[,!1:3]
dataMarSS_2<-dataMarSS[,1:3]

####with dummy with log

########All sp
set.seed(0)
NMDSMar<-metaMDS(log(dataMar1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSMar)
stressplot(NMDSMar)

#extract NMDS scores (x and y coordinates)
data.scoresMar = as.data.table(scores(NMDSMar))

#add columns to data frame 
data.scoresMar$coreID = dataMar2$coreID
data.scoresMar$site = dataMar2$site
data.scoresMar$month = factor(dataMar2$month)

head(data.scoresMar)

#p<-c("#E41Mar1C","#377EB8","#4DMarF4Mar","#984EMar3","#FF7F00","black")

PP_Mar<-ggplot(data.scoresMar, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Mar, with log, with dummy")+
  theme_bw()

###fit sps data in Bruce
fitMar<-envfit(NMDSMar,log(dataMar1+1),permutations=999)
arrowMar<-data.frame(fitMar$vectors$arrows,R = fitMar$vectors$r, P = fitMar$vectors$pvals)
arrowMar$FG <- rownames(arrowMar)
arrowMar.p<-arrowMar[arrowMar$P<=0.05,]

PP_Mar+
  geom_segment(data=arrowMar.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowMar.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)










###########################Apr##############################################################

dataApr<-data[data$month==4] ###filtrar para Outubro

dataApr1<-dataApr[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

dataApr2<-dataApr[,1:3] ###store aggregating variables to use latter

##############NMDS and bray curtis

table(is.na(dataBR1)) ###check if there is any NBRs

###subsetting for reduced sps numb

####### 75%
#nam<-m2[p75=="N",sp]
#str(nam)
dataAprS<-dataApr[,!..nam]
table(is.na(dataAprS))


##Remove cores with zero species after selection of p75, P95, P1P
dataAprSS<-dataAprS[apply(dataAprS[,!1:3],1,sum)!=0]

dataApr1$dummy<-88.49558
dataAprSS$dummy<-88.49558

dataAprSS_1<-dataAprSS[,!1:3]
dataAprSS_2<-dataAprSS[,1:3]

####with dummy with log

########All sp
set.seed(0)
NMDSApr<-metaMDS(log(dataApr1+1),distance="bray",k=3,trymax=1000,autotransform = F)
beep()

plot(NMDSApr)
stressplot(NMDSApr)

#extract NMDS scores (x and y coordinates)
data.scoresApr = as.data.table(scores(NMDSApr))

#add columns to data frame 
data.scoresApr$coreID = dataApr2$coreID
data.scoresApr$site = dataApr2$site
data.scoresApr$month = factor(dataApr2$month)

head(data.scoresApr)

#p<-c("#E41Apr1C","#377EB8","#4DAprF4Apr","#984EApr3","#FF7F00","black")

PP_Apr<-ggplot(data.scoresApr, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Month")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site))+
  #scale_colour_manual(values=p)+
  ggtitle("NMDS in Apr, with log, with dummy")+
  theme_bw()

###fit sps data in Bruce
fitApr<-envfit(NMDSApr,log(dataApr1+1),permutations=999)
arrowApr<-data.frame(fitApr$vectors$arrows,R = fitApr$vectors$r, P = fitApr$vectors$pvals)
arrowApr$FG <- rownames(arrowApr)
arrowApr.p<-arrowApr[arrowApr$P<=0.05,]

PP_Apr+
  geom_segment(data=arrowApr.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  ggrepel::geom_text_repel(data=arrowApr.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),cex=5,direction="both",segment.size=0.25)


