setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("vegan","ggplot2","viridis","RColorBrewer","psych","reshape2","beepr","data.table")
lapply(packs,require,character.only=T)

## Load database 
DB66<-fread("data_out/db/Final_DB_lowtaxa_density_polyexcl_20210121.csv") ### created in script called Database_cleanup_joining
str(DB66)


### names to change

# DONE Remove Corbula_sulcata

# DONE Verificar que nova base de dadso tem Skenidae e Rissoidae corrigidos dos typos do paulino

# DONE Phyllocida passa a chamar-se Phyllodocida1
# DONE Correct family name of Pseudopythina_nicklesi (it has a space) DONEEE


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
set.seed(1)
NMDS_d75<-metaMDS(dataSS_1,distance="bray",k=2,trymax=1000)
beep()

set.seed(2)
NMDS_d95<-metaMDS(dataSS1_1,distance="bray",k=2,trymax=1000)
beep()

set.seed(3)
NMDS_dp1p<-metaMDS(dataSS2_1,distance="bray",k=2,trymax=1000)
beep()

####with dummy with log
set.seed(4)
NMDS_d75log<-metaMDS(dataSS_1log,distance="bray",k=2,trymax=1000)
beep()

set.seed(5)
NMDS_d95log<-metaMDS(dataSS1_1log,distance="bray",k=2,trymax=1000)
beep()

set.seed(6)
NMDS_dp1plog<-metaMDS(dataSS2_1log,distance="bray",k=2,trymax=1000)
beep()


###ploting
plot(NMDS)
stressplot(NMDS)


#extract NMDS scores (x and y coordinates)
data.scores = as.data.table(scores(NMDS_d95))

#add columns to data frame 
data.scores$coreID = dataSS1_2$coreID
data.scores$site = dataSS1_2$site
data.scores$month = dataSS1_2$month

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
  ggtitle("NMDS d95, without log, with dummy")+
  theme_bw()


###fit sps data
fit<-envfit(NMDS_d95,dataSS1_1, permutations=999)
arrow<-data.frame(fit$vectors$arrows,R = fit$vectors$r, P = fit$vectors$pvals)
arrow$FG <- rownames(arrow)
arrow.p<-arrow[arrow$P<=0.05,]


ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Site")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site))+
  scale_colour_manual(values=p)+
  ggtitle("NMDS d95, without log, with dummy")+
  geom_segment(data=arrow.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5,label=FG),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  geom_text(data=arrow.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),size=5)+
  theme_bw()


################separating sites #########################################################################
#####Adonga

###NMDS with Bray curtis
##Prepare data and remove rows with zeros in all columns

dataAD<-data[data$site=="AD",] ###filtrar para adonga

data1<-data[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

data2<-data[,1:3] ###store aggregating variables to use latter
head(data2)


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
set.seed(1)
NMDS_d75<-metaMDS(dataSS_1,distance="bray",k=2,trymax=1000)
beep()

set.seed(2)
NMDS_d95<-metaMDS(dataSS1_1,distance="bray",k=2,trymax=1000)
beep()

set.seed(3)
NMDS_dp1p<-metaMDS(dataSS2_1,distance="bray",k=2,trymax=1000)
beep()

####with dummy with log
set.seed(4)
NMDS_d75log<-metaMDS(dataSS_1log,distance="bray",k=2,trymax=1000)
beep()

set.seed(5)
NMDS_d95log<-metaMDS(dataSS1_1log,distance="bray",k=2,trymax=1000)
beep()

set.seed(6)
NMDS_dp1plog<-metaMDS(dataSS2_1log,distance="bray",k=2,trymax=1000)
beep()


###ploting
plot(NMDS)
stressplot(NMDS)


#extract NMDS scores (x and y coordinates)
data.scores = as.data.table(scores(NMDS_d95))

#add columns to data frame 
data.scores$coreID = dataSS1_2$coreID
data.scores$site = dataSS1_2$site
data.scores$month = dataSS1_2$month

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
  ggtitle("NMDS d95, without log, with dummy")+
  theme_bw()


###fit sps data
fit<-envfit(NMDS_d95,dataSS1_1, permutations=999)
arrow<-data.frame(fit$vectors$arrows,R = fit$vectors$r, P = fit$vectors$pvals)
arrow$FG <- rownames(arrow)
arrow.p<-arrow[arrow$P<=0.05,]


ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4,aes(colour=site))+
  #geom_text(aes(label=month),size=7)+
  labs(x = "NMDS1", y = "NMDS2",colour="Site")+
  stat_ellipse(size=2, type="t",aes(group=site, colour=site))+
  scale_colour_manual(values=p)+
  ggtitle("NMDS d95, without log, with dummy")+
  geom_segment(data=arrow.p, aes(x=0,y=0,xend=NMDS1*R*4.5,yend=NMDS2*R*4.5,label=FG),arrow=arrow(length=unit(.2,"cm")),col="grey40",lwd=1)+
  geom_text(data=arrow.p,aes(x=NMDS1*R*4.5,y=NMDS2*R*4.5,label=FG),size=5)+
  theme_bw()


