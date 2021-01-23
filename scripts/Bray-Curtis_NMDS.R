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



###Limit database to the 75% most abundant species
perc<-function(x){x/sum(x)}
test<-DB1[,!c("coreID","month")][,lapply(.SD,mean),by=site]
test1<-sapply(test,perc)
#nam<-names(test)[2:ncol(test)]
test[,lapply(.SD,perc),by=site]
rows

###



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
####regular bc

table(is.na(data1)) ###check if there is any NAs

###subsetting for reduced sps numb
nam<-m2[p75=="N",sp]
str(nam)

dataS<-data[,!..nam]
table(is.na(dataS))


##Remove cores with zero species after selection of p75
dataS1<-dataS[apply(dataS[,!1:3],1,sum)!=0]



dataS1$dummy<-88.49558

dataS1_1<-dataS1[,!1:3]
dataS1_2<-dataS1[,1:3]

ff<-function(x){log(x+1)}
dataS1_1log<-dataS1_1[,lapply(.SD,ff)]


set.seed(2)
NMDS<-metaMDS(dataS1_1log,distance="bray",k=2,trymax=1000)
beep()

###ploting
plot(NMDS)
stressplot(NMDS)


#extract NMDS scores (x and y coordinates)
data.scores = as.data.table(scores(NMDS))

#add columns to data frame 
data.scores$coreID = dataS1_2$coreID
data.scores$site = dataS1_2$site
data.scores$month = dataS1_2$month

head(data.scores)


ggplot(data.scores, aes(x = NMDS1, y = NMDS2, colour=site,shape=factor(month),group=site)) + 
  geom_point(size = 4)+
  #stat_summary(geom="pointrange",size = 1, aes(colour=factor(month),shape=site))+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", shape = "Month", y = "NMDS2",colour="Site")+
  stat_ellipse(size=2)+
  scale_color_brewer(palette="Set1")+
  theme_bw()






###trial adjusted bray
data1$dummy<-127.3885

data1$dummy

table(is.na(data1)) ###check if there is any NAs

### log transform

ff<-function(x){log(x+1)}
data11<-data1[,lapply(.SD,ff)]



set.seed(100)
NMDS<-metaMDS(data1,distance="bray",k=2,trymax=1000)
beep()

plot(NMDS)
stressplot(NMDS)



#extract NMDS scores (x and y coordinates)
data.scores = as.data.table(scores(NMDS))

#add columns to data frame 
data.scores$coreID = data2$coreID
data.scores$site = data2$site
data.scores$month = data2$month

head(data.scores)


ggplot(data.scores, aes(x = NMDS1, y = NMDS2, colour=site,shape=factor(month),group=site)) + 
  geom_point(size = 4)+
  #stat_summary(geom="pointrange",size = 1, aes(colour=factor(month),shape=site))+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", shape = "Month", y = "NMDS2",colour="Site")+
  stat_ellipse(size=2)+
  scale_color_brewer(palette="Set1")+
  theme_bw()



ordiplot(NMDS,type="n")
orditorp(NMDS,display="species",col="red",air=0.01)
orditorp(NMDS,display="sites",cex=0.5,air=0.01)


ordiplot(NMDS,type="n")
ordihull(NMDS,groups=data2$site,draw="polygon",col="grey90",label=F)
orditorp(NMDS,display="species",col="red",air=0.01)
orditorp(NMDS,display="sites",col=c(rep("green",5),rep("blue",5)),
         air=0.01,cex=1)
