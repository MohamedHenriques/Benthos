setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("data.table","gridExtra","vegan","ggplot2","viridis","effects","RColorBrewer","xlsx","psych","reshape2","tidyr")
lapply(packs,require,character.only=T)


## MH database
db_mh<-fread("C:/Doutoramento1/Capitulos/Exclosure_experiments/Data/Cores_ID/DBMH/Cores_DB_all_20201214.csv")
db_mh0<-fread("C:/Doutoramento1/Capitulos/Exclosure_experiments/Data/Cores_ID/DBMH/Measures_Catarina/Cores_DB_all_20210131_CR.csv")
#db_mh<-read.table("D:/Work/FCUL/Doutoramento/Capitulos/Exclosure_experiments/Data/Cores_ID/DBMH/Cores_DB_all_20201126.csv",sep=";",header = T)
str(db_mh0)


unique(db_mh$low_taxa)

db_mh1<-db_mh0[treatment=="C",c(2:5,7:9,11:12,15,19:20,23,29)]
db_mh2<-db_mh1[done=="Y",]
str(db_mh2)

### Complete length column for rows with missing data
db_mh2[,unique(lenght)]
db_mh2[,as.numeric(unique(lenght))]
db_mh2[is.na(as.numeric(lenght)),.N]
db_mh2[,lenght:=as.numeric(lenght)]

db_mh2[,.(table(is.na(lenght)))]
db_mh2[,site:="AD"]
db_mh2[,coreID_ok:=paste(coreID,round,sep="")]
unique(db_mh2$coreID_ok)

db_mh2[is.na(lenght),.(table(low_taxa))]
db_mh2[,mean(lenght,na.rm=T),by=low_taxa]


X<-db_mh2[is.na(lenght),unique(low_taxa)]

#i=49
for(i in X){
  db_mh2[is.na(lenght)&low_taxa==i,lenght:=db_mh2[low_taxa==i,mean(lenght,na.rm=T)]]
}

db_mh2[is.na(lenght),.(table(low_taxa))]
db_mh2[is.na(lenght),unique(low_taxa)]



#### Correct names of taxa
unique(db_mh2$class)
db_mh2$class1<-as.character(db_mh2$class)
db_mh2$class1[which(db_mh2$subclass=="Sedentaria")]<-"Polychaeta_sedentaria"
db_mh2$class1[which(db_mh2$subclass=="Errantia")]<-"Polychaeta_errantia"
#db_mh2$class1[which(db_mh2$subclass=="NID")]<-"Polychaeta_NID" ## ja nao existem casos destes, mudei tudo paa sedentaria na base de dados original
db_mh2$class1[which(db_mh2$class1=="Enteropneusta")]<-"Other"
db_mh2$class1[which(db_mh2$class1=="NID")]<-"Other"
db_mh2$island<-as.character("Orango")

unique(db_mh2$class1)
db_mh2[db_mh2$class1=="Enteropneusta",]

#db_mh2[!(low_taxa== )]



###Ana's database
db_ac<-fread("C:/Doutoramento1/Capitulos/Spatial_temporal_variability_invertebrates/Data_Ana/Cores_V25_20201214.csv")
str(db_ac)
db_ac1<-db_ac[,c(2,4,5:7,9,11,15:16,27)]
db_ac1$numb<-ifelse(db_ac1$final_name=="empty",0,1)
db_ac1$coreID<-substr(db_ac1$point,1,nchar(as.character(db_ac1$point))-1)
db_ac1$class1<-ifelse(db_ac1$class=="BIV","Bivalvia",
                      ifelse(db_ac1$class=="SED","Polychaeta_sedentaria",
                             ifelse(db_ac1$class=="ERR","Polychaeta_errantia",
                                    ifelse(db_ac1$class=="GAS","Gastropoda",
                                           ifelse(db_ac1$class=="CRU","Malacostraca",
                                                  ifelse(db_ac1$class=="OTH","Other","empty"))))))
unique(db_ac1$site)
db_ac1$island<-ifelse(db_ac1$site=="AB"|db_ac1$site=="A","Formosa",ifelse(db_ac1$site=="E"|db_ac1$site=="BI"|db_ac1$site=="BR","Bubaque",NA))

length(unique(db_ac1$coreID))
unique(db_ac1$position)

### Check if there is any sps for which in MH database there are no records of size but are on the list for biomass estimates

##LIst of sps for biomass estimates
Z<-c("Amphipoda1","Cirratulidae","Diplodonta","Keletistes_sp","Maldanidae","Diopatra_sp","Orbiniidae","Paraonidae",
     "Moerella_distorta","Capitellidae","Goniadopsis_sp","Grandidierella_sp","Hyala_sp","Sigambra_sp","Solariella_sp",
     "Pelecyora_isocardia","Tagelus_adansonii","Austromacoma_nymphalis","Senilia_senilis","Nereididae","Glycera_sp",
     "Marphysa_sanguinea","Lumbrineridae","Afruca_tangeri","Balsscallichirus_balssi","Anthuridae")

match<-Z[Z%in%db_mh2[is.na(lenght),unique(low_taxa)]] ##only one match

str(db_ac1)
X1<-db_ac1[,mean(new.tl,na.rm=T),by=final_name]

##For each match case found, attribute the mean size found in Ana's samples to MH samples
for(i in match){
  db_mh2[is.na(lenght)&low_taxa==i,lenght:=X1[final_name==i,V1]]
}
db_mh2[low_taxa==match[1],lenght]
X1[final_name==match[1],V1]

db_mh2[low_taxa==match[2],lenght]
X1[final_name==match[2],V1]

Z[Z%in%db_mh2[is.na(lenght),unique(low_taxa)]] ##No more cases


####Paulino's database now

db_jp<-fread("C:/Doutoramento1/Capitulos/Spatial_temporal_variability_invertebrates/Data_Paulino/Invert_Adonga_2019_JPaulino_corrected_20201214.csv")
db_jp1<-db_jp[which(db_jp$Tipo =="Nuca")]

db_jp_a2<-db_jp1[which(db_jp1$Area ==2)]

db_jp2<-db_jp1[-which(db_jp1$Area ==2),]
db_jp3<-db_jp2[-which(db_jp2$Core>10),]
db_jp3$coreID<-paste("JP",db_jp3$Area,db_jp3$Core,sep="_")
db_jp3$site<-as.character("AD")
db_jp3$year<-2019
db_jp3$month<-2
db_jp3$island<-"Orango"
db_jp3$class1<-as.character(db_jp3$Class)
db_jp3$class1[which(db_jp3$class1=="Oligochaeta")]<-"Other"
db_jp3$class1[which(db_jp3$class1=="Polychaeta_incertae_sedis")]<-"Polychaeta_sedentaria"
db_jp3$depth<-ifelse(db_jp3$Profundidade=="Fundo","B","T")

unique(db_jp3$coreID)
unique(db_jp3$class1)
unique(db_jp3$Lowest_taxa)
unique(db_jp3$Family)


### Uniformizar
data.frame(names(db_mh2))
db_mh3<-db_mh2[,c(15,2:3,16,6,8,17,10:12,13)][,AFDW_mg:=NA]
colnames(db_mh3)[4]<-"coreID"

data.frame(names(db_ac1))
db_ac2<-db_ac1[,c(3,1,2,12,5,7,13,9,8,11,10)][,AFDW_mg:=NA]
names(db_ac2)<-names(db_mh3)

db_jp4_1<-db_jp3[,c("site","year","month","coreID","depth","Class","class1","Family","Lowest_taxa","n","AFDW_mg")]
db_jp4_1[,lenght:=NA]
db_jp4<-db_jp4_1[,c(c("site","year","month","coreID","depth","Class","class1","Family","Lowest_taxa","n","lenght","AFDW_mg"))]
names(db_jp4)<-names(db_mh3)

### aggregate to have one line per taxa per core per all other stuff
db_mh4<-as.data.table(aggregate(db_mh3$numb,by=list(site=db_mh3$site,year=db_mh3$year,month=db_mh3$month,coreID=db_mh3$coreID,depth=db_mh3$depth,
                               class1=db_mh3$class1,family=db_mh3$family,low_taxa=db_mh3$low_taxa,lenght=as.character(db_mh3$lenght)),FUN=sum))
str(db_mh4)
#db_mh4$coreID<-as.character(db_mh4$coreID)
# db_mh4$depth<-as.character(db_mh4$depth)
# db_mh4$family<-as.character(db_mh4$family)
# db_mh4$low_taxa<-as.character(db_mh4$low_taxa)
db_mh4$lenght<-as.numeric(db_mh4$lenght)
db_mh4$AFDW_mg<-NA

db_ac3<-as.data.table(aggregate(db_ac2$numb,by=list(site=db_ac2$site,year=db_ac2$year,month=db_ac2$month,coreID=db_ac2$coreID,depth=db_ac2$depth,
                                      class1=db_ac2$class1,family=db_ac2$family,low_taxa=db_ac2$low_taxa,lenght=as.character(db_ac2$lenght)),FUN=sum))

db_ac3$AFDW_mg<-NA

db_jp5<-as.data.table(aggregate(db_jp4$numb,by=list(site=db_jp4$site,year=db_jp4$year,month=db_jp4$month,coreID=db_jp4$coreID,depth=db_jp4$depth,
                                      class1=db_jp4$class1,family=db_jp4$family,low_taxa=db_jp4$low_taxa,AFDW_mg=as.character(db_jp4$AFDW_mg)),FUN=sum))

db_jp5[,AFDW_mg:=as.numeric(AFDW_mg)]
db_jp5[,lenght:=NA]


unique(db_ac3$low_taxa)
unique(db_mh4$low_taxa)
unique(db_jp5$low_taxa)

names(db_mh3)
names(db_ac2)
names(db_jp4)

length(unique(db_mh3$coreID))
length(unique(db_ac2$coreID))
length(unique(db_jp4$coreID))

DB<-rbind(db_mh3,db_ac2,db_jp4)
str(DB)
unique(DB$site)
unique(DB$year)
unique(DB$month)
unique(DB$depth)
unique(as.character(DB$class))
unique(DB$class1)
unique(DB$family)

DB[is.nan(lenght),lenght:=NA] ## Deal with NaN values
DB[is.nan(lenght),lenght]

write.table(DB,"data_out/db/FinalDB_lenghtinc_20210818.csv",sep=";",row.names = F)




############################### FIQUEI AQUI NA HISTORIA DA BIOMASSA ##################################





test<-aggregate(cbind(numb,AFDW_mg) ~ site+year+month+coreID+depth+class1+family+low_taxa+lenght,data=DB,FUN=sum,na.rm=T)

DB1<-aggregate(DB$numb,by=list(site=DB$site,year=DB$year,month=DB$month,coreID=DB$coreID,depth=DB$depth,
                               class1=DB$class1,family=DB$family,low_taxa=DB$low_taxa,lenght=DB$lenght,AFDW_mg=DB$AFDW_mg),FUN=sum)
str(DB1)
colnames(DB1)[9]<-"numb"
unique(as.character(DB1$low_taxa))
length(unique(DB1$coreID))

#write.table(DB1,"data_out/db/FinalDB_aggregated.csv",sep=";",row.names = F)

##### Complete database with zeros

DB12<-DB1[,c("low_taxa","year","month","coreID","depth","numb","site")]

DB13<-DB1[,c("low_taxa","class1","family")]
DB131<-aggregate(DB13$class1,by=list(low_taxa=DB13$low_taxa,class1=DB13$class1),FUN=length)
DB131<-DB131[,-3]

DB132<-aggregate(DB13$family,by=list(low_taxa=DB13$low_taxa,family=DB13$family),FUN=length)
DB132<-DB132[,-3]


DB12$low_taxa<-as.character(DB12$low_taxa)
length(unique(DB12$coreID))

DB2<-complete(DB12,low_taxa,nesting(year,month,site,coreID),fill = list(numb=0)) # this will leave depth with NA for now, since we are not interested in separating the different depths
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
write.table(DB5,"data_out/db/Final_DB_complete_20201208.csv",row.names=F,sep=";")


DB6<-aggregate(DB5$numb,by=list(year=DB5$year,month=DB5$month,site=DB5$site,coreID=DB5$coreID,class1=DB5$class1,
                                family=DB5$family,low_taxa=DB5$low_taxa),FUN=sum)
colnames(DB6)[8]<-"numb"
DB6$island<-ifelse(DB6$site=="AD","Orango",ifelse(DB6$site=="A"|DB6$site=="AB","Formosa","Bubaque"))

str(DB6)
length(unique(DB6$low_taxa))
length(unique(DB6$coreID))
length(unique(DB6$year))
length(unique(DB6$month))
length(unique(DB6$site))
length(unique(DB6$family))
length(unique(DB6$class1))
table(is.na(DB6))
#write.table(DB6,"data_out/db/Final_DB_complete_density_20201208.csv",row.names=F,sep=";")

## Excluir as linhas das poliquetas com cores de 2020 da base de dados, porque não foram amostradas, apenas os outros grupos
DB66<-DB6[-which(DB6$year==2020 & (DB6$class1=="Polychaeta_sedentaria"|DB6$class1=="Polychaeta_errantia")),]

length(unique(DB66$coreID[DB66$class1=="Polychaeta_sedentaria"]))
length(unique(DB66$coreID[DB66$class1=="Polychaeta_errantia"]))
length(unique(DB66$coreID[DB66$class1=="Bivalvia"]))
length(unique(DB66$coreID[DB66$class1=="Gastropoda"]))

#### Calcular densidades
DB66$dens<-ifelse(DB66$site=="AD",DB66$numb/0.0113,DB66$numb/0.00785)

#BD low_taxa
write.table(DB66,"data_out/db/Final_DB_lowtaxa_density_polyexcl_20210202.csv",row.names=F,sep=";")
DB66<-fread("data_out/db/Final_DB_lowtaxa_density_polyexcl_20210202.csv")

#BD family
DB67<-aggregate(DB66$numb,by=list(island=DB66$island,year=DB66$year,month=DB66$month,site=DB66$site,coreID=DB66$coreID,class1=DB66$class1,family=DB66$family),FUN=sum)
DB67$dens<-ifelse(DB67$site=="AD",DB67$x/0.0113,DB67$x/0.00785)
write.table(DB67,"data_out/db/Final_DB_family_density_polyexcl_20210202.csv",row.names=F,sep=";")

#BD class
DB68<-aggregate(DB66$numb,by=list(island=DB66$island,year=DB66$year,month=DB66$month,site=DB66$site,coreID=DB66$coreID,class1=DB66$class1),FUN=sum)
DB68$dens<-ifelse(DB68$site=="AD",DB68$x/0.0113,DB68$x/0.00785)
write.table(DB68,"data_out/db/Final_DB_class1_density_polyexcl_20210202.csv",row.names=F,sep=";")

#BD geral
DB69<-aggregate(DB66$numb[-which(DB66$class1=="Other")],by=list(island=DB66$island[-which(DB66$class1=="Other")],year=DB66$year[-which(DB66$class1=="Other")],month=DB66$month[-which(DB66$class1=="Other")],site=DB66$site[-which(DB66$class1=="Other")],coreID=DB66$coreID[-which(DB66$class1=="Other")]),FUN=sum)

DB69$dens<-ifelse(DB69$site=="AD",DB69$x/0.0113,DB69$x/0.00785)

write.table(DB69,"data_out/db/Final_DB_total_otherexc_density_polyexcl_20210202.csv",row.names=F,sep=";")
DB69<-fread("data_out/db/Final_DB_total_otherexc_density_polyexcl_20210202.csv")

