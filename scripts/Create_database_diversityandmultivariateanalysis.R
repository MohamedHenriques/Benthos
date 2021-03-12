setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("psych","reshape2","beepr","data.table")
lapply(packs,require,character.only=T)

## Load database 
DB66<-fread("data_out/db/Final_DB_lowtaxa_density_polyexcl_20210202.csv") ### created in script called Database_cleanup_joining
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
db2[,which(sum(numb)==0),by=taxaf] ##Check which species have zero individuals in the database combined outside 2020
db22<-db2[!c(taxaf=="Corbula_cadenati"|taxaf=="Cuspidaria"|taxaf=="Menippe_nodifrons")] ##Remove taxa identified in previous line
db22[,which(sum(numb)==0),by=taxaf] ##Check if it worked

write.table(db22,"Data_out/db/DB_community_analysis_island_20210312.csv",sep=";",row.names=F)

###aggregate and reshape database for analysis
db3<-db22[,lapply(.SD,sum,na.rm=T),.SDcols="numb",by=c("site","month","coreID","taxaf")]
write.table(db3,"Data_out/db/DB_community_analysis_20210312.csv",sep=";",row.names=F)

DB<-dcast.data.table(db3,coreID+site+month~taxaf,value.var="numb")
setkey(DB,coreID,site,month) ## isto define as variaveis core site e month como as variaveis de base para qualquer operação

###Calculate densities
dens1<-function(x){x/0.0113} #for Adonga
dens2<-function(x){x/0.00785} # for the rest of the sites
DB1<-DB[,lapply(.SD,ifelse(site=="AD",dens1,dens2)),by=c("coreID","site","month")]
setkey(DB1,coreID,month,site)

write.csv(DB1,"Data_out/db/DB_multianal_20210312.csv",row.names=F)