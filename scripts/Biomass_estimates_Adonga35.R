setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("xlsx","reshape2")
lapply(packs,require,character.only=T)

db_jp<-read.xlsx("D:/Work/FCUL/Doutoramento/Capitulos/Spatial_temporal_variability_invertebrates/Data_Paulino/Invert_Adonga_2019_JPaulino_corrected_20201108.xlsx",1)
db_jp1<-db_jp[which(db_jp$Tipo =="Nuca"),]

###Biomass
B_mg<-aggregate(db_jp1$AFDW_mg, by=list(lowest_taxa=db_jp1$Lowest_taxa),FUN=mean)
B_mg0<-aggregate(db_jp1$AFDW_mg_m2, by=list(lowest_taxa=db_jp1$Lowest_taxa),FUN=mean)
B_mg1<-aggregate(db_jp1$n,by=list(lowest_taxa=db_jp1$Lowest_taxa),FUN=sum)
B_mg2<-merge(B_mg,B_mg0, by="lowest_taxa")
names(B_mg2)[2:3]<-c("mean_AFDM_mg","mean_AFDM_mg_m2")
B_mg3<-merge(B_mg2,B_mg1, by="lowest_taxa")
names(B_mg3)[4]<-"numb_ind"

#write.xlsx(B_mg3,"data_out/db/Adonga_Biomass_estimates_35cores.xlsx",row.names=F)