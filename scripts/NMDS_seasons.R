setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("vegan","ggplot2","viridis","RColorBrewer","psych","reshape2","beepr","data.table","BiodiversityR","cluster")
lapply(packs,require,character.only=T)

DB1<-fread("Data_out/db/DB_multianal_20210127.csv")

##Prepare data and remove rows with zeros in all columns

data<-DB1[apply(DB1[,!(1:3)],1,sum)!=0] ### Bray curtis does not work with empty cores. So we have to remove all of them

data1<-data[,!(1:3)] ###NMDS requires a matrix of values only, so we need to remove the aggregating variables

data2<-data[,1:3] ###store aggregating variables to use latter