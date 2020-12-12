setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("vegan","ggplot2","viridis","RColorBrewer","xlsx","psych","reshape2")
lapply(packs,require,character.only=T)


DB66<-read.table("data_out/db/Final_DB_lowtaxa_density_polyexcl_20201208.csv",header=T,sep=";") ### created in script called Database_cleanup_joining
str(DB66)
unique(DB66$low_taxa[which(DB66$class1=="Other")])

db<-DB66[-which(DB66$class1=="Other"),]
unique(db$low_taxa)

db1<-db[-which(db$year==2020),]



###############################################################################################################################
### Checking differences in the mean total numb ind in each site based on boostrapping 71 samples 1000 times
t1<-aggregate(db1$numb,by=list(site=db1$site,coreID=db1$coreID,low_taxa=db1$low_taxa),FUN=sum)
t2<-dcast(t1,site+coreID~low_taxa)


set.seed(1)
r=1000
test<-array(NA,c(r,length(unique(t2$site))))
site<-as.character(unique(t2$site))
y=c("AB","A","E","BI","BR","AD")
site1<-site[order(match(site,y))]
colnames(test)<-site1

for (i in 1:length(site1)) {
  print(site1[i])
  t3<-t2[which(t2$site==site1[i]),]
  for (j in 1:r) {
    s<-t3[sample(1:nrow(t3),71,replace=T),]
    s1<-apply(s[,-c(1:2)],2,sum)
    test[j,i]<-sum(s1)
  }
}
tot<-data.frame(site=site1,mean_num_ind=apply(test,2,mean),sd_num_ind=apply(test,2,sd),min_num_ind=apply(test,2,min),max_num_ind=apply(test,2,max))
  

############## Calcular diversidade Shannon-Wienner limitando o bootstrap a 71 cores e/ou 800 individuos por site

dd<-aggregate(db1$numb,by=list(site=db1$site,coreID=db1$coreID,low_taxa=db1$low_taxa),FUN=sum)
dd1<-dcast(dd,site+coreID~low_taxa)

R=1000
set.seed(0)
site<-as.character(unique(dd1$site))
y=c("AB","A","E","BI","BR","AD")
site1<-site[order(match(site,y))]
d1<-array(NA, c(length(site1), 5))
colnames(d1)<-c("site","Mean_Shannon-Wiener","SD","N_cores","N_ind")

cores<-71
sum(tt1)
dens1<-function(x){x/0.00866}
dens2<-function(x){x/0.00817}

for (x in 1:length(site1))
{
  print(site1[x])
  d1[x,1]<-site1[x]
  tt<-dd1[dd1$site==site1[x],-(1:2)]
  tt1<-tt[sample(1:nrow(tt), 1, replace=T),]
  lim<-sum(tt1)
  for (i in 2:cores) {
    print(paste("starting core",i))
    if(lim<=800){
    tt1[i,]<-tt[sample(1:nrow(tt), 1, replace=T),]
    lim<-sum(tt1)
    } else{
    print(paste("sampling in",site1[x],"reached",lim,"ind at core",i))
    }
  }
  d1[x,4]<-nrow(tt1)
  d1[x,5]<-lim
  print(paste("starting diversity calculations in",site1[x]))
  ### calculating density
  if(site1[x]=="AD"){
    tt2<-sapply(tt1,dens1)
  } else {
    tt2<-sapply(tt1,dens2)
  }
  rr<-array(NA,R)
  for (y in 1:R)
  {
    sel<-sample(1:nrow(tt2), nrow(tt2), replace=T)
    f<-apply(tt2[sel,], 2, mean)
    f1<-f[-c(which(f==0))]
    rr[y]<-sum(f1/sum(f1)*log(f1/sum(f1)))*-1
  }
  d1[x,2]<-round(mean(rr),2)
  d1[x,3]<-round(sd(rr),2)
}

