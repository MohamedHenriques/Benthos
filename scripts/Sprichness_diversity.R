setwd("D:/Work/FCUL/Doutoramento/R/Benthos/GitHub/Benthos/Benthos")
graphics.off()
rm(list=ls())

## Pacotes
packs<-c("vegan","ggplot2","viridis","RColorBrewer","xlsx","psych","reshape2","beepr")
lapply(packs,require,character.only=T)


DB66<-read.table("data_out/db/Final_DB_lowtaxa_density_polyexcl_20201208.csv",header=T,sep=";") ### created in script called Database_cleanup_joining
str(DB66)
unique(DB66$low_taxa[which(DB66$class1=="Other")])

db<-DB66[-which(DB66$class1=="Other"),]
unique(db$low_taxa)

db1<-db[-which(db$year==2020),]



###############################################################################################################################
### Checking differences in the mean total numb ind in each site based on boostrapping 71 samples 1000 times
t1<-aggregate(db2$numb,by=list(site=db2$site,coreID=db2$coreID,low_taxa=db2$taxaf),FUN=sum)
t2<-dcast(t1,site+coreID~low_taxa)


set.seed(1)
r=1000
test<-array(NA,c(r,length(unique(t2$site))))
site<-as.character(unique(t2$site))
y=c("A","AB","BI","E","BR","AD")
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
########### Diversidade per site #################
dd<-aggregate(db2$numb,by=list(site=db2$site,coreID=db2$coreID,low_taxa=db2$taxaf),FUN=sum)
dd1<-reshape2::dcast(dd,site+coreID~low_taxa)

R=1000
set.seed(0)
site<-as.character(unique(dd1$site))
y=c("A","AB","BI","E","BR","AD")
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





### Checking differences in the mean total numb ind in each island based on boostrapping 108 samples 1000 times
i1<-aggregate(db1$numb,by=list(island=db1$island,coreID=db1$coreID,low_taxa=db1$low_taxa),FUN=sum)
i2<-dcast(i1,island+coreID~low_taxa)


set.seed(2)
r=1000
test<-array(NA,c(r,length(unique(i2$island))))
island<-as.character(unique(i2$island))
y=c("Orango","Formosa","Bubaque")
island1<-island[order(match(island,y))]
colnames(test)<-island1

for (i in 1:length(island1)) {
  print(island1[i])
  i3<-i2[which(i2$island==island1[i]),]
  for (j in 1:r) {
    s<-i3[sample(1:nrow(i3),108,replace=T),]
    s1<-apply(s[,-c(1:2)],2,sum)
    test[j,i]<-sum(s1)
  }
}
tot1<-data.frame(island=island1,mean_num_ind=apply(test,2,mean),sd_num_ind=apply(test,2,sd),min_num_ind=apply(test,2,min),max_num_ind=apply(test,2,max))



########### Diversidade per island #################
df<-aggregate(db1$numb,by=list(island=db1$island,coreID=db1$coreID,low_taxa=db1$low_taxa),FUN=sum)
df1<-dcast(df,island+coreID~low_taxa)
table(df1$island)

R=1000
set.seed(3)
island<-as.character(unique(df1$island))
y=c("Orango","Formosa","Bubaque")
island1<-island[order(match(island,y))]
g1<-array(NA, c(length(island1), 5))
colnames(g1)<-c("island","Mean_Shannon-Wiener","SD","N_cores","N_ind")

cores<-108
dens1<-function(x){x/0.00866}
dens2<-function(x){x/0.00817}

for (x in 1:length(island1))
{
  print(island1[x])
  g1[x,1]<-island1[x]
  tt<-df1[df1$island==island1[x],-(1:2)]
  tt1<-tt[sample(1:nrow(tt), 1, replace=T),]
  lim<-sum(tt1)
  for (i in 2:cores) {
    print(paste("starting core",i))
    if(lim<=1500){
      tt1[i,]<-tt[sample(1:nrow(tt), 1, replace=T),]
      lim<-sum(tt1)
    } else{
      print(paste("sampling in",island1[x],"reached",lim,"ind at core",i))
    }
  }
  g1[x,4]<-nrow(tt1)
  g1[x,5]<-lim
  print(paste("starting diversity calculations in",island1[x]))
  ### calculating density
  if(island1[x]=="Orango"){
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
  g1[x,2]<-round(mean(rr),2)
  g1[x,3]<-round(sd(rr),2)
}


############################### Now selecting only samples for which we have data for the same months #########################
unique(db1$month)
db2<-db1[which(db1$month=="2"|db1$month=="3"),]
unique(db2$month)

###############################################################################################################################
### Checking differences in the mean total numb ind in each site based on boostrapping 71 samples 1000 times
t1<-aggregate(db2$numb,by=list(site=db2$site,coreID=db2$coreID,low_taxa=db2$low_taxa),FUN=sum)
t2<-dcast(t1,site+coreID~low_taxa)
table(t2$site)

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
    s<-t3[sample(1:nrow(t3),28,replace=T),]
    s1<-apply(s[,-c(1:2)],2,sum)
    test[j,i]<-sum(s1)
  }
}
tot3<-data.frame(site=site1,mean_num_ind=apply(test,2,mean),sd_num_ind=apply(test,2,sd),min_num_ind=apply(test,2,min),max_num_ind=apply(test,2,max))


############## Calcular diversidade Shannon-Wienner limitando o bootstrap a 71 cores e/ou 800 individuos por site
########### Diversidade per site #################
dd<-aggregate(db2$numb,by=list(site=db2$site,coreID=db2$coreID,low_taxa=db2$low_taxa),FUN=sum)
dd1<-dcast(dd,site+coreID~low_taxa)

R=1000
set.seed(0)
site<-as.character(unique(dd1$site))
y=c("AB","A","E","BI","BR","AD")
site1<-site[order(match(site,y))]
d1<-array(NA, c(length(site1), 5))
colnames(d1)<-c("site","Mean_Shannon-Wiener","SD","N_cores","N_ind")

cores<-28
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
    if(lim<=1000){
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

###################### diversity estimates with increasing sample size (core) ###############################

############## Calcular diversidade Shannon-Wienner limitando o bootstrap a 1 até 71 cores e/ou 800 individuos por site
########### Diversidade per site #################
dd<-aggregate(db2$numb,by=list(site=db2$site,coreID=db2$coreID,low_taxa=db2$taxaf),FUN=sum)
dd1<-dcast(dd,site+coreID~low_taxa)
table(dd1$site)

R=1000
set.seed(0)
site<-as.character(unique(dd1$site))
y=c("A","AB","BI","E","BR","AD")
site1<-site[order(match(site,y))]
d1<-array(NA, c(length(site1), 5))
colnames(d1)<-c("site","Mean_Shannon_Wiener","SD","N_cores","N_ind")
myList<-vector("list",107)
cores<-c(2:108)
sum(tt1)
dens1<-function(x){x/0.00866}
dens2<-function(x){x/0.00817}

for (j in cores){
  
  for (x in 1:length(site1))
{
  print(site1[x])
  d1[x,1]<-site1[x]
  tt<-dd1[dd1$site==site1[x],-(1:2)]
  tt1<-tt[sample(1:nrow(tt), 1, replace=T),]
  lim<-sum(tt1)
  for (i in cores[1]:j) {
    print(paste("starting core",i))
    if(lim<=sum(db1$numb[db1$site=="AD"])){
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

myList[[j-1]]<-d1   
}
beep()


A1<-data.frame(site=NA,Mean_Shannon_Wiener=NA,SD=NA,N_cores=NA,N_ind=NA)
names(A1)

for (i in 1:length(myList)){
  A2<-data.frame(myList[[i]])
  #names(A2)<-names(A1)
  A1<-rbind(A1,A2)
}

str(A1)

A1<-A1[-1,]
unique(A1$N_cores[A1$site=="AD"])
A1$site1<-factor(A1$site,levels=c("A","AB","BI","E","BR","AD"))
A1$N_ind<-as.numeric(A1$N_ind)
A1$Mean_Shannon_Wiener<-as.numeric(A1$Mean_Shannon_Wiener)
A1$SD<-as.numeric(A1$SD)
A1$N_cores<-as.numeric(A1$N_cores)

##Plot diversity estimates varying with sample size (core) in each site

ggplot(A1,aes(x=N_cores,y=Mean_Shannon_Wiener,col=site1))+
  #geom_point()+
  #stat_summary(fun=median,geom="line",aes(group=site1))+
  stat_smooth(geom="smooth",method="loess",se=F)+
  theme_bw()+
  scale_colour_manual(values=c("red","steelblue2","royalblue3","darkblue","limegreen","darkgreen"))+
  scale_x_continuous(breaks = seq(2,108,10))

ggplot(A1,aes(x=N_ind,y=Mean_Shannon_Wiener,col=site1))+
  geom_point()+
  stat_smooth(geom="smooth",method="loess",se=F)+
  theme_bw()+
  scale_colour_manual(values=c("red","steelblue2","royalblue3","darkblue","limegreen","darkgreen"))+
  scale_x_continuous(breaks = seq(0,3030,200))





############## Calcular riqueza específica limitando o bootstrap a 1 até 71 cores por site
########### Riqueza específica por site #################
dd<-aggregate(db2$numb,by=list(site=db2$site,coreID=db2$coreID,low_taxa=db2$taxaf),FUN=sum)
dd1<-dcast(dd,site+coreID~low_taxa)
table(dd1$site)



R=1000
set.seed(5)
site<-as.character(unique(dd1$site))
y=c("A","AB","BI","E","BR","AD")
site1<-site[order(match(site,y))]
d1<-array(NA, c(length(site1), 5))
colnames(d1)<-c("site","Mean_Shannon_Wiener","SD","N_cores","N_ind")
myList1<-vector("list",107)
cores<-c(2:108)

presabs<-function(x){ifelse(x>0,1,0)}

for (j in cores){
  
  for (x in 1:length(site1))
  {
    print(site1[x])
    d1[x,1]<-site1[x]
    tt<-dd1[dd1$site==site1[x],-(1:2)]
    tt1<-tt[sample(1:nrow(tt), 1, replace=T),]
    lim<-sum(tt1)
    for (i in cores[1]:j) {
      print(paste("starting core",i))
      if(lim<=sum(db1$numb[db1$site=="AD"])){
        tt1[i,]<-tt[sample(1:nrow(tt), 1, replace=T),]
        lim<-sum(tt1)
      } else{
        print(paste("sampling in",site1[x],"reached",lim,"ind at core",i))
      }
    }
    d1[x,4]<-nrow(tt1)
    d1[x,5]<-lim
    print(paste("starting sp richness calculations in",site1[x]))
    
    rr<-array(NA,R)
    for (y in 1:R)
    {
      sel<-sample(1:nrow(tt1), nrow(tt1), replace=T)
      f<-apply(tt1[sel,],2,sum)
      f1<-sapply(f,presabs)
      rr[y]<-sum(f1)
    }
    d1[x,2]<-round(mean(rr),2)
    d1[x,3]<-round(sd(rr),2)
  }
  
  myList1[[j-1]]<-d1   
}
beep()

myList1[[1]]
B1<-data.frame(site=NA,Mean_Shannon_Wiener=NA,SD=NA,N_cores=NA,N_ind=NA)
names(B1)

for (i in 1:length(myList1)){
  B2<-data.frame(myList1[[i]])
  #names(A2)<-names(A1)
  B1<-rbind(B1,B2)
}

str(B1)
colnames(B1)[2]<-"Mean_sp_richness"

B1<-B1[-1,]
unique(B1$N_cores[B1$site=="AD"])
B1$site1<-factor(B1$site,levels=c("A","AB","BI","E","BR","AD"))
B1$N_ind<-as.numeric(B1$N_ind)
B1$Mean_sp_richness<-as.numeric(B1$Mean_sp_richness)
B1$SD<-as.numeric(B1$SD)
B1$N_cores<-as.numeric(B1$N_cores)

##Plot diversity estimates varying with sample size (core) in each site

ggplot(B1,aes(x=N_cores,y=Mean_sp_richness,col=site1))+
  geom_point()+
  stat_smooth(geom="smooth",method="loess",se=F)+
  theme_bw()+
  labs(x="Number of cores",y="Species richness")+
  scale_colour_manual(values=p)+
  scale_x_continuous(breaks = seq(2,108,10))

ggplot(B1,aes(x=N_ind,y=Mean_sp_richness,col=site1))+
  geom_point()+
  stat_smooth(geom="smooth",method="loess",se=F)+
  theme_bw()+
  labs(x="Number of individuals",y="Species richness")+
  scale_colour_manual(values=p)+
  scale_x_continuous(breaks = seq(0,3030,200))


#### Calculating species richness using 71 cores and a limit of 800 individuals
R=1000
set.seed(6)
site<-as.character(unique(dd1$site))
y=c("A","AB","BI","E","BR","AD")
site1<-site[order(match(site,y))]
s1<-array(NA, c(length(site1), 5))
colnames(s1)<-c("site","Mean_sp_richness","SD","N_cores","N_ind")

cores<-71

presabs<-function(x){ifelse(x>0,1,0)}

for (x in 1:length(site1))
  {
    print(site1[x])
    s1[x,1]<-site1[x]
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
    s1[x,4]<-nrow(tt1)
    s1[x,5]<-lim
    print(paste("starting sp richness calculations in",site1[x]))
    
    rr<-array(NA,R)
    for (y in 1:R)
    {
      sel<-sample(1:nrow(tt1), nrow(tt1), replace=T)
      f<-apply(tt1[sel,],2,sum)
      f1<-sapply(f,presabs)
      rr[y]<-sum(f1)
    }
    s1[x,2]<-round(mean(rr),2)
    s1[x,3]<-round(sd(rr),2)
}
beep()


