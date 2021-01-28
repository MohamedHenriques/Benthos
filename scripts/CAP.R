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
data2
str(data2)
data2[,month1:=factor(month, levels=c(10,11,12,1,2,3,4))]
data2[,site1:=factor(site)]

data1log<-log(data1+1)
datasqrt<-data1^.25

###
system.time(
  OM1<-CAPdiscrim(data1log~site1*month1,data=data2,dist="bray",axes=3,m=0,add=F, permutations = 1000)
)
beep()
summary(OM1)

###Add species weights
OM1<-add.spec.scores(OM1,data1log,method="cor.scores",multi=6.5)

###plot cap model on a linear discriminant axis
plot2 <- ordiplot(OM1, type="n",cex=1)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
ordisymbol(plot2, data2, "site1", legend=T,pchs = F,colors = T,col=palette("default"))
ordiellipse(plot2, site1,label=T,col=palette("default"),draw="polygon",lwd=3,alpha=0.3,border=palette("default"))
#ordiellipse(plot2, month1,label=T,col=palette("default"),draw="polygon",lwd=3,alpha=0.3,border=palette("default"))


# plot change in classification success against m
plot(seq(1:14), rep(-1000, 14), xlim=c(1, 14), ylim=c(0, 100), xlab="m", 
     ylab="classification success (percent)", type="n")
for (mseq in 1:14) {
  CAPdiscrim.result <- CAPdiscrim(data1log~site1*month1,data=data2,dist="bray",axes=2,add=F, permutations = 100,m=mseq)
  points(mseq, CAPdiscrim.result$percent)
}



OM2<-capscale(data1log~site1*month1,data=data2,distance = "bray",comm=data1log, metaMDSdist=T)
summary(OM2)
permutest(OM2,permutations=1000)

plot3 <- ordiplot(OM2, type="n")
ordisymbol(plot3, data2, "site1", legend=T,pchs = F,colors = T,col=p)
ordiellipse(plot3, site1,label=T,col=palette("default"),lwd=3)
ordispider(plot3, site1,label=T,col=palette("default"))





ordiplot(OM2, type="text",scaling=1)
cluster <- hclust(distmat, method="single")
ordicluster(plot4, cluster,col="green")


plot4 <- ordiplot(OM2, type="text")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
distdisplayed(data1log, plot4, distx="bray")

anova.cca(OM2, permutations = 1000)

plot4 <- ordiplot(Ordination.model3, type="n")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
ordisymbol(plot4, data2, "site1", legend=T)
ordiellipse(plot4, site1,label=T)

#Plotting categorical environmental variables onto an ordination graph
distmatrix <- vegdist(data1log, method="bray")
Ordination.model3 <- cmdscale(distmatrix, k=nrow(data1log)-1,eig=T, add=F)
plot4 <- ordiplot(Ordination.model3, type="n")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
ordisymbol(plot4, data2, "site1", legend=T)
# Click in the figure where the legend should be placed
attach(data2)
fitted2 <- envfit(plot4, data.frame(site1), permutations=100)
fitted2
plot(fitted2)
plot4 <- ordiplot(Ordination.model3, type="p")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
ordihull(plot4, site1)
plot4 <- ordiplot(Ordination.model3, type="p")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
ordispider(plot4, site1)
plot4 <- ordiplot(Ordination.model3, type="p")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
ordiellipse(plot4, site1)
