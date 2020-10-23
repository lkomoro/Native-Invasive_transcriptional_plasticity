#Intro:####
source('/users/Komo/Desktop/Rfiles/toolbox.functions.R')
setwd('/users/Komo/Desktop/Rfiles')
getwd()
library(ggplot2)
library(plyr)
data<-read.csv("ds.menidia.ctmax.csv")
str(data)
ggplot(data, aes(x=loe.temp.C,y=..density..)) + geom_histogram(fill="cornsilk", colour="grey60",size=.2)+ 
  geom_density()+ylab("Density")+  xlab("CTMax")
qqnorm(data$loe.temp.C)
qqline(data$loe.temp.C)#just adding line
ggplot(data, aes(x=loe.temp.C, fill=species)) + geom_density(alpha=.3)+ ylab("Density")+  xlab("CTMax")
#bimodal pattern/hump from difference between larval fish species

ggplot(data,aes(x=species,y=loe.temp.C))+geom_boxplot()+theme_bw()+xlab("Species")+
  ylab("CTMax (C)")+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14))+
  theme(axis.title.y=element_text(size=18, angle=90))

data.sumc <- summarySE(data, measurevar="loe.temp.C", groupvars=c("species"))

ggplot(data.sumc, aes(x=species, y=loe.temp.C)) +
  geom_bar(aes(fill=species),color="black", stat="identity")+ scale_fill_manual(values=c("#990000", "blue"))+
  theme_bw()+xlab("Species")+ coord_cartesian(ylim=c(20,38))+ 
  ylab("CTMax (°C)")+geom_errorbar(aes(ymin=loe.temp.C, ymax=loe.temp.C+se), width=.2,size=.5,position=position_dodge(.9))+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14))+
  theme(axis.title.y=element_text(size=18, angle=90))+guides(fill=FALSE)+
  annotate("text",x=2,y=34.5,label="*", size=20)

t.test(loe.temp.C ~ species, data, var.equal=TRUE)
str(data)
#just checking size covariate issues, menidia are a bit smaller as expected
ggplot(data,aes(x=species,y=fork.length.mm))+geom_boxplot()+theme_bw()+xlab("Species")+
  ylab("Fork Length (mm)")+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14))+
  theme(axis.title.y=element_text(size=18, angle=90))

#just checking size covariate issues, menidia are a bit smaller as expected
ggplot(data,aes(x=species,y=weight.g))+geom_boxplot()+theme_bw()+xlab("Species")+
  ylab("Weight (g)")+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14))+
  theme(axis.title.y=element_text(size=18, angle=90))

summary(data)
