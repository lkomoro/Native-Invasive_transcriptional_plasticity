#Intro:####
source('~/Box/Github_repositories/Native-Invasive_transcriptional_plasticity/Scripts/toolbox.functions.R')
setwd('~/Box/Github_repositories/Native-Invasive_transcriptional_plasticity/CTMax_data_&_analyses')
getwd()
library(ggplot2)
library(plyr)
data<-read.csv("ds.menidia.ctmax.csv")
str(data)
which(data$loe.temp.C==32.9)
data<-data[c(1:30,32:36),]#see experiment notes, taking out observation from fish that jumped out and got lost

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

pdf("Mb.vs.Ds.CTMax.bar.pdf",width=6,height=8) #Fig 1B
ggplot(data.sumc, aes(x=species, y=loe.temp.C)) +
  geom_bar(aes(fill=species),color="black", stat="identity")+ scale_fill_manual(values=c("#990000", "blue"))+
  theme_bw()+xlab("Species")+ coord_cartesian(ylim=c(20,38))+ 
  ylab("CTMax (°C)")+geom_errorbar(aes(ymin=loe.temp.C, ymax=loe.temp.C+se), width=.2,size=.5,position=position_dodge(.9))+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14))+
  theme(axis.title.y=element_text(size=18, angle=90))+guides(fill=FALSE)+
  annotate("text",x=2,y=34.5,label="*", size=20)
dev.off()

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
#just CTmax
t.test(loe.temp.C ~ species, data, var.equal=TRUE)

#linear model with size (various versions, all tell same story)
lm<-lm(loe.temp.C ~ species*fork.length.mm, data)
summary(lm)

anova(lm)#if use this output format reads like ancova, need to satisfy assumptions (which we don't)
lm1<-lm(loe.temp.C ~ species*weight.g, data)
summary(lm1)

lm2<-lm(loe.temp.C ~ species*weight.g*fork.length.mm, data)
summary(lm2)

#check model assumptions
data$resid<-resid(lm)
data$fitted<-fitted(lm)
##plot fitted vs. residuals to look for fanning, etc. to check assumption of homogeneity of variances
plot(data$resid~data$fitted)

ggplot(data,aes(x=fitted,y=resid))+geom_point()+theme_bw()

#since species is a factor, also should check equal variance assumption with boxplots of residuals for each group
plot(data$resid~data$species)
ggplot(data,aes(x=species,y=resid))+geom_boxplot()+theme_bw()

##this is a histogram of the residuals, which will tell you if you have normally distributed errors
ggplot(data, aes(x=resid,y=..density..)) + geom_histogram(binwidth=.1)+
  geom_density()
