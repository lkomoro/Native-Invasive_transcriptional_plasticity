library(ggplot2)

g<-read.table("~/Box/Github_repositories/Native-Invasive_transcriptional_plasticity/DDIG_data&analyses/Transdecoder&Orthofinderwork/Orthogroup_DE_analyses/Temptime_splitanalysisorthogroups_summary_sig_genes_LMK.txt",header=T)
s<-subset(g,Direction!="NotSig")
s$Time_point<-factor(s$Time_point)
s$Direction<-factor(s$Direction, levels = c("Up","Down"))
s1<-aggregate(Freq ~ Treat_Combined, s, sum)
s2<-s1 %>% separate(Treat_Combined, c("Temp_treat", "Time_point","Species","Temp_rounded"))
s2$Time_point<-factor(s2$Time_point)

#Temp_x-axis:
all.together<-ggplot(s2, aes(x=Temp_rounded, y=Freq,fill=Species, shape=Time_point)) +
  geom_point(size=4)  +  scale_shape_manual(values=c(24,21)) + scale_color_manual(values=c("red","blue"))+theme_bw()

s$Time_spp<-paste0(s$Species,s$Time_point)
Up_down_pos<-ggplot(s,aes(x=Temp_rounded,y=Freq, fill=Species, shape=Time_point))+
  geom_point(size=4,alpha=.8,position = position_dodge(.5))+geom_line()+theme_bw()+geom_line(aes(group =Time_spp, color=Species))+
  scale_shape_manual(values=c(24,21)) + 
  scale_fill_manual(values=c("red4","blue"))+
  scale_color_manual(values=c("red4","blue"))+
  facet_grid(Direction~.) 

pdf("Up_down_pos.pdf", 7, 8)
Up_down_pos
dev.off()

Up_down_neg<-ggplot(s,aes(x=Temp_rounded,y=Freq, fill=Species, shape=Time_point))+
  geom_point(size=4,alpha=.8,position = position_dodge(.5))+geom_line()+theme_bw()+geom_line(aes(group =Time_spp, color=Species))+
  scale_y_reverse()+
  scale_shape_manual(values=c(24,21)) + 
  scale_fill_manual(values=c("red4","blue"))+
  scale_color_manual(values=c("red4","blue"))+
  facet_grid(Direction~.)

pdf("Up_down_neg.pdf", 7, 8)
Up_down_neg
dev.off()

#Relative treatment x axis:
s$Temp_treat<-factor(s$Temp_treat, levels = c("HC","CTMax12","CTMax8","CTMax6","CTMax4","CTMax2"))
ggplot(s,aes(x=Temp_treat,y=Freq, color=Species,group=Time_point))+geom_point()+facet_grid(Direction~Time_point)
