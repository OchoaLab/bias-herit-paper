


setwd("E:/OneDrive - Duke University/D_disk/5.Duke_PhD/Research/alex_ochoa/1.Heritability/1.Code/Simulation/scripts/g=1")


library(tidyverse)


data<-read_tsv("herit.txt")

dim(data)

dat2<-data[data$kinship %in% c("true","popkin_rom","popkin_mor","std_rom","std_mor"),]

dat2$bias<-dat2$herit_est-dat2$herit

dat2$herit<-as.factor(dat2$herit)

dat2$kinship <- factor(dat2$kinship , levels=c("true","popkin_rom","popkin_mor","std_rom","std_mor"))

dat2$true<-ifelse(dat2$kinship=="popkin_rom" | dat2$kinship=="true", "unbiased","biased" )

dat2$true <- factor(dat2$true , levels=c("unbiased","biased"))


ggplot(dat2, aes(x=herit, y=herit_est, fill=kinship)) + 
  geom_boxplot()+
  labs(x = expression(h[true]^2) , y=expression(h[est]^2) )

ggplot(dat2, aes(x=herit, y=bias, fill=herit)) + 
  geom_boxplot()+
  facet_wrap(~kinship)+
  labs(x = expression(h[true]^2) )+
  theme(legend.position = "none")


ggplot(dat2, aes(x=true, y=herit_est, fill=kinship)) + 
  geom_boxplot()+
  facet_wrap(~herit, scale="free")


