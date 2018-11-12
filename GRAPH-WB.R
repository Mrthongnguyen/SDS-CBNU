data=read.csv("C:/Users/default.default-PC/Desktop/FINAL DATA/2018-DATA-THONG-FINAL-1.csv")
attach(data)
head(data)
names(data)
summary(data)
data$Group=as.factor(data$Group)
data$Genotype=as.factor(data$Genotype)
data$ID=as.factor(data$Group)
data$GroupxGenotype=as.factor(data$GroupxGenotype)

data=subset(data,GroupxGenotype=="Con-WT"|GroupxGenotype=="Uns-WT"|GroupxGenotype=="Sus-WT")
data$GroupxGenotype <- factor(data$GroupxGenotype, levels=c("Con-WT","Uns-WT","Sus-WT"))

data=subset(data,GroupxGenotype=="Con-KO"|GroupxGenotype=="Uns-KO"|GroupxGenotype=="Sus-KO")
data$GroupxGenotype <- factor(data$GroupxGenotype, levels=c("Con-KO","Uns-KO","Sus-KO"))


data2$GroupxGenotype <- factor(data2$GroupxGenotype, levels=c("Con-WT","Uns-WT","Sus-WT","Con-KO","Uns-KO","Sus-KO"))
data2$Group <- factor(data2$Group, levels=c("Con","Uns","Sus"))
require(reshape2)
data2= melt(data, id= 1:6,measure.vars=c("PFCD2SV","HPD2SV","AMYD2SV","STD2SV"))
data2= melt(data, id= 1:6,measure.vars=c("PFCD2LV","HPD2LV","AMYD2LV","STD2LV"))
data2= melt(data, id= 1:6,measure.vars=c("PFCP34V","HPP34V","AMYP34V","STP34V"))
data2= melt(data, id= 1:6,measure.vars=c("PFCP75V","HPP75V","AMYP75V","STP75V"))
data2= melt(data, id= 1:6,measure.vars=c("PFCD32V","HPD32V","AMYD32V","STD32V"))
data2= melt(data, id= 1:6,measure.vars=c("PFCP16","HPP16","AMYP16","STP16"))
data2= melt(data, id= 1:6,measure.vars=c("PFCS1","HPS1","AMYS1","STS1"))

data2= melt(data, id= 1:6,measure.vars=c("PFCD2S","STD2S","HPD2S","AMYD2S"))
data2= melt(data, id= 1:6,measure.vars=c("PFCD32","STD32","HPD32","AMYD32"))
data2= melt(data, id= 1:6,measure.vars=c("PFCP34","STP34","HPP34","AMYP34"))
data2= melt(data, id= 1:6,measure.vars=c("PFCP75","PFCP75","PFCP75","PFCP75"))




library(ggplot2) 
se <- function(value) sd(value)/sqrt(length(value))   # to calculate standard error in the mean
p = ggplot(data=data2, aes(x=variable, y=value, fill=GroupxGenotype))
p = p + stat_summary(fun.y = mean, geom = "bar", position = "dodge",colour="black")
p = p + labs(title="A")+ylab("D2S/¦Â - Actin (% Control)")+xlab("")+
  scale_x_discrete(name="Regions", 
                   breaks=c("PFCD2SV","HPD2SV","AMYD2SV","STD2SV"), 
                   labels=c("PFC","HP","AMY","ST"))+scale_fill_discrete(name="GroupxGenotype")
p= p + scale_fill_grey() + theme_classic()
p + stat_summary(geom="errorbar",position=position_dodge(width=0.9),
                 fun.data=function(value)c(ymin=mean(value)-se(value),ymax=mean(value)+se(value)), width=0.2)

library(ggplot2) 
se <- function(value) sd(value)/sqrt(length(value))   # to calculate standard error in the mean
p = ggplot(data=data2, aes(x=variable, y=value, fill=GroupxGenotype))
p = p + stat_summary(fun.y = mean, geom = "bar", position = "dodge",colour="black")
p = p + labs(title="B")+ylab("D2L/¦Â - Actin (% Control)")+xlab("")+
  scale_x_discrete(name="Region", 
                   breaks=c("PFCD2LV","HPD2LV","AMYD2LV","STD2LV"), 
                   labels=c("PFC","HP","AMY","ST"))+scale_fill_discrete(name="GroupxGenotype")
p= p + scale_fill_grey() + theme_classic()
p + stat_summary(geom="errorbar",position=position_dodge(width=0.9),
                 fun.data=function(value)c(ymin=mean(value)-se(value),ymax=mean(value)+se(value)), width=0.2)
library(ggplot2) 
se <- function(value) sd(value)/sqrt(length(value))   # to calculate standard error in the mean
p = ggplot(data=data2, aes(x=variable, y=value, fill=GroupxGenotype))
p = p + stat_summary(fun.y = mean, geom = "bar", position = "dodge",colour="black")
p = p + labs(title="C")+ylab("pDARPP-32 Thr34/¦Â - Actin (% Control)")+xlab("")+
  scale_x_discrete(name="Regions", 
                   breaks=c("PFCP34V","HPP34V","AMYP34V","STP34V"), 
                   labels=c("PFC","HP","AMY","ST"))+scale_fill_discrete(name="GroupxGenotype")
p= p + scale_fill_grey() + theme_classic()
p + stat_summary(geom="errorbar",position=position_dodge(width=0.9),
                 fun.data=function(value)c(ymin=mean(value)-se(value),ymax=mean(value)+se(value)), width=0.2)

library(ggplot2) 
se <- function(value) sd(value)/sqrt(length(value))   # to calculate standard error in the mean
p = ggplot(data=data2, aes(x=variable, y=value, fill=GroupxGenotype))
p = p + stat_summary(fun.y = mean, geom = "bar", position = "dodge",colour="black")
p = p + labs(title="D")+ylab("pDARPP-32 Thr75/¦Â - Actin (% Control)")+xlab("")+
  scale_x_discrete(name="Regions", 
                   breaks=c("PFCP75V","HPP75V","AMYP75V","STP75V"), 
                   labels=c("PFC","HP","AMY","ST"))+scale_fill_discrete(name="GroupxGenotype")
p= p + scale_fill_grey() + theme_classic()
p + stat_summary(geom="errorbar",position=position_dodge(width=0.9),
                 fun.data=function(value)c(ymin=mean(value)-se(value),ymax=mean(value)+se(value)), width=0.2)
library(ggplot2) 
se <- function(value) sd(value)/sqrt(length(value))   # to calculate standard error in the mean
p = ggplot(data=data2, aes(x=variable, y=value, fill=GroupxGenotype))
p = p + stat_summary(fun.y = mean, geom = "bar", position = "dodge",colour="black")
p = p + labs(title="E")+ylab("Total DARPP-32/¦Â - Actin (% Control)")+xlab("")+
  scale_x_discrete(name="Regions", 
                   breaks=c("PFCD32V","HPD32V","AMYD32V","STD32V"), 
                   labels=c("PFC","HP","AMY","ST"))+scale_fill_discrete(name="GroupxGenotype")
p= p + scale_fill_grey() + theme_classic()
p + stat_summary(geom="errorbar",position=position_dodge(width=0.9),
                 fun.data=function(value)c(ymin=mean(value)-se(value),ymax=mean(value)+se(value)), width=0.2)

library(ggplot2)

p = ggplot(data=d, aes(x=time, y=resp, group=subject, col=group))

p = p + geom_line() + facet_grid(. ~group) + stat_smooth(aes(group = 1, col="Mean"))

p = p + xlab("Time") + ylab("Resp") + theme(legend.position="none")



p = ggplot(data=data2, aes(x=variable, y=value, group=GroupxGenotype, col=GroupxGenotype))
p = p + geom_line() + facet_grid(. ~GroupxGenotype) + stat_smooth(aes(GroupxGenotype = 1, col="Mean"))
p = p + xlab("Time") + ylab("Resp") + theme(legend.position="none")
p

boxp1=boxp1+annotate("segment", x=c(0.68,0.72,0.93,1.05,1.17,1.32,1.68,1.68,1.93),
                     xend=c(1.68,1.93,1.93), y= c(1350,1400,1400), yend=c(1400,1400,1350))


