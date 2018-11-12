data=read.csv("C:/Users/default.default-PC/Desktop/FINAL DATA/18.10.14/2018-DATA-THONG-FINAL-1.csv")
attach(data)
head(data)
names(data)
summary(data)
#Subset
data=subset(data,GroupxGenotype=="Con-WT"|GroupxGenotype=="Uns-WT"|GroupxGenotype=="Sus-WT")
data=subset(data,GroupxGenotype=="Con-KO"|GroupxGenotype=="Uns-KO"|GroupxGenotype=="Sus-KO")
#reshape
require(reshape2)
data2= melt(data, id= 1:6,measure.vars=c("PFCD2SV","HPD2SV","AMYD2SV","STD2SV"))
data2= melt(data, id= 1:6,measure.vars=c("PFCD2LV","HPD2LV","AMYD2LV","STD2LV"))
data2= melt(data, id= 1:6,measure.vars=c("PFCP34V","HPP34V","AMYP34V","STP34V"))
data2= melt(data, id= 1:6,measure.vars=c("PFCP75V","HPP75V","AMYP75V","STP75V"))
data2= melt(data, id= 1:6,measure.vars=c("PFCD32V","HPD32V","AMYD32V","STD32V"))
data2= melt(data, id= 1:6,measure.vars=c("PFCP16","HPP16","AMYP16","STP16"))
data2= melt(data, id= 1:6,measure.vars=c("PFCS1","HPS1","AMYS1","STS1"))
data2= melt(data, id= 1:6,measure.vars=c("PFCD32","HPD32","AMYD32","STD32"))
data2= melt(data, id= 1:6,measure.vars=c("PFCP34","HPP34","AMYP34","STP34"))
data2= melt(data, id= 1:6,measure.vars=c("PFCP16V","HPP16V","AMYP16V","STP16V"))

data2$Group <- factor(data2$Group, levels=c("Con","Uns","Sus"))

# GRAPH
library(ggplot2) 
se <- function(value) sd(value)/sqrt(length(value))   # to calculate standard error in the mean
p = ggplot(data=data2, aes(x=variable, y=value, fill=Group))
p = p + stat_summary(fun.y = mean, geom = "bar", position = "dodge",colour="black")
p = p + labs(title="(A) WT")+ylab("D2S/GAPDH (% Control)")+xlab("")+
  scale_x_discrete(name="Regions", 
                   breaks=c("PFCD2SV","HPD2SV","AMYD2SV","STD2SV"), 
                   labels=c("PFC","HP","AMY","ST"))+scale_fill_discrete(name="Group")
p= p + scale_fill_grey() + theme_classic()
p= p+ annotate("text",x=1.3,y=1.45,label="*")
p + stat_summary(geom="errorbar",position=position_dodge(width=0.9),
                 fun.data=function(value)c(ymin=mean(value)-se(value),ymax=mean(value)+se(value)), width=0.2)
library(ggplot2) 
se <- function(value) sd(value)/sqrt(length(value))   # to calculate standard error in the mean
p = ggplot(data=data2, aes(x=variable, y=value, fill=Group))
p = p + stat_summary(fun.y = mean, geom = "bar", position = "dodge",colour="black")
p = p + labs(title="(B) KO")+ylab("D2S/GAPDH (% Control)")+xlab("")+
  scale_x_discrete(name="Regions", 
                   breaks=c("PFCD2SV","HPD2SV","AMYD2SV","STD2SV"), 
                   labels=c("PFC","HP","AMY","ST"))+scale_fill_discrete(name="Group")
p= p + scale_fill_grey() + theme_classic()
p= p+ annotate("text",x=1.3,y=0.9,label="*")
p= p+ annotate("text",x=4.3,y=1.45,label="*")
p= p+ annotate("text",x=4.32,y=1.55,label="#")
p + stat_summary(geom="errorbar",position=position_dodge(width=0.9),
                 fun.data=function(value)c(ymin=mean(value)-se(value),ymax=mean(value)+se(value)), width=0.2)
library(ggplot2) 
se <- function(value) sd(value)/sqrt(length(value))   # to calculate standard error in the mean
p = ggplot(data=data2, aes(x=variable, y=value, fill=Group))
p = p + stat_summary(fun.y = mean, geom = "bar", position = "dodge",colour="black")
p = p + labs(title="(A) WT")+ylab("D2L/GAPDH (% Control)")+xlab("")+
  scale_x_discrete(name="Region", 
                   breaks=c("PFCD2LV","HPD2LV","AMYD2LV","STD2LV"), 
                   labels=c("PFC","HP","AMY","ST"))+scale_fill_discrete(name="Group")
p= p + scale_fill_grey() + theme_classic()
p= p+ annotate("text",x=3.3,y=1.6,label="*")
p + stat_summary(geom="errorbar",position=position_dodge(width=0.9),
                 fun.data=function(value)c(ymin=mean(value)-se(value),ymax=mean(value)+se(value)), width=0.2)

library(ggplot2) 
se <- function(value) sd(value)/sqrt(length(value))   # to calculate standard error in the mean
p = ggplot(data=data2, aes(x=variable, y=value, fill=Group))
p = p + stat_summary(fun.y = mean, geom = "bar", position = "dodge",colour="black")
p = p + labs(title="(B) KO")+ylab("D2L/GAPDH (% Control)")+xlab("")+
  scale_x_discrete(name="Region", 
                   breaks=c("PFCD2LV","HPD2LV","AMYD2LV","STD2LV"), 
                   labels=c("PFC","HP","AMY","ST"))+scale_fill_discrete(name="Group")
p= p + scale_fill_grey() + theme_classic()
p + stat_summary(geom="errorbar",position=position_dodge(width=0.9),
                 fun.data=function(value)c(ymin=mean(value)-se(value),ymax=mean(value)+se(value)), width=0.2)

#P34
library(ggplot2) 
se <- function(value) sd(value)/sqrt(length(value))   # to calculate standard error in the mean
p = ggplot(data=data2, aes(x=variable, y=value, fill=Group))
p = p + stat_summary(fun.y = mean, geom = "bar", position = "dodge",colour="black")
p = p + labs(title="(A) WT")+ylab("pDARPP-32 Thr34/¦Â - Actin(% Control)")+xlab("")+
  scale_x_discrete(name="Regions", 
                   breaks=c("PFCP34V","HPP34V","AMYP34V","STP34V"), 
                   labels=c("PFC","HP","AMY","ST"))+scale_fill_discrete(name="Group")
p= p + scale_fill_grey() + theme_classic()
p= p+ annotate("text",x=3.3,y=1.6,label="T")
p= p+ annotate("text",x=3.3,y=1.5,label="$")
p + stat_summary(geom="errorbar",position=position_dodge(width=0.9),
                 fun.data=function(value)c(ymin=mean(value)-se(value),ymax=mean(value)+se(value)), width=0.2)

se <- function(value) sd(value)/sqrt(length(value))   # to calculate standard error in the mean
p = ggplot(data=data2, aes(x=variable, y=value, fill=Group))
p = p + stat_summary(fun.y = mean, geom = "bar", position = "dodge",colour="black")
p = p + labs(title="(C) KO")+ylab("pDARPP-32 Thr34/Total DARPP-32(% Control)")+xlab("")+
  scale_x_discrete(name="Regions", 
                   breaks=c("PFCP34","HPP34","AMYP34","STP34"),  
                   labels=c("PFC","HP","AMY","ST"))+scale_fill_discrete(name="Group")
p= p + scale_fill_grey() + theme_classic()
p= p+ annotate("text",x=3.30,y=0.9,label="**")
p + stat_summary(geom="errorbar",position=position_dodge(width=0.9),
                 fun.data=function(value)c(ymin=mean(value)-se(value),ymax=mean(value)+se(value)), width=0.2)

library(ggplot2) 
se <- function(value) sd(value)/sqrt(length(value))   # to calculate standard error in the mean
p = ggplot(data=data2, aes(x=variable, y=value, fill=Group))
p = p + stat_summary(fun.y = mean, geom = "bar", position = "dodge",colour="black")
p = p + labs(title="(A) WT")+ylab("pDARPP-32 Thr75/¦Â - Actin (% Control)")+xlab("")+
  scale_x_discrete(name="Regions", 
                   breaks=c("PFCP75V","HPP75V","AMYP75V","STP75V"), 
                   labels=c("PFC","HP","AMY","ST"))+scale_fill_discrete(name="Group")
p= p + scale_fill_grey() + theme_classic()
p= p+ annotate("text",x=1.0,y=1.35,label="**")
p= p+ annotate("text",x=1.30,y=1.25,label="T")
p= p+ annotate("text",x=4.0,y=1.45,label="T")
p= p+ annotate("text",x=4.30,y=1.45,label="*")
p + stat_summary(geom="errorbar",position=position_dodge(width=0.9),
                 fun.data=function(value)c(ymin=mean(value)-se(value),ymax=mean(value)+se(value)), width=0.2)


library(ggplot2) 
se <- function(value) sd(value)/sqrt(length(value))   # to calculate standard error in the mean
p = ggplot(data=data2, aes(x=variable, y=value, fill=Group))
p = p + stat_summary(fun.y = mean, geom = "bar", position = "dodge",colour="black")
p = p + labs(title="(B) KO")+ylab("pDARPP-32 Thr75/Total DARPP-32(% Control)")+xlab("")+
  scale_x_discrete(name="Regions", 
                   breaks=c("PFCP75V","HPP75V","AMYP75V","STP75V"), 
                   labels=c("PFC","HP","AMY","ST"))+scale_fill_discrete(name="Group")
p= p + scale_fill_grey() + theme_classic()
p= p+ annotate("text",x=1.0,y=2.0,label="**")
p= p+ annotate("text",x=1.30,y=2.0,label="**")
p + stat_summary(geom="errorbar",position=position_dodge(width=0.9),
                 fun.data=function(value)c(ymin=mean(value)-se(value),ymax=mean(value)+se(value)), width=0.2)
library(ggplot2) 
se <- function(value) sd(value)/sqrt(length(value))   # to calculate standard error in the mean
p = ggplot(data=data2, aes(x=variable, y=value, fill=Group))
p = p + stat_summary(fun.y = mean, geom = "bar", position = "dodge",colour="black")
p = p + labs(title="(B) WT")+ylab("Total DARPP-32/¦Â - Actin (% Control)")+xlab("")+
  scale_x_discrete(name="Regions", 
                   breaks=c("PFCD32","HPD32","AMYD32","STD32"), 
                   labels=c("PFC","HP","AMY","ST"))+scale_fill_discrete(name="Group")
p= p+ annotate("text",x=3.0,y=1.4,label="T")
p= p+ annotate("text",x=3.30,y=1.4,label="*")
p= p+ annotate("text",x=4.30,y=1.3,label="*")
p= p + scale_fill_grey() + theme_classic()
p= p+ annotate("text",x=1.0,y=1.45,label="T")
p + stat_summary(geom="errorbar",position=position_dodge(width=0.9),
                 fun.data=function(value)c(ymin=mean(value)-se(value),ymax=mean(value)+se(value)), width=0.2)

library(ggplot2) 
se <- function(value) sd(value)/sqrt(length(value))   # to calculate standard error in the mean
p = ggplot(data=data2, aes(x=variable, y=value, fill=Group))
p = p + stat_summary(fun.y = mean, geom = "bar", position = "dodge",colour="black")
p = p + labs(title="(C) KO")+ylab("pDARPP-32 Thr75/Total DARPP-32(% Control)")+xlab("")+
  scale_x_discrete(name="Regions", 
                   breaks=c("PFCP75","HPP75","AMYP75","STP75"), 
                   labels=c("PFC","HP","AMY","ST"))+scale_fill_discrete(name="Group")
p= p+ annotate("text",x=3.3,y=0.9,label="*")
p= p + scale_fill_grey() + theme_classic()
p + stat_summary(geom="errorbar",position=position_dodge(width=0.9),
                 fun.data=function(value)c(ymin=mean(value)-se(value),ymax=mean(value)+se(value)), width=0.2)




