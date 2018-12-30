data=read.csv("C:/Users/default.default-PC/Desktop/FINAL DATA/18.12.21/2018-DATA-THONG-FINAL-2.csv")
names(data)
summary(data)
data$Group=as.factor(data$Group)
data$Genotype=as.factor(data$Genotype)
data$GroupxGenotype=as.factor(data$GroupxGenotype)
data$GroupxGenotype <- factor(data$GroupxGenotype, levels=c("Con-WT","Uns-WT","Sus-WT","Con-KO","Uns-KO","Sus-KO"))

#library
library(FSA) 
library(DescTools)
require(reshape2)
#summary
Sum = Summarize(IZA~ GroupxGenotype, data=data, digits=3)
Sum$se = Sum$sd / sqrt(Sum$n)
Sum$se = signif(Sum$se, digits=3)
Sum
#Combine 2 column to 1 column
data2= melt(data, id= 1:6,measure.vars=c("DT","DTA"))
#change name column
head(data2)

#three way repeated anova interaction
model = lm(value ~ Group*Genotype*variable, data=data2)
anova(model)
model1<-aov(value ~ Group*Genotype*variable, data=data2)
summary(model1)
TukeyHSD(model1)
TukeyHSD(model1, "Group")
#two way anova interaction
model1<-aov(DT~ Group * Genotype, data=data)
summary(model1)
TukeyHSD(model1)
TukeyHSD(model1, "Group")
#One way anova(6groups)
fm1 <- aov(DT~ GroupxGenotype,data=data)
summary(fm1)
PostHocTest(fm1, method = "hsd")
#One way anova(3groups)
data1=subset(data,GroupxGenotype=="Con-WT"|GroupxGenotype=="Uns-WT"|GroupxGenotype=="Sus-WT")
data1=subset(data,GroupxGenotype=="Con-KO"|GroupxGenotype=="Uns-KO"|GroupxGenotype=="Sus-KO")
fm1 <- aov(DT~ GroupxGenotype,data=data)
summary(fm1)
PostHocTest(fm1, method = "hsd")
#T-test
data1=subset(data,GroupxGenotype=="Con-WT"|GroupxGenotype=="Uns-WT")
t.test(DT ~ GroupxGenotype, data = data1)
#PAIR T-TEST
pairwise.t.test(DT,GroupxGenotype,data=data)
#summary normality test
shapiro.test(data$DT)
#nonparametric test
group_by(data, GroupxGenotype) %>%summarise(count = n(),median = median(SUC, na.rm = TRUE),IQR = IQR(DOC, na.rm = TRUE))
group_by(data, GroupxGenotype) %>%summarise(count = n(),median = median(DOC, na.rm = TRUE),IQR = IQR(DOC, na.rm = TRUE))
group_by(data, GroupxGenotype) %>%summarise(count = n(),median = median(DOFC, na.rm = TRUE),IQR = IQR(DOC, na.rm = TRUE))
group_by(data2, variable) %>%summarise(count = n(),median = median(value, na.rm = TRUE),IQR = IQR(DOC, na.rm = TRUE))
wilcox.test(DT ~ GroupxGenotype, data = data1,exact = FALSE)
kruskal.test(DT ~ GroupxGenotype, data = data)
NemenyiTest(x=data$DT, g=data$GroupxGenotype, dist="tukey")
