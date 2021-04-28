#Title-------------------------------------------------------------------------------
#Data analysis for Ova produced - RFD 2021 paper R. muscosa 

#Packages & Functions-------------------------------------------------------------------------------
library(Rmisc)

#Read in data---------------------------------------------------
setwd("C:/Users/LJacobs/Desktop")
dir()
rfd<-read.csv("rfd_ova_fert_data.csv")

#Clean up data---------------------------------------------------
rfd$f.year<-factor(rfd$year)


#Create table of mean eggs produced per female and NPF
mean.table<- summarySE(rfd, measurevar="ova.produced", groupvars = c("Sex", "X."))

#Visualize mean table for both females and NPF
hist(mean.table[mean.table$Sex=="F",]$ova.produced)
 hist(mean.table[mean.table$Sex=="NPF",]$ova.produced)
 
 
#Run Shapiro-Wild normality test to examine data normalcy
shapiro.test(mean.table[mean.table$Sex=="F",]$ova.produced)

#Run Wilcoxon rank sum test to compare means of females to NPF
wilcox.test(mean.table$ova.produced~mean.table$Sex)
