#Title-------------------------------------------------------------------------------
#Data analysis for RFD 2021 paper R. muscosa hermaphroditism

#Packages & Functions-------------------------------------------------------------------------------
library(lme4)
library(lmerTest)
library(DHARMa)
library(visreg)
library(ggplot2)
library(gridExtra)

#DHARMa function for model checks
#based on: https://aosmith.rbind.io/2017/12/21/using-dharma-for-residual-checks-of-unsupported-models/
dharma.fun <- function(model, data, responseVar){
  sim_nbz = simulate(model, nsim = 150)
  sim_nbz = do.call(cbind, sim_nbz)
  sim_res_nbz = createDHARMa(simulatedResponse = sim_nbz, 
                             observedResponse = data[,responseVar],
                             fittedPredictedResponse = predict(model),integerResponse = F)
  return(print(testResiduals(sim_res_nbz)))
}

#Read in data---------------------------------------------------
setwd("~/Box/MYLF_CBARP/Manuscripts/Thumbpads Females/hormone R script_data")
d <- read.csv("rfd_hormone_data_for_github.csv")


#log transform hormone values:
d.l <- d
d.l$T <- log10(as.numeric(as.character(d.l$T))+1)
d.l$E3 <- log10(as.numeric(as.character(d.l$E3))+1)
d.l$BE2 <- log10(as.numeric(as.character(d.l$BE2))+1)



#Models -------------------------------------------------------------------------------
d.l$sex <- factor(d.l$sex, levels=c("tpf","f","m"))#re-level factor
d.l$Accession <- factor(d.l$Accession) #store Accession (ID) as factor
d.l$f.month <- factor(d.l$f.month)

mt.fm.i <- (lmer(T ~ sex+f.month + (1|Accession), data=d.l))#
me3.fm.i <- (lmer(E3 ~ sex*f.month + (1|Accession), data=d.l))#
me2.fm.i <- (lmer(BE2 ~ sex*f.month + (1|Accession), data=d.l))#
#we get singularity warnings here because identity (Accession) explains 0 variance
#in our data. While we could run a GLM instead without the random effect to eliminate that warning,
#we chose to leave it as-is  to show that our assumption going in 
#was that identity (Accession) might explain some variation, 
#and that we did test/control for this possibility


#Model summaries
summary(mt.fm.i)
summary(me3.fm.i)
summary(me2.fm.i)

#DHARMa checks
dharma.fun(mt.fm.i, d.l, "T")
dharma.fun(me3.fm.i, d.l, "E3")
dharma.fun(me2.fm.i, d.l, "BE2")

#Visreg plots
visreg(mt.fm.i, "sex", by="f.month",partial=T)
visreg(me3.fm.i, "sex", by="f.month", partial=T)
visreg(me2.fm.i, "sex", by="f.month", partial=T)


#Figure 1 -------------------------------------------------------------------------------
colorVals <- c("salmon","dodgerblue","mediumpurple")
d.l$sex <- factor(d.l$sex, levels=c("f","m","tpf"))

p1 <- ggplot(d.l, aes(factor(month), T, fill=sex))+geom_boxplot(outlier.size=0)+
  geom_point(position=position_jitterdodge(),size=2.5,alpha=0.8)+theme_minimal(base_size=20)+
  scale_fill_manual(values=colorVals)
p2 <- ggplot(d.l, aes(factor(month), BE2, fill=sex))+geom_boxplot(outlier.size=0)+
  geom_point(position=position_jitterdodge(),size=2.5,alpha=0.8)+theme_minimal(base_size=20)+
  scale_fill_manual(values=colorVals)
p3 <- ggplot(d.l, aes(factor(month), E3, fill=sex))+geom_boxplot(outlier.size=0)+
  geom_point(position=position_jitterdodge(),size=2.5,alpha=0.8)+theme_minimal(base_size=20)+
  scale_fill_manual(values=colorVals)
grid.arrange(p1,p2,p3)

