#working directory and library paths
getwd() 
.libPaths() 
options(stringsAsFactors=FALSE)

#load libraries
library(survival)
library(ppcor)

library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(gridExtra) 

library(limma) #geneSetTest
library(globaltest)
library(GSA) #Efron and Tibshirani (2006)

#######################################################################################
#load metabolite data
#Variables: id, psatest, fastcat, matchid, chol_ab, fastingr, fastingm, fasting1, 
#fasting2, fasting3, timecatm, timecat2, timecat3, timecat4, yearcat2, yearcat3, 
#ageblood, time_bet, advcase, t3b_up, t4_up, dist_at_dx, bmi, lowgrade, higrade,
#higrd43, logrd34, 243 metabolites
metabData <- read.table("~/metabolomics/data/metab.for.ericka.2017.withgrade.csv", header=TRUE, sep=",")
colnames(metabData)

#load class information
classData <- read.table("~/metabolomics/data/metabolite_classes.2017.csv", header=TRUE, sep=",")
colnames(classData)

#load annotation information
annotData <- read.table("~/metabolomics/data/annotation.2017.csv", header=TRUE, sep=",")
colnames(annotData)

########################################################################################
#process metabolite data

#remove observations without matches
mids<-names(which(table(metabData$matchid)==2))
metabData<-metabData[metabData$matchid %in% mids,]
table(metabData$advcase)

#metabolite data
probeinds <- 5:247
data<-metabData[,probeinds]
names(data)

########################################################################################
#identify differentially expressed metabolites by case status using clogit
metabNames<-names(data)
x <- data.frame(data,
                id = metabData$id,
                psatest = metabData$psatest,
                fastcat = factor(metabData$fastcat, levels=c(1,2,3,4,NA), exclude=NULL),
                matchid = metabData$matchid,
                chol_ab = metabData$chol_ab,
                fastingr = metabData$fastingr,
                fastingm = metabData$fastingm,
                fasting1 = metabData$fasting1,
                fasting2 = metabData$fasting2,
                fasting3 = metabData$fasting3,
                timecatm = metabData$timecatm,
                timecat2 = metabData$timecat2,
                timecat3 = metabData$timecat3,
                timecat4 = metabData$timecat4,
                yearcat2 = metabData$yearcat2,
                yearcat3 = metabData$yearcat3,
                ageblood = metabData$ageblood,
                time_bet = metabData$time_bet,
                advcase = metabData$advcase,
                t3b_up = metabData$t3b_up,
                t4_up = metabData$t4_up,
                dist_at_dx = metabData$dist_at_dx,
                bmi = metabData$bmi,
                lowgrade = metabData$lowgrade,
                higrade = metabData$higrade,
                higrd43 = metabData$higrd43, 
                logrd34 = metabData$logrd34)

########################################################################################
#main analysis
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))

for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}
adj.pvals_clogit <- p.adjust(pvals_clogit, method="fdr")
ordered <- order(pvals_clogit)
clogitResults<-cbind(name=metabNames[ordered], estimate=format(estimates_clogit[ordered], digits=2),pValue=format(pvals_clogit[ordered], digits=4), 
                 adjpValue=format(adj.pvals_clogit[ordered], digits=4))  
matchinds <- as.numeric(sapply(clogitResults[,"name"], function(x) which(annotData[,"name"] == x)))
clogitResults <- cbind(clogitResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
clogitResults[1:10,] #top 10

#save results
write.table(clogitResults,file="~/metabolomics/results/diffexp_clogit_041417.csv",sep=',', col.names=NA)

#main analysis adjusted for bmi
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm+bmi"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat+bmi"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}
adj.pvals_clogit <- p.adjust(pvals_clogit, method="fdr")
ordered <- order(pvals_clogit)
clogitResults<-cbind(name=metabNames[ordered], estimate=format(estimates_clogit[ordered], digits=2),pValue=format(pvals_clogit[ordered], digits=4), 
                     adjpValue=format(adj.pvals_clogit[ordered], digits=4))  
matchinds <- as.numeric(sapply(clogitResults[,"name"], function(x) which(annotData[,"name"] == x)))
clogitResults <- cbind(clogitResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
clogitResults[1:10,] #top 10

#save results
write.table(clogitResults,file="~/metabolomics/results/diffexp_clogit_adjbmi_041417.csv",sep=',', col.names=NA)

#main analysis adjusted for cholesterol
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm+chol_ab"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat+chol_ab"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}
adj.pvals_clogit <- p.adjust(pvals_clogit, method="fdr")
ordered <- order(pvals_clogit)
clogitResults<-cbind(name=metabNames[ordered], estimate=format(estimates_clogit[ordered], digits=2),pValue=format(pvals_clogit[ordered], digits=4), 
                     adjpValue=format(adj.pvals_clogit[ordered], digits=4))  
matchinds <- as.numeric(sapply(clogitResults[,"name"], function(x) which(annotData[,"name"] == x)))
clogitResults <- cbind(clogitResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
clogitResults[1:10,] #top 10

#save results
write.table(clogitResults,file="~/metabolomics/results/diffexp_clogit_adjchol_041417.csv",sep=',', col.names=NA)

########################################################################################
#by time to dx
x_early<-x[x$matchid %in% x$matchid[which(x$time_bet<66)],]
x_late<-x[x$matchid %in% x$matchid[which(x$time_bet>=66)],]

#early cases (diagnosed <66 months after blood draw, 106 cases)
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_early)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}
adj.pvals_clogit <- p.adjust(pvals_clogit, method="fdr")
ordered <- order(pvals_clogit)
clogitResults<-cbind(name=metabNames[ordered], estimate=format(estimates_clogit[ordered], digits=2),pValue=format(pvals_clogit[ordered], digits=4), 
                     adjpValue=format(adj.pvals_clogit[ordered], digits=4))  
matchinds <- as.numeric(sapply(clogitResults[,"name"], function(x) which(annotData[,"name"] == x)))
clogitResults <- cbind(clogitResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
clogitResults[1:10,] #top 10

#save results
write.table(clogitResults,file="~/metabolomics/results/diffexp_clogit_early_041417.csv",sep=',', col.names=NA)

#late cases (diagnosed >=66 months after blood draw, 106 cases)
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_late)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}
adj.pvals_clogit <- p.adjust(pvals_clogit, method="fdr")
ordered <- order(pvals_clogit)
clogitResults<-cbind(name=metabNames[ordered], estimate=format(estimates_clogit[ordered], digits=2),pValue=format(pvals_clogit[ordered], digits=4), 
                     adjpValue=format(adj.pvals_clogit[ordered], digits=4))  
matchinds <- as.numeric(sapply(clogitResults[,"name"], function(x) which(annotData[,"name"] == x)))
clogitResults <- cbind(clogitResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
clogitResults[1:10,] #top 10

#save results
write.table(clogitResults,file="~/metabolomics/results/diffexp_clogit_late_041417.csv",sep=',', col.names=NA)

########################################################################################
#early, advanced cases
table(x_early$t3b_up,useNA="ifany") #of 106 early cases, 50 T3b or higher
table(x_early$t4_up,useNA="ifany") #of 106 early cases, 26 T4 or higher
table(x_early$dist_at_dx,useNA="ifany") #of 106 early cases, 38 T4 or higher or PSA>20

#early advanced cases (diagnosed <66 months after blood draw with stage t3b or higher, 50 cases)
x_early_t3b<-x[x$matchid %in% x$matchid[which(x$time_bet<66 & x$t3b_up==1)],]

fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_early_t3b)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}
adj.pvals_clogit <- p.adjust(pvals_clogit, method="fdr")
ordered <- order(pvals_clogit)
clogitResults<-cbind(name=metabNames[ordered], estimate=format(estimates_clogit[ordered], digits=2),pValue=format(pvals_clogit[ordered], digits=4), 
                     adjpValue=format(adj.pvals_clogit[ordered], digits=4))  
matchinds <- as.numeric(sapply(clogitResults[,"name"], function(x) which(annotData[,"name"] == x)))
clogitResults <- cbind(clogitResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
clogitResults[1:10,] #top 10

#save results
write.table(clogitResults,file="~/metabolomics/results/diffexp_clogit_earlyt3b_041417.csv",sep=',', col.names=NA)

#early advanced cases (diagnosed <66 months after blood draw with stage t4 or higher or PSA>20, 38 cases)
x_early_t4psa<-x[x$matchid %in% x$matchid[which(x$time_bet<66 & x$dist_at_dx==1)],]

fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_early_t4psa)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}
adj.pvals_clogit <- p.adjust(pvals_clogit, method="fdr")
ordered <- order(pvals_clogit)
clogitResults<-cbind(name=metabNames[ordered], estimate=format(estimates_clogit[ordered], digits=2),pValue=format(pvals_clogit[ordered], digits=4), 
                     adjpValue=format(adj.pvals_clogit[ordered], digits=4))  
matchinds <- as.numeric(sapply(clogitResults[,"name"], function(x) which(annotData[,"name"] == x)))
clogitResults <- cbind(clogitResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
clogitResults[1:10,] #top 10

#save results
write.table(clogitResults,file="~/metabolomics/results/diffexp_clogit_earlyt4psa_041417.csv",sep=',', col.names=NA)

########################################################################################
#stratified by bmi
x_low <- x[which(x$bmi<25),]
x_high <- x[which(x$bmi>=25),]

#low bmi
fits_glm <- vector("list", length(metabNames))
pvals_glm <- rep(0, length(metabNames))
estimates_glm <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+ fasting1 + fasting2 + fasting3 + fastingm + ageblood + psatest + timecat2 + timecat3 + timecat4 + timecatm + yearcat2 + yearcat3"))
  fits_glm[[i]] <- glm(metabFormula, data = x_low, family=binomial())
  estimates_glm[i] <- summary(fits_glm[[i]])$coef[2,"Estimate"]
  pvals_glm[i] <- summary(fits_glm[[i]])$coef[2,"Pr(>|z|)"]
}
adj.pvals_glm <- p.adjust(pvals_glm, method="fdr")
ordered <- order(pvals_glm)
glmResults<-cbind(name=metabNames[ordered], estimate=format(estimates_glm[ordered], digits=2),pValue=format(pvals_glm[ordered], digits=4), 
                     adjpValue=format(adj.pvals_glm[ordered], digits=4))  
matchinds <- as.numeric(sapply(glmResults[,"name"], function(x) which(annotData[,"name"] == x)))
glmResults <- cbind(glmResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
glmResults[1:10,] #top 10

#save results
write.table(glmResults,file="~/metabolomics/results/diffexp_glm_lowbmi_041417.csv",sep=',', col.names=NA)

#low bmi adjusted for bmi
fits_glm <- vector("list", length(metabNames))
pvals_glm <- rep(0, length(metabNames))
estimates_glm <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+ fasting1 + fasting2 + fasting3 + fastingm + ageblood + psatest + timecat2 + timecat3 + timecat4 + timecatm + yearcat2 + yearcat3 + bmi"))
  fits_glm[[i]] <- glm(metabFormula, data = x_low, family=binomial())
  estimates_glm[i] <- summary(fits_glm[[i]])$coef[2,"Estimate"]
  pvals_glm[i] <- summary(fits_glm[[i]])$coef[2,"Pr(>|z|)"]
}
adj.pvals_glm <- p.adjust(pvals_glm, method="fdr")
ordered <- order(pvals_glm)
glmResults<-cbind(name=metabNames[ordered], estimate=format(estimates_glm[ordered], digits=2),pValue=format(pvals_glm[ordered], digits=4), 
                  adjpValue=format(adj.pvals_glm[ordered], digits=4))  
matchinds <- as.numeric(sapply(glmResults[,"name"], function(x) which(annotData[,"name"] == x)))
glmResults <- cbind(glmResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
glmResults[1:10,] #top 10

#save results
write.table(glmResults,file="~/metabolomics/results/diffexp_glm_lowbmi_adjbmi_041417.csv",sep=',', col.names=NA)

#high bmi
fits_glm <- vector("list", length(metabNames))
pvals_glm <- rep(0, length(metabNames))
estimates_glm <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+ fasting1 + fasting2 + fasting3 + fastingm + ageblood + psatest + timecat2 + timecat3 + timecat4 + timecatm + yearcat2 + yearcat3"))
  fits_glm[[i]] <- glm(metabFormula, data = x_high, family=binomial())
  estimates_glm[i] <- summary(fits_glm[[i]])$coef[2,"Estimate"]
  pvals_glm[i] <- summary(fits_glm[[i]])$coef[2,"Pr(>|z|)"]
}
adj.pvals_glm <- p.adjust(pvals_glm, method="fdr")
ordered <- order(pvals_glm)
glmResults<-cbind(name=metabNames[ordered], estimate=format(estimates_glm[ordered], digits=2),pValue=format(pvals_glm[ordered], digits=4), 
                  adjpValue=format(adj.pvals_glm[ordered], digits=4))  
matchinds <- as.numeric(sapply(glmResults[,"name"], function(x) which(annotData[,"name"] == x)))
glmResults <- cbind(glmResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
glmResults[1:10,] #top 10

#save results
write.table(glmResults,file="~/metabolomics/results/diffexp_glm_highbmi_041417.csv",sep=',', col.names=NA)

#high bmi adjusted for bmi
fits_glm <- vector("list", length(metabNames))
pvals_glm <- rep(0, length(metabNames))
estimates_glm <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+ fasting1 + fasting2 + fasting3 + fastingm + ageblood + psatest + timecat2 + timecat3 + timecat4 + timecatm + yearcat2 + yearcat3 + bmi"))
  fits_glm[[i]] <- glm(metabFormula, data = x_high, family=binomial())
  estimates_glm[i] <- summary(fits_glm[[i]])$coef[2,"Estimate"]
  pvals_glm[i] <- summary(fits_glm[[i]])$coef[2,"Pr(>|z|)"]
}
adj.pvals_glm <- p.adjust(pvals_glm, method="fdr")
ordered <- order(pvals_glm)
glmResults<-cbind(name=metabNames[ordered], estimate=format(estimates_glm[ordered], digits=2),pValue=format(pvals_glm[ordered], digits=4), 
                  adjpValue=format(adj.pvals_glm[ordered], digits=4))  
matchinds <- as.numeric(sapply(glmResults[,"name"], function(x) which(annotData[,"name"] == x)))
glmResults <- cbind(glmResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
glmResults[1:10,] #top 10

#save results
write.table(glmResults,file="~/metabolomics/results/diffexp_glm_highbmi_adjbmi_041417.csv",sep=',', col.names=NA)

########################################################################################
#corr with age, bmi, cholesterol

#correlation with bmi
fits_cor <- vector("list", length(metabNames))
pvals_cor <- rep(0, length(metabNames))
estimates_cor <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  fits_cor[[i]] <- cor.test(x$bmi,data[,i])
  estimates_cor[i] <- fits_cor[[i]]$estimate
  pvals_cor[i] <- fits_cor[[i]]$p.value
}
adj.pvals_cor <- p.adjust(pvals_cor, method="fdr")
ordered <- order(pvals_cor)
corResults<-cbind(name=metabNames[ordered], estimate=format(estimates_cor[ordered], digits=2),pValue=format(pvals_cor[ordered], digits=4), 
                  adjpValue=format(adj.pvals_cor[ordered], digits=4))  
matchinds <- as.numeric(sapply(corResults[,"name"], function(x) which(annotData[,"name"] == x)))
corResults <- cbind(corResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
corResults[1:10,] #top 10

#save results
write.table(corResults,file="~/metabolomics/results/cor_bmi_041417.csv",sep=',', col.names=NA)

#correlation with age
fits_cor <- vector("list", length(metabNames))
pvals_cor <- rep(0, length(metabNames))
estimates_cor <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  fits_cor[[i]] <- cor.test(x$ageblood,data[,i])
  estimates_cor[i] <- fits_cor[[i]]$estimate
  pvals_cor[i] <- fits_cor[[i]]$p.value
}
adj.pvals_cor <- p.adjust(pvals_cor, method="fdr")
ordered <- order(pvals_cor)
corResults<-cbind(name=metabNames[ordered], estimate=format(estimates_cor[ordered], digits=2),pValue=format(pvals_cor[ordered], digits=4), 
                  adjpValue=format(adj.pvals_cor[ordered], digits=4))  
matchinds <- as.numeric(sapply(corResults[,"name"], function(x) which(annotData[,"name"] == x)))
corResults <- cbind(corResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
corResults[1:10,] #top 10

#save results
write.table(corResults,file="~/metabolomics/results/cor_age_041417.csv",sep=',', col.names=NA)

#correlation with cholesterol
fits_cor <- vector("list", length(metabNames))
pvals_cor <- rep(0, length(metabNames))
estimates_cor <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  fits_cor[[i]] <- pcor.test(x$chol_ab[which(x$chol_ab!="NA")],data[which(x$chol_ab!="NA"),i],x$ageblood[which(x$chol_ab!="NA")])
  estimates_cor[i] <- fits_cor[[i]]$estimate
  pvals_cor[i] <- fits_cor[[i]]$p.value
}
adj.pvals_cor <- p.adjust(pvals_cor, method="fdr")
ordered <- order(pvals_cor)
corResults<-cbind(name=metabNames[ordered], estimate=format(estimates_cor[ordered], digits=2),pValue=format(pvals_cor[ordered], digits=4), 
                  adjpValue=format(adj.pvals_cor[ordered], digits=4))  
matchinds <- as.numeric(sapply(corResults[,"name"], function(x) which(annotData[,"name"] == x)))
corResults <- cbind(corResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
corResults[1:10,] #top 10

#save results
write.table(corResults,file="~/metabolomics/results/cor_chol_041417.csv",sep=',', col.names=NA)

########################################################################################
#by grade
table(x$logrd34,useNA="ifany") #of 212 cases, 78 grade 3+4 and below
table(x$higrd43,useNA="ifany") #of 212 cases, 104 grade 4+3 and above

#low grade, advanced cases (grade 3+4 and below, 78 cases)
x_logrd34<-x[x$matchid %in% x$matchid[which(x$logrd34==1)],]

fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_logrd34)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}
adj.pvals_clogit <- p.adjust(pvals_clogit, method="fdr")
ordered <- order(pvals_clogit)
clogitResults<-cbind(name=metabNames[ordered], estimate=format(estimates_clogit[ordered], digits=2),pValue=format(pvals_clogit[ordered], digits=4), 
                     adjpValue=format(adj.pvals_clogit[ordered], digits=4))  
matchinds <- as.numeric(sapply(clogitResults[,"name"], function(x) which(annotData[,"name"] == x)))
clogitResults <- cbind(clogitResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
clogitResults[1:10,] #top 10

#save results
write.table(clogitResults,file="~/metabolomics/results/diffexp_clogit_logrd34_112517.csv",sep=',', col.names=NA)

#high grade, advanced cases (grade 4+3 and above, 104 cases)
x_higrd43<-x[x$matchid %in% x$matchid[which(x$higrd43==1)],]

fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_higrd43)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}
adj.pvals_clogit <- p.adjust(pvals_clogit, method="fdr")
ordered <- order(pvals_clogit)
clogitResults<-cbind(name=metabNames[ordered], estimate=format(estimates_clogit[ordered], digits=2),pValue=format(pvals_clogit[ordered], digits=4), 
                     adjpValue=format(adj.pvals_clogit[ordered], digits=4))  
matchinds <- as.numeric(sapply(clogitResults[,"name"], function(x) which(annotData[,"name"] == x)))
clogitResults <- cbind(clogitResults,annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
clogitResults[1:10,] #top 10

#save results
write.table(clogitResults,file="~/metabolomics/results/diffexp_clogit_higrd43_112517.csv",sep=',', col.names=NA)

########################################################################################
#pheatmap analysis - unsupervised

#specify heatmap options - all data (N=424)
mat <- t(data)
colnames(mat)<-metabData$id

annotation <- data.frame(case=factor(metabData$advcase))
levels(annotation$case) <- c("no", "yes")
names(annotation) <- c("advCase")
rownames(annotation) = metabData$id

samplecols <- c("orange", "green")
names(samplecols) <- c("no", "yes")
ann_colors <- list(advCase = samplecols)

#output heatmap
cols <- colorRampPalette(colors = c("#CA0020","#CA0020","#CA0020","#CA0020",brewer.pal(5,"RdBu"),"#0571B0","#0571B0","#0571B0","#0571B0"))
pheatmap(mat, color = rev(cols(100)), cluster_rows=TRUE, cluster_cols=TRUE, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", show_colnames=FALSE, show_rownames=FALSE, treeheight_row=0, treeheight_col=0,scale="row", annotation=annotation, annotation_colors=ann_colors)

#specify heatmap options - low bmi (N=188)
mat_low <- t(data10[which(metabData$bmi<25),])
colnames(mat_low)<-metabData$id[which(metabData$bmi<25)]

annotation_low <- data.frame(case=factor(metabData$advcase[which(metabData$bmi<25)]))
levels(annotation_low$case) <- c("no", "yes")
names(annotation_low) <- c("advCase")
rownames(annotation_low) = metabData$id[which(metabData$bmi<25)]

samplecols <- c("orange", "green")
names(samplecols) <- c("no", "yes")
ann_colors <- list(advCase = samplecols)

#output heatmap
cols <- colorRampPalette(colors = c("#CA0020","#CA0020","#CA0020","#CA0020",brewer.pal(5,"RdBu"),"#0571B0","#0571B0","#0571B0","#0571B0"))
pheatmap(mat_low, color = rev(cols(100)), cluster_rows=TRUE, cluster_cols=TRUE, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", show_colnames=FALSE, show_rownames=FALSE, treeheight_row=0, treeheight_col=0,scale="row", annotation=annotation_low, annotation_colors=ann_colors)

#specify heatmap options - high bmi (N=222)
mat_high <- t(data10[which(metabData$bmi>=25),])
colnames(mat_high)<-metabData$id[which(metabData$bmi>=25)]

annotation_high <- data.frame(case=factor(metabData$advcase[which(metabData$bmi>=25)]))
levels(annotation_high$case) <- c("no", "yes")
names(annotation_high) <- c("advCase")
rownames(annotation_high) = metabData$id[which(metabData$bmi>=25)]

samplecols <- c("orange", "green")
names(samplecols) <- c("no", "yes")
ann_colors <- list(advCase = samplecols)

#output heatmap
cols <- colorRampPalette(colors = c("#CA0020","#CA0020","#CA0020","#CA0020",brewer.pal(5,"RdBu"),"#0571B0","#0571B0","#0571B0","#0571B0"))
pheatmap(mat_high, color = rev(cols(100)), cluster_rows=TRUE, cluster_cols=TRUE, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", show_colnames=FALSE, show_rownames=FALSE, treeheight_row=0, treeheight_col=0,scale="row", annotation=annotation_high, annotation_colors=ann_colors)

#specify heatmap options - controls (N=212)
mat_controls <- t(data10[which(metabData$advcase==0),])
colnames(mat_controls)<-metabData$id[which(metabData$advcase==0)]

annotation_controls <- data.frame(case=factor(metabData$advcase[which(metabData$advcase==0)]))
levels(annotation_controls$case) <- c("no")
names(annotation_controls) <- c("advCase")
rownames(annotation_controls) = metabData$id[which(metabData$advcase==0)]

samplecols <- c("orange", "green")
names(samplecols) <- c("no", "yes")
ann_colors <- list(advCase = samplecols)

#output heatmap
cols <- colorRampPalette(colors = c("#CA0020","#CA0020","#CA0020","#CA0020",brewer.pal(5,"RdBu"),"#0571B0","#0571B0","#0571B0","#0571B0"))
pheatmap(mat_controls, color = rev(cols(100)), cluster_rows=TRUE, cluster_cols=TRUE, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", show_colnames=FALSE, show_rownames=FALSE, treeheight_row=0, treeheight_col=0,scale="row", annotation=annotation_controls, annotation_colors=ann_colors)

#specify heatmap options - cases (N=212)
mat_cases <- t(data10[which(metabData$advcase==1),])
colnames(mat_cases)<-metabData$id[which(metabData$advcase==1)]

annotation_cases <- data.frame(case=factor(metabData$advcase[which(metabData$advcase==1)]))
levels(annotation_cases$case) <- c("yes")
names(annotation_cases) <- c("advCase")
rownames(annotation_cases) = metabData$id[which(metabData$advcase==1)]

samplecols <- c("orange", "green")
names(samplecols) <- c("no", "yes")
ann_colors <- list(advCase = samplecols)

#output heatmap
cols <- colorRampPalette(colors = c("#CA0020","#CA0020","#CA0020","#CA0020",brewer.pal(5,"RdBu"),"#0571B0","#0571B0","#0571B0","#0571B0"))
pheatmap(mat_cases, color = rev(cols(100)), cluster_rows=TRUE, cluster_cols=TRUE, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", show_colnames=FALSE, show_rownames=FALSE, treeheight_row=0, treeheight_col=0,scale="row", annotation=annotation_cases, annotation_colors=ann_colors)

########################################################################################
#pheatmap analysis - supervised
data10_order<-data10[order(metabData$advcase),]
metabData_order<-metabData[order(metabData$advcase),]

#specify heatmap options - all data, ordered (N=424)
mat_order <- t(data10_order)
colnames(mat_order)<-metabData_order$id

annotation_order <- data.frame(case=factor(metabData_order$advcase))
levels(annotation_order$case) <- c("no", "yes")
names(annotation_order) <- c("advCase")
rownames(annotation_order) = metabData_order$id

samplecols <- c("orange", "green")
names(samplecols) <- c("no", "yes")
ann_colors <- list(advCase = samplecols)

#output heatmap
cols <- colorRampPalette(colors = c("#CA0020","#CA0020","#CA0020","#CA0020",brewer.pal(5,"RdBu"),"#0571B0","#0571B0","#0571B0","#0571B0"))
pheatmap(mat_order, color = rev(cols(100)), cluster_rows=TRUE, cluster_cols=FALSE, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", show_colnames=FALSE, show_rownames=FALSE, treeheight_row=0, treeheight_col=0,scale="row", annotation=annotation_order, annotation_colors=ann_colors)

########################################################################################
#principal component analysis
pca <- prcomp(x = data[complete.cases(data),], center = TRUE, scale. = TRUE) #removes 0 samples with missing values
scores = as.data.frame(pca$x)

#scores plot
p1 <- ggplot(scores, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color=factor(metabData$advcase)[complete.cases(data)]), size=3) +
  scale_color_manual(values = c("red", "blue"), name = "case status") +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  labs(x="PC1", y = "PC2", title = "scores plot of 424 samples using 295 metabolites") +
  theme_bw() +
  theme(legend.position = c(0.85, 0.125))

p2 <- ggplot(scores, aes(x = PC2, y = PC3)) + 
  geom_point(aes(color=factor(metabData$advcase)[complete.cases(data)]), size=3) +
  scale_color_manual(values = c("red", "blue"), name = "case status") +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  labs(x="PC2", y = "PC3", title = "scores plot of 424 samples using 295 metabolites") +
  theme_bw() +
  theme(legend.position = c(0.85, 0.125))

p3 <- ggplot(scores, aes(x = PC3, y = PC4)) + 
  geom_point(aes(color=factor(metabData$advcase)[complete.cases(data)]), size=3) +
  scale_color_manual(values = c("red", "blue"), name = "case status") +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  labs(x="PC3", y = "PC4", title = "scores plot of 424 samples using 295 metabolites") +
  theme_bw() +
  theme(legend.position = c(0.85, 0.125))

p4 <- ggplot(scores, aes(x = PC4, y = PC5)) + 
  geom_point(aes(color=factor(metabData$advcase)[complete.cases(data)]), size=3) +
  scale_color_manual(values = c("red", "blue"), name = "case status") +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  labs(x="PC4", y = "PC5", title = "scores plot of 424 samples using 295 metabolites") +
  theme_bw() +
  theme(legend.position = c(0.85, 0.125))

grid.arrange(p1, p2, p3, p4, nrow=2)

#summary
summary(pca)
dev.off(); plot(pca, type = "l") #plot variance

#clogit analysis
pcaNames<-names(scores)[1:10]
x <- data.frame(scores,
                id = metabData$id,
                psatest = metabData$psatest,
                fastcat = factor(metabData$fastcat, levels=c(1,2,3,4,NA), exclude=NULL),
                matchid = metabData$matchid,
                chol_ab = metabData$chol_ab,
                fastingr = metabData$fastingr,
                fastingm = metabData$fastingm,
                fasting1 = metabData$fasting1,
                fasting2 = metabData$fasting2,
                fasting3 = metabData$fasting3,
                timecatm = metabData$timecatm,
                timecat2 = metabData$timecat2,
                timecat3 = metabData$timecat3,
                timecat4 = metabData$timecat4,
                yearcat2 = metabData$yearcat2,
                yearcat3 = metabData$yearcat3,
                ageblood = metabData$ageblood,
                time_bet = metabData$time_bet,
                advcase = metabData$advcase,
                t3b_up = metabData$t3b_up,
                t4_up = metabData$t4_up,
                dist_at_dx = metabData$dist_at_dx,
                bmi = metabData$bmi,
                lowgrade = metabData$lowgrade,
                higrade = metabData$higrade,
                higrd43 = metabData$higrd43, 
                logrd34 = metabData$logrd34)

#main analysis
fits_clogit <- vector("list", length(pcaNames))
pvals_clogit <- rep(0, length(pcaNames))
estimates_clogit <- rep(0, length(pcaNames))

for (i in 1:length(pcaNames)) {
  #metabFormula<-as.formula(paste("advcase~",pcaNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",pcaNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}
adj.pvals_clogit <- p.adjust(pvals_clogit, method="fdr")
ordered <- order(pvals_clogit)
clogitResults<-cbind(name=pcaNames[ordered], estimate=format(estimates_clogit[ordered], digits=2),pValue=format(pvals_clogit[ordered], digits=4), 
                     adjpValue=format(adj.pvals_clogit[ordered], digits=4))  
clogitResults[1:10,] #top 10

#save results
write.table(clogitResults,file="~/metabolomics/results/pca_clogit_112617.csv",sep=',', col.names=NA)

#PC2
loadings<-pca$rotation
absloadings <- abs(pca$rotation)
perccontr<-sweep(absloadings, 2, colSums(absloadings), "/")
matchinds <- as.numeric(sapply(rownames(perccontr), function(x) which(annotData[,"name"] == x)))
PC2results <- cbind(PC2=perccontr[,2],annotData[matchinds,c("molecule","m_z","method","class_1","class_2","class_3")])
PC2results<-PC2results[order(PC2results[,1],decreasing=TRUE),]

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(rownames(loadings) %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, loadings[,2], alternative="up", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#######################################################################################
#geneSetTest (mean-rank gene set test, Michaud et al (2008))
metabNames<-names(data)
x <- data.frame(data,
                id = metabData$id,
                psatest = metabData$psatest,
                fastcat = factor(metabData$fastcat, levels=c(1,2,3,4,NA), exclude=NULL),
                matchid = metabData$matchid,
                chol_ab = metabData$chol_ab,
                fastingr = metabData$fastingr,
                fastingm = metabData$fastingm,
                fasting1 = metabData$fasting1,
                fasting2 = metabData$fasting2,
                fasting3 = metabData$fasting3,
                timecatm = metabData$timecatm,
                timecat2 = metabData$timecat2,
                timecat3 = metabData$timecat3,
                timecat4 = metabData$timecat4,
                yearcat2 = metabData$yearcat2,
                yearcat3 = metabData$yearcat3,
                ageblood = metabData$ageblood,
                time_bet = metabData$time_bet,
                advcase = metabData$advcase,
                t3b_up = metabData$t3b_up,
                t4_up = metabData$t4_up,
                dist_at_dx = metabData$dist_at_dx,
                bmi = metabData$bmi,
                lowgrade = metabData$lowgrade,
                higrade = metabData$higrade,
                higrd43 = metabData$higrd43, 
                logrd34 = metabData$logrd34)

#main analysis
#get estimates and p-values from clogit
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_clogit), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                  adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#main analysis adjusted for bmi
#get estimates and p-values from clogit
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))

for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm+bmi"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat+bmi"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_clogit), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#main analysis adjusted for cholesterol
#get estimates and p-values from clogit
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm+chol_ab"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat+chol_ab"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_clogit), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#by time to dx
x_early<-x[x$matchid %in% x$matchid[which(x$time_bet<66)],]
x_late<-x[x$matchid %in% x$matchid[which(x$time_bet>=66)],]

#early cases (diagnosed <66 months after blood draw, 106 cases)
#get estimates and p-values from clogit
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_early)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_clogit), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#early cases adjusted for bmi
#get estimates and p-values from clogit
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm+bmi"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat+bmi"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_early)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_clogit), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))


#late cases (diagnosed >=66 months after blood draw, 106 cases)
#get estimates and p-values from clogit
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_late)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_clogit), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#late cases adjusted for bmi
#get estimates and p-values from clogit
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm+bmi"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat+bmi"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_late)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_clogit), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#early, advanced cases
table(x_early$t3b_up,useNA="ifany") #of 106 early cases, 50 T3b or higher
table(x_early$t4_up,useNA="ifany") #of 106 early cases, 26 T4 or higher
table(x_early$dist_at_dx,useNA="ifany") #of 106 early cases, 38 T4 or higher or PSA>20

#early advanced cases (diagnosed <66 months after blood draw with stage t3b or higher, 50 cases)
x_early_t3b<-x[x$matchid %in% x$matchid[which(x$time_bet<66 & x$t3b_up==1)],]

#get estimates and p-values from clogit
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_early_t3b)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_clogit), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#early advanced cases adjusted for bmi
x_early_t3b<-x[x$matchid %in% x$matchid[which(x$time_bet<66 & x$t3b_up==1)],]

#get estimates and p-values from clogit
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm+bmi"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat+bmi"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_early_t3b)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_clogit), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#early advanced cases (diagnosed <66 months after blood draw with stage t4 or higher or PSA>20, 38 cases)
x_early_t4psa<-x[x$matchid %in% x$matchid[which(x$time_bet<66 & x$dist_at_dx==1)],]

#get estimates and p-values from clogit
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_early_t4psa)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_clogit), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#early advanced cases adjusted for bmi
x_early_t4psa<-x[x$matchid %in% x$matchid[which(x$time_bet<66 & x$dist_at_dx==1)],]

#get estimates and p-values from clogit
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm+bmi"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat+bmi"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_early_t4psa)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_clogit), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#stratified by bmi
x_low <- x[which(metabData$bmi<25),]
x_high <- x[which(metabData$bmi>=25),]

#low bmi
#get estimates and p-values from glm
fits_glm <- vector("list", length(metabNames))
pvals_glm <- rep(0, length(metabNames))
estimates_glm <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+ fasting1 + fasting2 + fasting3 + fastingm + ageblood + psatest + timecat2 + timecat3 + timecat4 + timecatm + yearcat2 + yearcat3"))
  fits_glm[[i]] <- glm(metabFormula, data = x_low, family=binomial())
  estimates_glm[i] <- summary(fits_glm[[i]])$coef[2,"Estimate"]
  pvals_glm[i] <- summary(fits_glm[[i]])$coef[2,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_glm), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#low bmi adjusted for bmi
#get estimates and p-values from glm
fits_glm <- vector("list", length(metabNames))
pvals_glm <- rep(0, length(metabNames))
estimates_glm <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+ fasting1 + fasting2 + fasting3 + fastingm + ageblood + psatest + timecat2 + timecat3 + timecat4 + timecatm + yearcat2 + yearcat3 + bmi"))
  fits_glm[[i]] <- glm(metabFormula, data = x_low, family=binomial())
  estimates_glm[i] <- summary(fits_glm[[i]])$coef[2,"Estimate"]
  pvals_glm[i] <- summary(fits_glm[[i]])$coef[2,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_glm), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#high bmi
#get estimates and p-values from glm
fits_glm <- vector("list", length(metabNames))
pvals_glm <- rep(0, length(metabNames))
estimates_glm <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+ fasting1 + fasting2 + fasting3 + fastingm + ageblood + psatest + timecat2 + timecat3 + timecat4 + timecatm + yearcat2 + yearcat3"))
  fits_glm[[i]] <- glm(metabFormula, data = x_high, family=binomial())
  estimates_glm[i] <- summary(fits_glm[[i]])$coef[2,"Estimate"]
  pvals_glm[i] <- summary(fits_glm[[i]])$coef[2,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_glm), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#high bmi adjusted for bmi
#get estimates and p-values from glm
fits_glm <- vector("list", length(metabNames))
pvals_glm <- rep(0, length(metabNames))
estimates_glm <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+ fasting1 + fasting2 + fasting3 + fastingm + ageblood + psatest + timecat2 + timecat3 + timecat4 + timecatm + yearcat2 + yearcat3 + bmi"))
  fits_glm[[i]] <- glm(metabFormula, data = x_high, family=binomial())
  estimates_glm[i] <- summary(fits_glm[[i]])$coef[2,"Estimate"]
  pvals_glm[i] <- summary(fits_glm[[i]])$coef[2,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_glm), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#corr with age, bmi, cholesterol

#correlation with bmi
#get estimates and p-values from cor.test
fits_cor <- vector("list", length(metabNames))
pvals_cor <- rep(0, length(metabNames))
estimates_cor <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  fits_cor[[i]] <- cor.test(x$bmi,data[,i])
  estimates_cor[i] <- fits_cor[[i]]$estimate
  pvals_cor[i] <- fits_cor[[i]]$p.value
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_cor), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#correlation with age
#get estimates and p-values from cor.test
fits_cor <- vector("list", length(metabNames))
pvals_cor <- rep(0, length(metabNames))
estimates_cor <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  fits_cor[[i]] <- cor.test(x$ageblood,data[,i])
  estimates_cor[i] <- fits_cor[[i]]$estimate
  pvals_cor[i] <- fits_cor[[i]]$p.value
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_cor), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#correlation with cholesterol
#get estimates and p-values from pcor.test
fits_cor <- vector("list", length(metabNames))
pvals_cor <- rep(0, length(metabNames))
estimates_cor <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  fits_cor[[i]] <- pcor.test(x$chol_ab[which(x$chol_ab!="NA")],data[which(x$chol_ab!="NA"),i],x$ageblood[which(x$chol_ab!="NA")])
  estimates_cor[i] <- fits_cor[[i]]$estimate
  pvals_cor[i] <- fits_cor[[i]]$p.value
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_cor), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#by grade
table(x$logrd34,useNA="ifany") #of 212 cases, 78 grade 3+4 and below
table(x$higrd43,useNA="ifany") #of 212 cases, 104 grade 4+3 and above

#low grade, advanced cases (grade 3+4 and below, 78 cases)
x_logrd34<-x[x$matchid %in% x$matchid[which(x$logrd34==1)],]

#get estimates and p-values from clogit
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_logrd34)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_clogit), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#low grade advanced cases adjusted for bmi
x_logrd34<-x[x$matchid %in% x$matchid[which(x$logrd34==1)],]

#get estimates and p-values from clogit
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm+bmi"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat+bmi"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_logrd34)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_clogit), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#high grade, advanced cases (grade 4+3 and above, 104 cases)
x_higrd43<-x[x$matchid %in% x$matchid[which(x$higrd43==1)],]

#get estimates and p-values from clogit
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_higrd43)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_clogit), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#high grade advanced cases adjusted for bmi
x_higrd43<-x[x$matchid %in% x$matchid[which(x$higrd43==1)],]

#get estimates and p-values from clogit
fits_clogit <- vector("list", length(metabNames))
pvals_clogit <- rep(0, length(metabNames))
estimates_clogit <- rep(0, length(metabNames))
for (i in 1:length(metabNames)) {
  #metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fasting1+fasting2+fasting3+fastingm+bmi"))
  metabFormula<-as.formula(paste("advcase~",metabNames[i],"+strata(matchid)+fastcat+bmi"))
  fits_clogit[[i]] <- clogit(metabFormula, data = x_higrd43)
  estimates_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"exp(coef)"]
  pvals_clogit[i] <- summary(fits_clogit[[i]])$coef[1,"Pr(>|z|)"]
}

#enrichment analysis
pvals_wilcoxon <- rep(0, dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  inds <- which(metabNames %in% classData[,i])
  pvals_wilcoxon[i] <- geneSetTest(inds, (1-pvals_clogit), alternative="mixed", ranks.only=TRUE)
}
adj.pvals_wilcoxon <- p.adjust(pvals_wilcoxon, method="fdr")
ordered <- order(pvals_wilcoxon)
wilcoxonResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_wilcoxon[ordered], digits=4), 
                       adjpValue=format(adj.pvals_wilcoxon[ordered], digits=4))

#######################################################################################
#globaltest
pvals_gt <- rep(0, dim(classData)[2])
fits_gt <- vector("list", dim(classData)[2])
for (i in 1:dim(classData)[2]) {
  class<-classData[,i]
  x <- subset(data, select=which(names(data) %in% class))
  fits_gt[[i]] <- gt(metabData$advcase, x, model="logistic")
  pvals_gt[i] <- as.numeric(result(fits_gt[[i]])[1])
}
adj.pvals_gt <- p.adjust(pvals_gt, method="fdr")
ordered <- order(pvals_gt)
gtResults<-cbind(name=names(classData)[ordered], pValue=format(pvals_gt[ordered], digits=4), 
      adjpValue=format(adj.pvals_gt[ordered], digits=4))

#maxmean
dataNA = metabData[complete.cases(metabData[,probeinds]),probeinds]
x<-t(dataNA)
y<-(metabData$advcase+1)[complete.cases(metabData[,probeinds])]
genesets<-classData
genenames<-rownames(x)
fit_gsa<-GSA(x,y, genesets, genenames,
    method=c("maxmean"), resp.type=c("Two class unpaired"))

gsaResultsNeg<-GSA.listsets(fit_gsa, geneset.names=names(genesets),FDRcut=1)$negative
gsaResultsPos<-GSA.listsets(fit_gsa, geneset.names=names(genesets),FDRcut=1)$positive

dataNA = metabData[complete.cases(metabData[,probeinds]),probeinds]
x<-t(dataNA)
y<-(metabData$advcase+1)[complete.cases(metabData[,probeinds])]
genesets<-classData
genenames<-rownames(x)
fit_gsa<-GSA(x,y, genesets, genenames,
             method=c("absmean"), resp.type=c("Two class unpaired"))

gsaResultsNeg<-GSA.listsets(fit_gsa, geneset.names=names(genesets),FDRcut=1)$negative
gsaResultsPos<-GSA.listsets(fit_gsa, geneset.names=names(genesets),FDRcut=1)$positive

#######################################################################################
# save data set for summer students
#load metabolite data
#Variables: id, psatest, fastcat, matchid, chol_ab, fastingr, fastingm, fasting1, 
#fasting2, fasting3, timecatm, timecat2, timecat3, timecat4, yearcat2, yearcat3, 
#ageblood, time_bet, advcase, t3b_up, t4_up, dist_at_dx, bmi, lowgrade, higrade,
#higrd43, logrd34, 243 metabolites
metabData <- read.table("~/metabolomics/data/metab.for.ericka.2017.withgrade.csv", header=TRUE, sep=",")
colnames(metabData)

set.seed(423)
metabData$FAKEid<-sample(1:424, 424, replace=FALSE)
metabData$FAKEmatchid<-NA

for (i in 1:424) {
metabData$FAKEmatchid[i]<-metabData$FAKEid[which(metabData$id%in%metabData$matchid[i])]
}

metabData_forstudents<-metabData[,c("HMDB05396","id","FAKEid","psatest","fastcat", "matchid", "FAKEmatchid", "chol_ab", "timecatm",
                                 "timecat2", "timecat3", "timecat4", "yearcat2", "yearcat3", 
                                 "ageblood", "time_bet", "advcase", "t3b_up", "t4_up", "dist_at_dx", "bmi",
                                 "higrd43", "logrd34")]

write.table(metabData_forstudents, file="metabolomics/data/metabData_forstudents.csv",sep=",")
