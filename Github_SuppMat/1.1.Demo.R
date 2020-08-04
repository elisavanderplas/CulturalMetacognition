## Analyse demographics differences in PKU vs. UCL dataset for Table 1 Supplementary Material
## EVDP 2019 elisa.plas.18@ucl.ac.uk

rm(list=ls())
require(R.matlab) 
require(lme4)
require(car)
require(optimx)

bigData = NULL
demoData = NULL
nat = c('PKU', 'UCL')

j=1
for (d in 1:2) {
  if (d == 1) {
    dataset = "PKU"
    dataDir = "~/Dropbox/CulturalMetacognition_2020/DATA/EXP2/PKU_data/PKU_data/"
    filePrefix = "fMRI_pilotData_sub_"
    suffix = "_2"
    subjects = c(403, seq(404,418), seq(420,432), seq(434,435), seq(437,443), seq(445,459))
  } else if ( d == 2) {
    dataset = "UCL"
    dataDir = "~/Dropbox/CulturalMetacognition_2020/DATA/EXP2/UCL_data/UCL_data/"
    filePrefix = "fMRI_pilotData_sub_"
    suffix = "_2"
    subjects =c(seq(25,76), 79) 
  }
  setwd(paste(dataDir,nat[d], '_betas', sep = ""))
  demoData1<- read.csv(paste(nat[d], '_subject_log.csv',sep=""), header = T, sep = ",", na.strings = "EMPTY")
  demoData = rbind(demoData, demoData1)
  }

PKU = demoData[-c(length(subjects)+1:nrow(demoData)),]
UCL = demoData[-c(1:length(subjects)),]

mean(PKU$age, na.rm = T)
sd(PKU$age, na.rm = TRUE)/sqrt(length(PKU$age))
mean(UCL$age, na.rm = T)
t.test(PKU$age, UCL$age, var.equal =T)
mean(PKU$gender, na.rm = T)
mean(UCL$gender, na.rm = T)
t.test(PKU$gender, UCL$gender, var.equal = T)

#relative to PPP difference UK/China
UCL$income = UCL$income/1.71

mean(PKU$income, na.rm = TRUE)
sd(PKU$income, na.rm = TRUE)/sqrt(length(PKU$income))

mean(UCL$income, na.rm = TRUE)
sd(UCL$income, na.rm = TRUE)/sqrt(length(UCL$income))

t.test(PKU$income, UCL$income, var.equal= T)

#IQ
mean(PKU$IQ, na.rm = T)
mean(UCL$IQ, na.rm = T)
t.test(PKU$IQ, UCL$IQ, var.equal = T)

#Holism
mean(PKU$hol_tot, na.rm = T)
mean(UCL$hol_tot, na.rm = T)
t.test(PKU$hol_tot, UCL$hol_tot, var.equal = T)

#Independent construals 
mean(PKU$constr2, na.rm = T)
mean(UCL$constr2, na.rm = T)
t.test(PKU$constr2, UCL$constr2, var.equal = T)

#Dependent construal
mean(PKU$constr1, na.rm = TRUE)
mean(UCL$constr1, na.rm = TRUE)
t.test(PKU$constr1, UCL$constr1, var.equal = T)

#BCIS
mean(PKU$BCIS, na.rm = TRUE)
mean(UCL$BCIS, na.rm = TRUE)
t.test(PKU$BCIS, UCL$BCIS, var.equal =T)

sd(UCL$BCIS, na.rm = TRUE)/sqrt(length(UCL$BCIS))



