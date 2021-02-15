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
    dataDir = "~/Dropbox/CulturalMetacognition-master/DATA/EXP1/PKU_data/PKU_data/"
    filePrefix = "Data_sub_"
    suffix = "_2"
    subjects =c(seq(101,109), seq(111,115), seq(117,141))
  } else if (d == 2) {
    dataDir = "~/Dropbox/CulturalMetacognition-master/DATA/EXP1/UCL_data/UCL_data/"
    filePrefix = "Data_sub_"
    suffix = "_2"
    subjects = c(seq(201,204), seq(206, 227), seq(229, 234), seq(236,242))
  } 
  # Extract out data from mat files
  for (s in 1:length(subjects)){
    
    setwd(dataDir)
    DATA = readMat(paste(filePrefix,subjects[s],suffix,'.mat',sep=""))
    
    dat = DATA$locDATA
    
    precoh = NULL
    postcoh = NULL
    RT = NULL
    conf = NULL
    confRT = NULL
    response = NULL
    dir = NULL
    
    precoh = dat[,,1]$dots.coherence[1,]
    postcoh = dat[,,1]$post.coherence[1,]
    
    coh = sort(unique(precoh))
    
    precoh_cat = ifelse(precoh %in% coh[1], -0.5,
                        ifelse(precoh %in% coh[2], 0, 0.5))
    postcoh_cat = ifelse(postcoh %in% coh[1], -0.5,
                         ifelse(postcoh %in% coh[2], 0, 0.5))
    
    RT = dat[,,1]$reaction.time.button[1,]
    conf = dat[,,1]$mouse.response[1,]
    confRT = dat[,,1]$reaction.time.mouse[1,]
    response = dat[,,1]$button.response[1,] - 1
    response[response == 0] = -1
    dir = dat[,,1]$dots.direction[1,]/360
    dir[dir == 0.5] = -1
    accuracy = response == dir
    accuracy_dv = accuracy #accuracy when used as dependent variable: 0/1
    accuracy[accuracy == 0] <- -1 #accuracy when used as predictor variable: -1/1
    
    logRT = scale(log(RT))
    logConfRT = scale(log(confRT))
    subj = rep(j, length(accuracy))
    country = rep(d, length(accuracy))
    subData1 = data.frame("country"=country, "subj"=subj, "dir"=dir, "precoh"=precoh, "postcoh"=postcoh, "precoh_cat"=precoh_cat, "postcoh_cat"=postcoh_cat,"conf"=conf, "logConfRT"=logConfRT, "response"=response, "dir"=dir, "logRT"=logRT, "accuracy"=accuracy, "accuracy_dv" = accuracy_dv)
    
    #add to larger file 
    bigData = rbind(bigData, subData1)
    
    j=j+1 # total subject counter
  }
  setwd(paste(dataDir,nat[d], '_betas', sep = ""))
  demoData1<- read.csv(paste(nat[d], '_subject_log.csv',sep=""), header = T, sep = ",", na.strings = "EMPTY")
  demoData = rbind(demoData, demoData1)
}

## Experiment 1
PKU = c(demoData[1:length(subjects),])
UCL = c(demoData[length(subjects)+1:nrow(demoData),])

#age
mean(PKU$age, na.rm = T)
sd(PKU$age, na.rm = T)/sqrt(length(subjects))
mean(UCL$age, na.rm = T)
sd(UCL$age, na.rm = T)/sqrt(length(subjects))
t.test(PKU$age, UCL$age, var.equal = T)
#gender
mean(PKU$gender, na.rm = T)
mean(UCL$gender, na.rm = T)
t.test(PKU$gender, UCL$gender, var.equal = T)
#income: remove outlier in UCL
UCL$income[UCL$income == 3000000] <- NA
#relative to PPP difference UK/China
UCL$income = UCL$income/1.71
mean(PKU$income, na.rm = T)
sd(PKU$income, na.rm = T)/sqrt(length(subjects))
mean(UCL$income, na.rm = T)
sd(UCL$income, na.rm = T)/sqrt(length(subjects))
t.test(PKU$income, UCL$income, var.equal = T)
#IQ
mean(PKU$IQ, na.rm = T)
sd(PKU$IQ, na.rm = T)/sqrt(length(subjects))
mean(UCL$IQ, na.rm = T)
sd(UCL$IQ, na.rm = T)/sqrt(length(subjects))
t.test(PKU$IQ, UCL$IQ, var.equal = T)
#AHS
mean(PKU$hol_tot, na.rm = T) ##CHECK THIS
sd(PKU$hol_tot, na.rm = T)/sqrt(length(subjects))
mean(UCL$hol_tot, na.rm = T)
sd(UCL$hol_tot, na.rm = T)/sqrt(length(subjects))
t.test(PKU$hol_tot, UCL$hol_tot,var.equal = T)
#SCS independent
mean(PKU$constr2, na.rm = T)
sd(PKU$constr2, na.rm = T)/sqrt(length(subjects))
mean(UCL$constr2, na.rm = T)
sd(UCL$constr2, na.rm = T)/sqrt(length(subjects))
t.test(PKU$constr2, UCL$constr2, var.equal = T)
#SCS inter-dependent
mean(PKU$constr1, na.rm = T)
sd(PKU$constr1, na.rm = T)/sqrt(length(subjects))
mean(UCL$constr1, na.rm = T)
sd(UCL$constr1, na.rm = T)/sqrt(length(subjects))
t.test(PKU$constr1, UCL$constr1, var.equal = T)


##Experiment 2
rm(list=ls())
bigData = NULL
demoData = NULL
nat = c('PKU', 'UCL')

j=1
for (d in 1:2) {
  if (d == 1) {
    dataset = "PKU"
    dataDir = "~/Dropbox/CulturalMetacognition-master/DATA/EXP2/PKU_data/PKU_data/"
    filePrefix = "Data_sub_"
    suffix = "_2"
    subjects = c(403, seq(404,418), seq(420,432), seq(434,435), seq(437,443), seq(445,459))
  } else if ( d == 2) {
    dataset = "UCL"
    dataDir = "~/Dropbox/CulturalMetacognition-master/DATA/EXP2/UCL_data/UCL_data/"
    filePrefix = "Data_sub_"
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
sd(UCL$age, na.rm = TRUE)/sqrt(length(UCL$age))
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
sd(PKU$IQ, na.rm = TRUE)/sqrt(length(PKU$IQ))
mean(UCL$IQ, na.rm = T)
sd(UCL$IQ, na.rm = TRUE)/sqrt(length(UCL$IQ))
t.test(PKU$IQ, UCL$IQ, var.equal = T)

#Holism
mean(PKU$hol_tot, na.rm = T)
sd(PKU$hol_tot, na.rm = TRUE)/sqrt(length(PKU$hol_tot))
mean(UCL$hol_tot, na.rm = T)
sd(UCL$hol_tot, na.rm = TRUE)/sqrt(length(UCL$hol_tot))
t.test(PKU$hol_tot, UCL$hol_tot, var.equal = T)

#Independent construals 
mean(PKU$constr2, na.rm = T)
sd(PKU$constr2, na.rm = TRUE)/sqrt(length(PKU$constr2))
mean(UCL$constr2, na.rm = T)
sd(UCL$constr2, na.rm = TRUE)/sqrt(length(UCL$constr2))
t.test(PKU$constr2, UCL$constr2, var.equal = T)

#Dependent construal
mean(PKU$constr1, na.rm = TRUE)
sd(PKU$constr1, na.rm = TRUE)/sqrt(length(PKU$constr1))
mean(UCL$constr1, na.rm = TRUE)
sd(UCL$constr1, na.rm = TRUE)/sqrt(length(UCL$constr1))
t.test(PKU$constr1, UCL$constr1, var.equal = T)

#BCIS
mean(PKU$BCIS, na.rm = TRUE)
sd(PKU$BCIS, na.rm = TRUE)/sqrt(length(PKU$BCIS))
mean(UCL$BCIS, na.rm = TRUE)
sd(UCL$BCIS, na.rm = TRUE)/sqrt(length(UCL$BCIS))
t.test(PKU$BCIS, UCL$BCIS, var.equal =T)

#BCIS Self-reflective
mean(PKU$BCIS_SR, na.rm = TRUE)
sd(PKU$BCIS_SR, na.rm = TRUE)/sqrt(length(PKU$BCIS_SR))
mean(UCL$BCIS_SR, na.rm = TRUE)
sd(UCL$BCIS_SR, na.rm = TRUE)/sqrt(length(UCL$BCIS_SR))
t.test(PKU$BCIS_SR, UCL$BCIS_SR, var.equal =T)

#BCIS Self-confidence
mean(PKU$BCIS_SC, na.rm = TRUE)
sd(PKU$BCIS_SC, na.rm = TRUE)/sqrt(length(PKU$BCIS_SC))
mean(UCL$BCIS_SC, na.rm = TRUE)
sd(UCL$BCIS_SC, na.rm = TRUE)/sqrt(length(UCL$BCIS_SC))
t.test(PKU$BCIS_SC, UCL$BCIS_SC, var.equal =T)



