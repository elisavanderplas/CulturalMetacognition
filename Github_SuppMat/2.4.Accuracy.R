## EVDP 2019 elisa.plas.18@ucl.ac.uk
## Analyses section 2.4 from Supplementary Material

rm(list=ls())
require(R.matlab) 
require(lme4)
require(car)
require(optimx)
bigData = NULL
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
    # Extract out data from mat files
    for (s in 1:length(subjects)){
      
      ## EXTRACT DATA FILE
      setwd(dataDir)
      DATA = readMat(paste(filePrefix,subjects[s],suffix,'.mat',sep=""))
      
      dat = DATA$locDATA
      
      precoh = NULL
      task = NULL
      pde_level = NULL
      
      precoh = dat[,,1]$dots.coherence[1,]
      postcoh = dat[,,1]$post.coherence[1,]
      task = dat[,,1]$condition[1,]
      coh = sort(unique(precoh))
      
      precoh_cat = ifelse(precoh %in% coh[1], -0.5,
                          ifelse(precoh %in% coh[2], 0, 0.5))
      
      response = dat[,,1]$button.response[1,] - 1
      response[response == 0] = -1
      dir = dat[,,1]$dots.direction[1,]/360
      dir[dir == 0.5] = -1
      accuracy = response == dir
      
      subj = rep(subjects[s], length(accuracy))
      task = dat[,,1]$condition[1,]
      country = rep(d, length(accuracy))
      
      subData1 = data.frame("country" = country, "subj"=subj, "task"=task, "dir"=dir, "precoh"=precoh, "postcoh"=postcoh,"precoh_cat"=precoh_cat,"accuracy"=accuracy)
      bigData = rbind(bigData, subData1)
      
      j=j+1 # total subject counter
    }
}

# Fit regression models
bigData$subj <- factor(bigData$subj)
bigData$country <- factor(bigData$country, levels = c(1,2), labels = c("PKU", "UCL"))
bigData$task <- factor(bigData$task, levels = c(0,1), labels = c("perceptual", "social"))

## FIRST-ORDER PERFORMANCE
accModel_task = glmer(accuracy ~ country*(task*(precoh_cat)) + (1 + country + task + precoh_cat | subj), data=bigData, family="binomial")
fix <- fixef(accModel_task)
fix.se <- sqrt(diag(vcov(accModel_task)))

# P-values and estimates 
print(paste("Accuracy analysis for social vs. nonsocial PDE", sep=""))
print(summary(accModel_task))
print(Anova(accModel_task, type =3))
hist(resid(accModel_task))
coef(summary(accModel_task))




