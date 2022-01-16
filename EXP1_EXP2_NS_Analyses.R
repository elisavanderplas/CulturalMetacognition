## EVDP 2019 elisa.plas.18@ucl.ac.uk
# This script replicates the Results section of Experiment 2 in chronological order

rm(list=ls())
require(R.matlab) 
require(lme4)
require(car)
require(optimx)

bigData1 = NULL
bigData2 = NULL
demoData = NULL
nat = c('PKU', 'UCL')
j=1
for (d in 1:2) {
  if (d == 1) {
    dataDir = "~/Dropbox/Github/CulturalMetacognition/DATA/EXP1/PKU_data/PKU_data/"
    filePrefix = "Data_sub_"
    suffix = "_2"
    subjects =c(seq(101,109), seq(111,115), seq(117,141))
  } else if (d == 2) {
    dataDir = "~/Dropbox/Github/CulturalMetacognition/DATA/EXP1/UCL_data/UCL_data/"
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
    experiment = rep(1, length(accuracy))
    subData1 = data.frame("country"=country, "subj"=subj, "dir"=dir, "precoh"=precoh, "postcoh"=postcoh, "precoh_cat"=precoh_cat, "postcoh_cat"=postcoh_cat,"conf"=conf, "logConfRT"=logConfRT, "response"=response, "dir"=dir, "logRT"=logRT, "accuracy"=accuracy,"exp" = experiment)
    
    #add to larger file 
    bigData1 = rbind(bigData1, subData1)
    
    j=j+1 # total subject counter
  }
}

j=1
for (d in 1:2) {
  if (d == 1) {
    dataset = "PKU"
    dataDir = "~/Dropbox/Github/CulturalMetacognition/DATA/EXP2/PKU_data/PKU_data/"
    filePrefix = "Data_sub_"
    suffix = "_2"
    subjects = c(403, seq(404,418), seq(420,432), seq(434,435), seq(437,443), seq(445,459))
  } else if ( d == 2) {
    dataset = "UCL"
    dataDir = "~/Dropbox/Github/CulturalMetacognition/DATA/EXP2/UCL_data/UCL_data/"
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
    accuracy[accuracy == 0] <- -1
    condition = dat[,,1]$condition[1,]

    logRT = scale(log(RT))
    logConfRT = scale(log(confRT))
    subj = rep(j, length(accuracy))
    country = rep(d, length(accuracy))
    experiment = rep(2, length(accuracy))
    subData1 = data.frame("country"=country, "subj"=subj, "dir"=dir, "precoh"=precoh, "postcoh"=postcoh, "precoh_cat"=precoh_cat, "postcoh_cat"=postcoh_cat,"conf"=conf, "logConfRT"=logConfRT, "response"=response, "dir"=dir, "logRT"=logRT, "accuracy"=accuracy,"exp" = experiment)
    
    #add to larger file 
    bigData2 = rbind(bigData2, subData1)
    
    j=j+1 # total subject counter
  }
}

bigData <- rbind(bigData1, bigData2)

# Factors
bigData$subj <- factor(bigData$subj)
bigData$country <- factor(bigData$country, levels = c(1,2), labels = c("PKU","UCL"))
bigData$exp <- factor(bigData$exp, levels = c(1,2), labels = c("exp1","exp2"))

## SECOND-ORDER PERFORMANCE
confModel_wcountry = lmer(conf ~ exp*(country*(accuracy + precoh_cat + postcoh_cat + precoh_cat:postcoh_cat + precoh_cat:accuracy + postcoh_cat:accuracy + precoh_cat:postcoh_cat:accuracy + logRT)) + (1 + accuracy + precoh_cat + postcoh_cat + precoh_cat:postcoh_cat + precoh_cat:accuracy + postcoh_cat:accuracy + precoh_cat:postcoh_cat:accuracy + logRT|subj), data=bigData
                          , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE, REML = FALSE)))

fix <- fixef(confModel_wcountry)
print(summary(confModel_wcountry))
print(Anova(confModel_wcountry, type = 3))
coef(summary(confModel_wcountry)) #get the contrast statistics
fix.se <- sqrt(diag(vcov(confModel_wcountry)))

1