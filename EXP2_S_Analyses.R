## EVDP 2019 elisa.plas.18@ucl.ac.uk

rm(list=ls())
require(R.matlab) 
require(lme4)
require(car)
require(optimx)
require(arules) #used to discretize conf adviser on line 94
require(plyr)#used to mapvalues on line 96

bigData = NULL
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
      conf_dir = NULL 
      confRT = NULL
      response = NULL
      dir = NULL
      conf_adv = NULL
      task = NULL
      pde_level = NULL
      
      precoh = dat[,,1]$dots.coherence[1,]
      postcoh = dat[,,1]$post.coherence[1,]
      coh = sort(unique(precoh))
      
      precoh_cat = ifelse(precoh %in% coh[1], -0.5,
                          ifelse(precoh %in% coh[2], 0, 0.5))
      postcoh_cat = ifelse(postcoh %in% coh[1], -0.5,
                           ifelse(postcoh %in% coh[2], 0, 0.5))
      
      RT = dat[,,1]$reaction.time.button[1,]
      conf = dat[,,1]$mouse.response[1,]
      conf_dir = dat[,,1]$mouse.response.dir[1,]
      confRT = dat[,,1]$reaction.time.mouse[1,]
      response = dat[,,1]$button.response[1,] - 1
      response[response == 0] = -1
      dir = dat[,,1]$dots.direction[1,]/360
      dir[dir == 0.5] = -1
      accuracy = response == dir
      accuracy[accuracy == 0] <- -1 
      logRT = scale(log(RT))
      logConfRT = scale(log(confRT))
      task = dat[,,1]$condition[1,]
      conf_adv = dat[,,1]$conf.adv[1,]
      
      conf_adv_dir = NA
      for(i in 1:length(precoh_cat)){
        if(conf_adv[i] < 0.5){
        conf_adv_dir[i] <- (1 - conf_adv[i])} 
        else{
          conf_adv_dir[i] <- conf_adv[i]
        }
      }
      conf_adv_dir[conf_adv_dir == 99] <- NA
      
      # get each ID's binned conf
      conf_adv_cat = quantile(conf_adv_dir, probs = c(0.33, 0.66), na.rm = T, names = F)
      conf_adv_cat = discretize(conf_adv_dir, method = "fixed", breaks = c(0.5, conf_adv_cat[1], conf_adv_cat[2], 1), include.lowest = TRUE, labels = FALSE)
      
      conf_adv_cat <- mapvalues(conf_adv_cat, from = c(1,2,3), to = c(-0.5, 0, 0.5))
      
      conf_inv = NA
      for (i in 1:length(precoh_cat)){
        if(dir[i] == -1){
          conf_inv[i] <- (1 - conf_dir[i])}
        else{
          conf_inv[i] <- conf_dir[i]
        }
      }
      
      #get PDE-level dependent on the task condition
      for(i in 1:length(postcoh_cat)){
        if(task[i] == 1){
          pde_level[i] <- (conf_adv_cat[i])} 
        else{
          pde_level[i] <- postcoh_cat[i]
        }
      }
     
      a_adv = dat[,,1]$a.adv[1,]
      a_adv[a_adv == 99] <- NA
      a_adv = a_adv-1    # get from 1/2 to 0/1
      agree = a_adv == response
      agree = agree -0.5#to get -0.5 or 0.5
      
      subj = rep(subjects[s], length(accuracy))
      country = rep(d, length(accuracy))
      
      subData1 = data.frame("country" = country, "subj"=subj, "task"=task, "dir"=dir, "precoh"=precoh, "postcoh"=postcoh,"precoh_cat"=precoh_cat, "postcoh_cat"=postcoh_cat, "conf"=conf, "conf_inv" = conf_inv, "conf_adv" = conf_adv, "conf_adv_dir" = conf_adv_dir, "conf_adv_cat" = conf_adv_cat, "conf_dir" = conf_dir, "logConfRT"=logConfRT, "response"=response, "logRT"=logRT, "accuracy"=accuracy, "a_adv" = a_adv, "agree" = agree, "pde_level"= pde_level)
      bigData = rbind(bigData, subData1)
      
      j=j+1 # total subject counter
    }
}

# Fit regression models
bigData$subj <- factor(bigData$subj)
bigData$country <- factor(bigData$country, levels = c(1,2), labels = c("PKU", "UCL"))

## distinguish between social/nonsocial trials
bigData_social <- bigData[bigData$task == 1, ]
bigData_nonsocial <- bigData[bigData$task == 0, ]

## distinguish between error/correct trials
bigData_err <- bigData_social[bigData_social$accuracy == -1, ]
bigData_corr <- bigData_social[bigData_social$accuracy == 1, ]

## distinguish between PKU/UCL trials
bigData_PKU <- bigData_social[bigData_social$country == "PKU",]
bigData_UCL <- bigData_social[bigData_social$country == "UCL",]

## distinguish between both PKU/UCL trials and correct/error trials on non-social condition
bigData_PKU_correct <- bigData_social[bigData_social$country == "PKU" & bigData_social$accuracy == 1,]
bigData_UCL_correct <- bigData_social[bigData_social$country == "UCL" & bigData_social$accuracy == 1,]
bigData_PKU_incorrect <- bigData_social[bigData_social$country == "PKU" & bigData_social$accuracy == -1,]
bigData_UCL_incorrect <- bigData_social[bigData_social$country == "UCL" & bigData_social$accuracy == -1,]

#Main Model with pre*post interaction - across countries 
confModel = lmer(conf ~ country*(accuracy + precoh_cat + conf_adv_cat + precoh_cat:conf_adv_cat + precoh_cat:accuracy + conf_adv_cat:accuracy + precoh_cat:conf_adv_cat:accuracy + logRT) + (1 + accuracy + precoh_cat + conf_adv_cat + precoh_cat:conf_adv_cat + precoh_cat:accuracy + conf_adv_cat:accuracy + precoh_cat:conf_adv_cat:accuracy + logRT|subj), data=bigData_social
                          , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE, REML = FALSE)))
print(summary(confModel))
print(Anova(confModel, type = 3))
hist(resid(confModel))
coef(summary(confModel))

## Get the beta coefficients per country (Fig 3A plotted)
# error, PKU
error_coef_PKU = lmer(conf ~ precoh_cat*postcoh_cat + logRT + (1 +  precoh_cat*postcoh_cat + logRT|subj), data=bigData_PKU_incorrect
                      , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(error_coef_PKU)
fix.se <- sqrt(diag(vcov(error_coef_PKU)))
betas <- c(fix, fix.se)
setwd('~/Dropbox/CulturalMetacognition_2020/DATA/EXP2/PKU_data/PKU_data/PKU_betas/')
write.csv(betas, file = paste('regression_betas_EXP2_S_err_PKU.csv'))

# correct, PKU 
corr_coef_PKU = lmer(conf ~ precoh_cat*postcoh_cat + logRT + (1 +  precoh_cat*postcoh_cat + logRT|subj), data=bigData_PKU_correct
                     , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(corr_coef_PKU)
fix.se <- sqrt(diag(vcov(corr_coef_PKU)))
betas <- c(fix, fix.se)
write.csv(betas, file = paste('regression_betas_EXP2_S_corr_PKU.csv'))

# error, UCL
error_coef_UCL = lmer(conf ~ precoh_cat*postcoh_cat + logRT + (1 +  precoh_cat*postcoh_cat + logRT|subj), data=bigData_UCL_incorrect
                      , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(error_coef_UCL)
fix.se <- sqrt(diag(vcov(error_coef_UCL)))
betas <- c(fix, fix.se)
setwd('~/Dropbox/CulturalMetacognition_2020/DATA/EXP2/UCL_data/UCL_data/UCL_betas/')
write.csv(betas, file = paste('regression_betas_EXP2_S_err_UCL.csv'))

# correct, UCL
corr_coef_UCL = lmer(conf ~ precoh_cat*postcoh_cat + logRT + (1 +  precoh_cat*postcoh_cat + logRT|subj), data=bigData_UCL_correct
                     , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(corr_coef_UCL)
fix.se <- sqrt(diag(vcov(corr_coef_UCL)))
betas <- c(fix, fix.se)
write.csv(betas, file = paste('regression_betas_EXP2_S_corr_UCL.csv'))


## get the statistical significance for the country-interaction on error or correct trials (Fig. 3B)
confModel_error = lmer(conf ~ country * (precoh_cat*conf_adv_cat + logRT) + (1 +  precoh_cat*conf_adv_cat + logRT|subj), data=bigData_err
                       , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
print(summary(confModel_error))
print(Anova(confModel_error, type =3))

confModel_correct = lmer(conf ~ country * (precoh_cat*conf_adv_cat + logRT) + (1 +  precoh_cat*conf_adv_cat + logRT|subj), data=bigData_corr
                         , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
print(summary(confModel_correct))
print(Anova(confModel_correct, type =3))


## Analyse comparison between condition and countries (line 845 - 850 manuscript) 
confModel_sitetask = lmer(conf ~ country*(task*(accuracy + precoh_cat + pde_level + accuracy:precoh_cat + accuracy:pde_level + precoh_cat:pde_level + accuracy:precoh_cat:pde_level + logRT)) + (1 + accuracy + precoh_cat + pde_level + accuracy:precoh_cat + accuracy:pde_level + precoh_cat:pde_level + accuracy:precoh_cat:pde_level + logRT|subj), data=bigData
                          , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
print(summary(confModel_sitetask))
print(Anova(confModel_sitetask, type = 3))
hist(resid(confModel_sitetask))
coef(summary(confModel_sitetask))

## Analyse the social condition across both countries (line 851-856 manuscript)
confModel_social = lmer(conf ~ accuracy + precoh_cat + conf_adv_cat + precoh_cat:conf_adv_cat + precoh_cat:accuracy + conf_adv_cat:accuracy + precoh_cat:conf_adv_cat:accuracy + logRT + (1 + accuracy + precoh_cat + conf_adv_cat + precoh_cat:conf_adv_cat + precoh_cat:accuracy + conf_adv_cat:accuracy + precoh_cat:conf_adv_cat:accuracy + logRT|subj), data=bigData_social
                 , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE, REML = FALSE)))
print(summary(confModel_social))
print(Anova(confModel_social, type = 3))
hist(resid(confModel_social))
coef(summary(confModel_social))




