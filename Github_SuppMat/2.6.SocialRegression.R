## Analyse covariate differences in PKU vs. UCL dataset 
## EVDP 2019 elisa.plas.18@ucl.ac.uk
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
    dataDir = "~/Dropbox/CulturalMetacognition/DATA/EXP2/PKU_data/PKU_data/"
    filePrefix = "Data_sub_"
    suffix = "_2"
    subjects = c(403, seq(404,418), seq(420,432), seq(434,435), seq(437,443), seq(445,459))
  } else if ( d == 2) {
    dataset = "UCL"
    dataDir = "~/Dropbox/CulturalMetacognition/DATA/Exp2/UCL_data/UCL_data/"
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
      conf_dir = NULL 
      confRT = NULL
      response = NULL
      dir = NULL
      conf_adv = NULL
      task = NULL
      
      precoh = dat[,,1]$dots.coherence[1,]
      postcoh = dat[,,1]$post.coherence[1,]
      coh = sort(unique(precoh))
      
      precoh_cat = ifelse(precoh %in% coh[1], -0.5,
                          ifelse(precoh %in% coh[2], 0, 0.5))
      postcoh_cat = ifelse(postcoh %in% coh[1], -0.5,
                           ifelse(postcoh %in% coh[2], 0, 0.5))
      
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
      conf_adv_cat = quantile(conf_adv_dir, probs = c(0.33, 0.66), na.rm = T, names = F)
      conf_adv_cat = discretize(conf_adv_dir, method = "fixed", breaks = c(0.5, conf_adv_cat[1], conf_adv_cat[2], 1), include.lowest = TRUE, labels = FALSE)
      
      conf_adv_cat <- mapvalues(conf_adv_cat, from = c(1,2,3), to = c(-0.5, 0, 0.5))
      
      RT = dat[,,1]$reaction.time.button[1,]
      conf = dat[,,1]$mouse.response[1,]
      conf_dir = dat[,,1]$mouse.response.dir[1,]
      confRT = dat[,,1]$reaction.time.mouse[1,]
      response = dat[,,1]$button.response[1,] - 1
      response[response == 0] = -1
      dir = dat[,,1]$dots.direction[1,]/360
      dir[dir == 0.5] = -1
      accuracy = response == dir
      
      conf_inv = NA
      for (i in 1:length(precoh_cat)){
        if(dir[i] == -1){
          conf_inv[i] <- (1 - conf_dir[i])}
        else{
          conf_inv[i] <- conf_dir[i]
        }
      }
      accuracy = accuracy -0.5 #to get -0.5 or 0.5
      a_adv = dat[,,1]$a.adv[1,]
      a_adv[a_adv == 99] <- NA
      a_adv = a_adv-1
      agree = a_adv == response
      agree = agree -0.5#to get -0.5 or 0.5
      
      logRT = scale(log(RT))
      logConfRT = scale(log(confRT))
      subj = rep(subjects[s], length(accuracy))
      task = dat[,,1]$condition[1,]
      country = rep(d, length(accuracy))
      
      subData1 = data.frame("country" = country, "subj"=subj, "task"=task, "dir"=dir, "precoh"=precoh, "postcoh"=postcoh,"precoh_cat"=precoh_cat, "postcoh_cat"=postcoh_cat, "conf"=conf, "conf_inv" = conf_inv, "conf_adv" = conf_adv, "conf_adv_dir" = conf_adv_dir, "conf_adv_cat" = conf_adv_cat, "conf_dir" = conf_dir, "logConfRT"=logConfRT, "response"=response, "logRT"=logRT, "accuracy"=accuracy, "a_adv" = a_adv, "agree" = agree)
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

bigData$task <- factor(bigData$task,levels = c(0,1), labels = c("nonsocial", "social"))

## distinguish between correct/incorrect trials
bigData_correct <- bigData[bigData$accuracy == 0.5, ]
bigData_incorrect <- bigData[bigData$accuracy == -0.5, ]

## distinguish between social/nonsocial trials
bigData_social_corr <- bigData_correct[bigData_correct$task == 'social', ]
bigData_nonsocial_corr <- bigData_correct[bigData_correct$task == 'nonsocial', ]
bigData_social_err <- bigData_incorrect[bigData_incorrect$task == 'social', ]
bigData_nonsocial_err <- bigData_incorrect[bigData_incorrect$task == 'nonsocial', ]

## distinguish between social/nonsocial trials, agree vs disagree
bigData_PKU_corr_da <- bigData_social_corr[bigData_social_corr$agree == -0.5 & bigData_social_corr$country == 'PKU', ]
bigData_PKU_corr_a <- bigData_social_corr[bigData_social_corr$agree == 0.5 & bigData_social_corr$country == 'PKU', ]
bigData_PKU_err_da <- bigData_social_err[bigData_social_err$agree == -0.5 & bigData_social_err$country == 'PKU', ]
bigData_PKU_err_a <- bigData_social_err[bigData_social_err$agree == 0.5 & bigData_social_err$country == 'PKU', ]

bigData_UCL_corr_da <- bigData_social_corr[bigData_social_corr$agree == -0.5 & bigData_social_corr$country == 'UCL', ]
bigData_UCL_corr_a <- bigData_social_corr[bigData_social_corr$agree == 0.5 & bigData_social_corr$country == 'UCL', ]
bigData_UCL_err_da <- bigData_social_err[bigData_social_err$agree == -0.5 & bigData_social_err$country == 'UCL', ]
bigData_UCL_err_a <- bigData_social_err[bigData_social_err$agree == 0.5 & bigData_social_err$country == 'UCL', ]

## Conf - agree x accuracy
confModel_social = lmer(conf_inv ~ country*(accuracy + precoh_cat + conf_adv_cat + agree + accuracy:precoh_cat + accuracy:conf_adv_cat + accuracy:agree + precoh_cat:conf_adv_cat + precoh_cat:agree + conf_adv_cat:agree + accuracy:precoh_cat:conf_adv_cat + accuracy:conf_adv_cat:agree + accuracy:precoh_cat:agree + accuracy:precoh_cat:conf_adv_cat:agree + logRT) + (1 + accuracy + precoh_cat + conf_adv_cat + agree + accuracy:precoh_cat + accuracy:conf_adv_cat + accuracy:agree + precoh_cat:conf_adv_cat + precoh_cat:agree + conf_adv_cat:agree + accuracy:precoh_cat:conf_adv_cat + accuracy:conf_adv_cat:agree + accuracy:precoh_cat:agree + accuracy:precoh_cat:conf_adv_cat:agree + logRT|subj), data=bigData_social
                        , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))

fix <- fixef(confModel_social)
fix.se <- sqrt(diag(vcov(confModel_social)))
print(summary(confModel_social))
print(Anova(confModel_social, type = 3))
hist(resid(confModel_social))
coef(summary(confModel_social))

## Get the beta coefficients per country as a function of agreement and accuracy (Supplementary Fig 5B)
# error- agree, PKU
err_a_PKU = lmer(conf_inv ~ precoh_cat*conf_adv_cat + logRT + (1 +  precoh_cat*conf_adv_cat  + logRT|subj), data=bigData_PKU_err_a
                      , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(err_a_PKU)
fix.se <- sqrt(diag(vcov(err_a_PKU)))
betas <- c(fix, fix.se)
setwd('~/Dropbox/CulturalMetacognition/DATA/EXP2/PKU_data/PKU_data/PKU_betas/')
write.csv(betas, file = paste('regression_betas_err_a_PKU.csv'))

# correct-agree PKU 
corr_a_PKU = lmer(conf_inv ~ precoh_cat*conf_adv_cat  + logRT + (1 +  precoh_cat*conf_adv_cat + logRT|subj), data=bigData_PKU_corr_a
                     , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(corr_a_PKU)
fix.se <- sqrt(diag(vcov(corr_a_PKU)))
betas <- c(fix, fix.se)
write.csv(betas, file = paste('regression_betas_corr_a_PKU.csv'))

# error- disagree, PKU
err_da_PKU = lmer(conf_inv ~ precoh_cat*conf_adv_cat + logRT + (1 +  precoh_cat*conf_adv_cat + logRT|subj), data=bigData_PKU_err_da
                 , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(err_da_PKU)
fix.se <- sqrt(diag(vcov(err_da_PKU)))
betas <- c(fix, fix.se)
setwd('~/Dropbox/CulturalMetacognition/DATA/EXP2/PKU_data/PKU_data/PKU_betas/')
write.csv(betas, file = paste('regression_betas_err_da_PKU.csv'))

# correct-disagree PKU 
corr_da_PKU = lmer(conf_inv ~ precoh_cat*conf_adv_cat  + logRT + (1 +  precoh_cat*conf_adv_cat + logRT|subj), data=bigData_PKU_corr_da
                  , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(corr_da_PKU)
fix.se <- sqrt(diag(vcov(corr_da_PKU)))
betas <- c(fix, fix.se)
write.csv(betas, file = paste('regression_betas_corr_da_PKU.csv'))

# error, agree UCL
err_a_UCL = lmer(conf_inv ~ precoh_cat*conf_adv_cat + logRT + (1 +  precoh_cat*conf_adv_cat + logRT|subj), data=bigData_UCL_err_a
                      , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(err_a_UCL)
fix.se <- sqrt(diag(vcov(err_a_UCL)))
betas <- c(fix, fix.se)
setwd('~/Dropbox/CulturalMetacognition/DATA/EXP2/UCL_data/UCL_data/UCL_betas/')
write.csv(betas, file = paste('regression_betas_err_a_UCL.csv'))

# correct, agree UCL
corr_a_UCL = lmer(conf_inv ~ precoh_cat*conf_adv_cat + logRT + (1 +  precoh_cat*conf_adv_cat + logRT|subj), data=bigData_UCL_corr_a
                     , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(corr_a_UCL)
fix.se <- sqrt(diag(vcov(corr_a_UCL)))
betas <- c(fix, fix.se)
write.csv(betas, file = paste('regression_betas_corr_a_UCL.csv'))

# error, disagree UCL
err_da_UCL = lmer(conf_inv ~ precoh_cat*conf_adv_cat + logRT + (1 +  precoh_cat*conf_adv_cat + logRT|subj), data=bigData_UCL_err_da
                 , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(err_da_UCL)
fix.se <- sqrt(diag(vcov(err_da_UCL)))
betas <- c(fix, fix.se)
setwd('~/Dropbox/CulturalMetacognition/DATA/EXP2/UCL_data/UCL_data/UCL_betas/')
write.csv(betas, file = paste('regression_betas_err_da_UCL.csv'))

# correct, agree UCL
corr_da_UCL = lmer(conf_inv ~ precoh_cat*conf_adv_cat + logRT + (1 +  precoh_cat*conf_adv_cat + logRT|subj), data=bigData_UCL_corr_da
                  , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(corr_da_UCL)
fix.se <- sqrt(diag(vcov(corr_da_UCL)))
betas <- c(fix, fix.se)
write.csv(betas, file = paste('regression_betas_corr_da_UCL.csv'))





