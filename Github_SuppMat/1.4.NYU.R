## NYU vs. PKU vs. UCL comparison (Supplementary Material)
rm(list=ls())
require(R.matlab) 
require(lme4)
require(car)
require(optimx)
options(contrasts = c("contr.treatment", "contr.poly")) # This is R defaults but set it anyway to be safe

bigData = NULL
j=1
for (d in 1:3) {
  if (d == 1) {
    dataDir = "~/Dropbox/CulturalMetacognition/DATA/EXP1/PKU_data/PKU_data/"
    filePrefix = "Data_sub_"
    suffix = "_2"
    subjects =c(seq(101,109), seq(111,115), seq(117,141))
  } else if ( d == 2) {
    dataDir = "~/Dropbox/CulturalMetacognition/DATA/EXP1/NYU_data/NYU_data/"
    filePrefix = "Data_sub_"
    suffix = "_2"
    subjects = c(seq(12,28), seq(30,37))
  } else if (d == 3) {
    dataDir = "~/Dropbox/CulturalMetacognition/DATA/EXP1/UCL_data/UCL_data/"
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
    precoh_cat = NULL
    postcoh_cat = NULL
    precoh = dat[,,1]$dots.coherence[1,]
    postcoh = dat[,,1]$post.coherence[1,]
    RT = dat[,,1]$reaction.time.button[1,]
    conf = dat[,,1]$mouse.response[1,]
    confRT = dat[,,1]$reaction.time.mouse[1,]
    response = dat[,,1]$button.response[1,] - 1
    response[response == 0] = -1
    dir = dat[,,1]$dots.direction[1,]/360
    dir[dir == 0.5] = -1
    accuracy = response == dir
    accuracy[accuracy == 0] <- -1
    
    # get on every trial a dummy-coded coherence level
    coh = sort(unique(precoh))
    
    precoh_cat = ifelse(precoh %in% coh[1], -0.5,
                        ifelse(precoh %in% coh[2], 0, 0.5))
    postcoh_cat = ifelse(postcoh %in% coh[1], -0.5,
                         ifelse(postcoh %in% coh[2], 0, 0.5))
    
    logRT = scale(log(RT))
    logConfRT = scale(log(confRT))
    subj = rep(j, length(accuracy))
    country = rep(d, length(accuracy))
    subData = data.frame("country"=country, "subj"=subj, "dir"=dir, "precoh"=precoh, "postcoh"=postcoh, "conf"=conf, "logConfRT"=logConfRT, "response"=response, "dir"=dir, "logRT"=logRT, "accuracy"=accuracy, "precoh_cat"=precoh_cat, "postcoh_cat"=postcoh_cat)
    
    bigData = rbind(bigData, subData)
    
    j=j+1 # total subject counter
  }
}

# Fit regression models
bigData$subj <- factor(bigData$subj)
bigData$country <- factor(bigData$country, levels = c(1,2,3), labels = c("PKU", "NYU", "UCL"))
bigData_PKU <- bigData[bigData$country == "PKU",]
bigData_UCL <- bigData[bigData$country == "UCL",]
bigData_NYU <- bigData[bigData$country == "NYU",]

## distinguish between correct/incorrect & PKU/UCL trials
bigData_correct <- bigData[bigData$accuracy == 1, ]
bigData_incorrect <- bigData[bigData$accuracy == -1, ] 
bigData_NYU_correct <- bigData_correct[bigData_correct$country == "NYU",]
bigData_NYU_incorrect <- bigData_incorrect[bigData_incorrect$country == "NYU",]

#Main replication in the new datasets (Supplementary Material 1.4)
confModel_PKU = lmer(conf ~ accuracy + precoh_cat + postcoh_cat + precoh_cat:accuracy + postcoh_cat:accuracy + precoh_cat:postcoh_cat + precoh_cat:postcoh_cat:accuracy + logRT + (1 + accuracy + precoh_cat + postcoh_cat + precoh_cat:accuracy + postcoh_cat:accuracy + precoh_cat:postcoh_cat + precoh_cat:postcoh_cat:accuracy + logRT|subj), data=bigData_PKU
                 , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
print(paste("Confidence analysis for PKU", sep=""))
print(summary(confModel_PKU))
print(Anova(confModel_PKU, type = 3))
hist(resid(confModel_PKU))
coef(summary(confModel_PKU))

fix <- fixef(confModel_PKU)
fix.se <- sqrt(diag(vcov(confModel_PKU)))
betas <- c(fix, fix.se)
setwd("~/Dropbox/CulturalMetacognition/DATA/EXP1/PKU_data/PKU_data/PKU_betas/")
write.csv(betas, file = paste('regression_betas_PKU.csv'))

confModel_UCL = lmer(conf ~ accuracy + precoh_cat + postcoh_cat + precoh_cat:accuracy + postcoh_cat:accuracy + precoh_cat:postcoh_cat + precoh_cat:postcoh_cat:accuracy + logRT + (1 + accuracy + precoh_cat + postcoh_cat + precoh_cat:accuracy + postcoh_cat:accuracy + precoh_cat:postcoh_cat + precoh_cat:postcoh_cat:accuracy + logRT|subj), data=bigData_UCL
                     , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
print(paste("Confidence analysis for UCL", sep=""))
print(summary(confModel_UCL))
print(Anova(confModel_UCL, type = 3))
hist(resid(confModel_UCL))
coef(summary(confModel_UCL))

fix <- fixef(confModel_UCL)
fix.se <- sqrt(diag(vcov(confModel_UCL)))
betas <- c(fix, fix.se)
setwd("~/Dropbox/CulturalMetacognition/DATA/EXP1/UCL_data/UCL_data/UCL_betas/")
write.csv(betas, file = paste('regression_betas_UCL.csv'))

#Main replication in the new datasets (Supplementary Material 1.4)
confModel_NYU = lmer(conf ~ accuracy + precoh_cat + postcoh_cat + precoh_cat:accuracy + postcoh_cat:accuracy + precoh_cat:postcoh_cat + precoh_cat:postcoh_cat:accuracy + logRT + (1 + accuracy + precoh_cat + postcoh_cat + precoh_cat:accuracy + postcoh_cat:accuracy + precoh_cat:postcoh_cat + precoh_cat:postcoh_cat:accuracy + logRT|subj), data=bigData_NYU
                     , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(confModel_NYU)
fix.se <- sqrt(diag(vcov(confModel_NYU)))
betas <- c(fix, fix.se)
setwd("~/Dropbox/CulturalMetacognition/DATA/EXP1/NYU_data/NYU_data/NYU_betas/")
write.csv(betas, file = paste('regression_betas_NYU.csv'))

#contrasts between 3 countries with PKU as baseline (Supplementary Material 1.5)
confModel = lmer(conf ~ country*(accuracy + precoh_cat + postcoh_cat + precoh_cat:accuracy + postcoh_cat:accuracy + precoh_cat:postcoh_cat + precoh_cat:postcoh_cat:accuracy + logRT) + (1 + accuracy + precoh_cat + postcoh_cat + precoh_cat:accuracy + postcoh_cat:accuracy + precoh_cat:postcoh_cat + precoh_cat:postcoh_cat:accuracy + logRT|subj), data=bigData
                 , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
print(paste("Confidence analysis across three countries", sep=""))
print(summary(confModel))
print(Anova(confModel, type = 3))
hist(resid(confModel))
coef(summary(confModel))

# Get the coefficients for Supplementary Fig. 2B

# error, NYU
error_coef_NYU = lmer(conf ~ precoh_cat*postcoh_cat + logRT + (1 +  precoh_cat*postcoh_cat + logRT|subj), data=bigData_NYU_incorrect
                      , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(error_coef_NYU)
fix.se <- sqrt(diag(vcov(error_coef_NYU)))
betas <- c(fix, fix.se)
setwd("~/Dropbox/CulturalMetacognition/DATA/EXP1/NYU_data/NYU_data/NYU_betas/")
write.csv(betas, file = paste('regression_betas_err_NYU.csv'))

# correct, NYU
corr_coef_NYU = lmer(conf ~ precoh_cat*postcoh_cat + logRT + (1 +  precoh_cat*postcoh_cat + logRT|subj), data=bigData_NYU_correct
                     , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(corr_coef_NYU)
fix.se <- sqrt(diag(vcov(corr_coef_NYU)))
betas <- c(fix, fix.se)
write.csv(betas, file = paste('regression_betas_corr_NYU.csv'))

## Country interaction on error trials (Suppfig 2C/D)
confModel_error = lmer(conf ~ country * (precoh_cat*postcoh_cat + logRT) + (1 +  precoh_cat*postcoh_cat + logRT|subj), data=bigData_incorrect
                       , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
print(summary(confModel_error))
print(Anova(confModel_error, type =3))
hist(resid(confModel_error))
coef(summary(confModel_error))

## Country interaction on correct trials (Suppfig 2 C/D)
confModel_correct = lmer(conf ~ country * (precoh_cat*postcoh_cat + logRT) + (1 +  precoh_cat*postcoh_cat + logRT|subj), data=bigData_correct
                         , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
print(summary(confModel_correct))
print(Anova(confModel_correct, type =3))

