## EVDP 2019 elisa.plas.18@ucl.ac.uk

rm(list=ls())
require(R.matlab) 
require(lme4)
require(car)
require(optimx)
options(contrasts = c("contr.treatment", "contr.poly")) # This is R defaults but set it anyway to be safe

bigData = NULL
demoData = NULL
nat = c('PKU', 'UCL')

j=1
for (d in 1:2) {
  if (d == 1) {
    dataDir = "~/Dropbox/CulturalMetacognition_2020/DATA/EXP1/PKU_data/PKU_data/"
    filePrefix = "fMRI_pilotData_sub_"
    suffix = "_2"
    subjects =c(seq(101,109), seq(111,115), seq(117,141))
  } else if (d == 2) {
    dataDir = "~/Dropbox/CulturalMetacognition_2020/DATA/EXP1/UCL_data/UCL_data/"
    filePrefix = "fMRI_pilotData_sub_"
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
# Factors
bigData$subj <- factor(bigData$subj)
bigData$country <- factor(bigData$country, levels = c(1,2), labels = c("PKU","UCL"))

## distinguish between correct/incorrect & PKU/UCL trials
bigData_correct <- bigData[bigData$accuracy == 1, ]
bigData_incorrect <- bigData[bigData$accuracy == -1, ] 
bigData_PKU <- bigData[bigData$country == "PKU",]
bigData_UCL <- bigData[bigData$country == "UCL",]
bigData_PKU_correct <- bigData_correct[bigData_correct$country == "PKU",]
bigData_UCL_correct <- bigData_correct[bigData_correct$country == "UCL",]
bigData_PKU_incorrect <- bigData_incorrect[bigData_incorrect$country == "PKU",]
bigData_UCL_incorrect <- bigData_incorrect[bigData_incorrect$country == "UCL",]

## DEMOGRAPHICS
PKU = c(demoData[1:length(subjects),])
UCL = c(demoData[length(subjects)+1:nrow(demoData),])

mean(PKU$age, na.rm = T)
mean(UCL$age, na.rm = T)
t.test(PKU$age, UCL$age, var.equal = T)
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

# questionnaires
mean(PKU$constr2, na.rm = T)
mean(UCL$constr2, na.rm = T)
t.test(PKU$constr2, UCL$constr2, var.equal = T)
mean(PKU$hol_tot, na.rm = T)
mean(UCL$hol_tot, na.rm = T)
t.test(PKU$hol_tot, UCL$hol_tot,var.equal = T)
mean(PKU$IQ, na.rm = T)
mean(UCL$IQ, na.rm = T)
t.test(PKU$IQ, UCL$IQ, var.equal = T)

## FIRST-ORDER PERFORMANCE
accModel = glmer(accuracy_dv ~ country*precoh_cat + (1 + precoh_cat | subj), data=bigData, family="binomial")
print(summary(accModel))
print(Anova(accModel, type =3))

## SECOND-ORDER PERFORMANCE
confModel_nocountry = lmer(conf ~ accuracy + precoh_cat + postcoh_cat + precoh_cat:postcoh_cat + precoh_cat:accuracy + postcoh_cat:accuracy + precoh_cat:postcoh_cat:accuracy + logRT + (1 + accuracy + precoh_cat + postcoh_cat + precoh_cat:postcoh_cat + precoh_cat:accuracy + postcoh_cat:accuracy + precoh_cat:postcoh_cat:accuracy + logRT|subj), data=bigData
                        , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))

confModel_wcountry = lmer(conf ~ country*(accuracy + precoh_cat + postcoh_cat + precoh_cat:postcoh_cat + precoh_cat:accuracy + postcoh_cat:accuracy + precoh_cat:postcoh_cat:accuracy + logRT) + (1 + accuracy + precoh_cat + postcoh_cat + precoh_cat:postcoh_cat + precoh_cat:accuracy + postcoh_cat:accuracy + precoh_cat:postcoh_cat:accuracy + logRT|subj), data=bigData
                 , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE, REML = FALSE)))

fix <- fixef(confModel_wcountry)
print(summary(confModel_wcountry))
print(Anova(confModel_wcountry, type = 3))
coef(summary(confModel_wcountry)) #get the contrast statistics
fix.se <- sqrt(diag(vcov(confModel_wcountry))) ## 2-way interaction, no 3 way interaction, line 463 & 471

## check if including country as interaction improved the fit of the model (line 457)
anova(confModel_nocountry,confModel_wcountry)


## Get the beta coefficients per country (Fig 1D plotted)
# error, PKU
error_coef_PKU = lmer(conf ~ precoh_cat*postcoh_cat + logRT + (1 +  precoh_cat*postcoh_cat + logRT|subj), data=bigData_PKU_incorrect
                       , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(error_coef_PKU)
fix.se <- sqrt(diag(vcov(error_coef_PKU)))
betas <- c(fix, fix.se)
setwd("~/Dropbox/CulturalMetacognition_2020/DATA/EXP1/PKU_data/PKU_data/PKU_betas/")
write.csv(betas, file = paste('regression_betas_err_PKU.csv'))

# correct, PKU 
corr_coef_PKU = lmer(conf ~ precoh_cat*postcoh_cat + logRT + (1 +  precoh_cat*postcoh_cat + logRT|subj), data=bigData_PKU_correct
                      , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(corr_coef_PKU)
fix.se <- sqrt(diag(vcov(corr_coef_PKU)))
betas <- c(fix, fix.se)
write.csv(betas, file = paste('regression_betas_corr_PKU.csv'))

# error, UCL
error_coef_UCL = lmer(conf ~ precoh_cat*postcoh_cat + logRT + (1 +  precoh_cat*postcoh_cat + logRT|subj), data=bigData_UCL_incorrect
                      , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(error_coef_UCL)
fix.se <- sqrt(diag(vcov(error_coef_UCL)))
betas <- c(fix, fix.se)
setwd("~/Dropbox/CulturalMetacognition_2020/DATA/EXP1/UCL_data/UCL_data/UCL_betas/")
write.csv(betas, file = paste('regression_betas_err_UCL.csv'))

# correct, UCL
corr_coef_UCL = lmer(conf ~ precoh_cat*postcoh_cat + logRT + (1 +  precoh_cat*postcoh_cat + logRT|subj), data=bigData_UCL_correct
                     , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
fix <- fixef(corr_coef_UCL)
fix.se <- sqrt(diag(vcov(corr_coef_UCL)))
betas <- c(fix, fix.se)
write.csv(betas, file = paste('regression_betas_corr_UCL.csv'))

## Country interaction on error trials (Fig 1D star)
confModel_error = lmer(conf ~ country * (precoh_cat*postcoh_cat + logRT) + (1 +  precoh_cat*postcoh_cat + logRT|subj), data=bigData_incorrect
                       , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
print(summary(confModel_error))
print(Anova(confModel_error, type =3))

## Country interaction on correct trials (Fig 1D 'n.s.')
confModel_correct = lmer(conf ~ country * (precoh_cat*postcoh_cat + logRT) + (1 +  precoh_cat*postcoh_cat + logRT|subj), data=bigData_correct
                       , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
print(summary(confModel_correct))
print(Anova(confModel_correct, type =3))



