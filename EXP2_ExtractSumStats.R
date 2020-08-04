## EVDP 2019 elisa.plas.18@ucl.ac.uk

rm(list=ls())
require(R.matlab) 
require(lme4)
require(car)
require(optimx)

nat = c('PKU', 'UCL')
bigData = NULL

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
       # replace nonsocial data w NA
      conf_adv_dir[conf_adv_dir == 99] <- NA
      
      # get each inidiv's binned conf
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
      accuracy = accuracy -0.5 #to get -0.5 or 0.5
      
      a_adv = dat[,,1]$a.adv[1,]
      # replace nonsocial data w NA
      a_adv[a_adv == 99] <- NA
      a_adv = a_adv-1# get from 1/2 to 0/1
      agree = a_adv == response
      agree = agree -0.5 #to get -0.5 or 0.5

      logRT = scale(log(RT))
      logConfRT = scale(log(confRT))
      subj = rep(subjects[s], length(accuracy))
      task = dat[,,1]$condition[1,]
      country = rep(d, length(accuracy))
      
      subData1 = data.frame("country" = country, "subj"=subj, "task"=task, "dir"=dir, "precoh"=precoh, "postcoh"=postcoh,"precoh_cat"=precoh_cat, "postcoh_cat"=postcoh_cat, "conf"=conf, "conf_adv" = conf_adv, "conf_adv_dir" = conf_adv_dir, "conf_adv_cat" = conf_adv_cat, "conf_dir" = conf_dir, "logConfRT"=logConfRT, "response"=response, "logRT"=logRT, "accuracy"=accuracy, "a_adv" = a_adv, "agree" = agree)
  
      #get the summary stat social and non-social
      subData_social <- subData1[subData1$task == 1,]
      subData_nonsocial <- subData1[subData1$task == 0,]
      
      COM_model_social = coef(lm(data = subData_social, formula = conf ~ precoh_cat + conf_adv_cat + logRT + accuracy + accuracy * conf_adv_cat))
      COM_model_nonsocial = coef(lm(data = subData_nonsocial, formula = conf ~ precoh_cat + postcoh_cat + logRT + accuracy + accuracy * postcoh_cat))
      
      PDE_social_pre = abs(COM_model_social[2])
      PDE_social_post = abs(COM_model_social[3])
      PDE_social_accpost = abs(COM_model_social[6])
      PDE_nonsocial_pre = abs(COM_model_nonsocial[2])
      PDE_nonsocial_post = abs(COM_model_nonsocial[3])
      PDE_nonsocial_accpost = abs(COM_model_nonsocial[6])
      
      subj = subjects[s]
      country = d
      subData2 = data.frame("country" = country, "subj" = subj, "PDE_social_pre" = PDE_social_pre, "PDE_social_post" = PDE_social_post, "PDE_nonsocial_pre" = PDE_nonsocial_pre, "PDE_nonsocial_post" = PDE_nonsocial_post, "PDE_social_accpost" = PDE_social_accpost, "PDE_nonsocial_accpost" = PDE_nonsocial_accpost)
      
      bigData = rbind(bigData, subData2)
      j=j+1 # total subject counter
    }
}

# Fit regression models
bigData$subj <- factor(bigData$subj)
bigData$country <- factor(bigData$country, levels = c(1,2), labels = c("PKU", "UCL"))

setwd("~/Dropbox/CulturalMetacognition_2020/DATA/EXP2/")
write.csv(bigData,file = paste('regression_betas_IDs_EXP2.csv'))
