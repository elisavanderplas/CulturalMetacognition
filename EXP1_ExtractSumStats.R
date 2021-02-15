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
      accuracy = accuracy -0.5 #to get -0.5 or 0.5

      logRT = scale(log(RT))
      logConfRT = scale(log(confRT))
      subj = rep(subjects[s], length(accuracy))
      country = rep(d, length(accuracy))

      subData1 = data.frame("country" = country, "subj"=subj, "dir"=dir, "precoh"=precoh, "postcoh"=postcoh,"precoh_cat"=precoh_cat, "postcoh_cat"=postcoh_cat, "conf"=conf, "response"=response, "logRT"=logRT, "accuracy"=accuracy)
  
      COM_model_nonsocial = coef(lm(data = subData1, formula = conf ~ precoh_cat + postcoh_cat + logRT + accuracy + accuracy * postcoh_cat))
    
      PDE_nonsocial_pre = abs(COM_model_nonsocial[2])
      PDE_nonsocial_post = abs(COM_model_nonsocial[3])
      PDE_nonsocial_accpost = abs(COM_model_nonsocial[6])
      
      subj = subjects[s]
      country = d
      subData2 = data.frame("country" = country, "subj" = subj, "PDE_nonsocial_pre" = PDE_nonsocial_pre, "PDE_nonsocial_post" = PDE_nonsocial_post, "PDE_nonsocial_accpost" = PDE_nonsocial_accpost)
      
      bigData = rbind(bigData, subData2)
      j=j+1 # total subject counter
    }
}

# Fit regression models
bigData$subj <- factor(bigData$subj)
bigData$country <- factor(bigData$country, levels = c(1,2), labels = c("PKU", "UCL"))

setwd("~/Dropbox/CulturalMetacognition-master/DATA/EXP1/")
write.csv(bigData,file = paste('regression_betas_IDs_EXP1.csv'))
