################################################
################################################
######     R script for analyzing success of diagnoses for TFP project



################ read in the raw threat diagnosis results for processing

setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\R Workspace")
load("DiagnosisResults_Simulated.RData")

ls()




################################################
################# generate master data frame with specs and results

#######################################
###############  determine if any model specifications did not run properly...

badmodels <- numeric(0)
for(i in 1:nmodels){
  temp <- length(bugs.runs[[i]][[1]])
  if(temp<10) badmodels <- c(badmodels,i)
}

badmodels
  #filenames[badmodels]

goodmodels <- c(1:nmodels)[-badmodels]



names(specs)

master <- specs[goodmodels,]     # start master file (carry over all specs for simulated time series generating models)

 # names(bugs.runs[[1]]$details)
 # names(bugs.runs[[1]]$summary)

master$mostSelected <- "0"
master$freqSelected <- 0
master$sdSelected <- 0
master$veryRight <- 0
master$weaklyRight <- 0
master$veryWrong <- 0
master$weaklyWrong <- 0
master$selected_null <- 0
master$selected_habloss <- 0
master$selected_harv <- 0
master$pw_null <- 0
master$K_null <- 0
master$pw_habloss <- 0
master$hlrate_ricker <- 0
master$hlrate_ceiling <- 0
master$pw_harvest <- 0
master$K_harvest <- 0
master$harvrate_frac <- 0
master$harvrate_const <- 0
master$variability  <- 0
master$growthrate_null <- 0
master$growthrate_hl_ricker <- 0
master$growthrate_hl_ceiling <- 0
master$growthrate_harv_frac <- 0
master$growthrate_harv_const  <- 0
master$growthrate_hl_ceiling2 <- 0

###########  loop through diagnosis runs and extract this information
ls()
names(bugs.runs[[1]]$summary)

temp <- bugs.runs[[1]]$summary

  # i=1
for(i in 1:length(goodmodels)){
  temp <- bugs.runs[[goodmodels[i]]]$summary
  details <- bugs.runs[[goodmodels[i]]]$details

  counter = array(0,dim=c(3))
     # r=1
  for(r in 1:nreps){
     mmm <- which.max(details$modselect_array[r,])
     counter[mmm] <- counter[mmm] + 1
  }
  counter = counter/nreps

  master$selected_null[i] <- counter[1]
  master$selected_habloss[i] <- counter[2]
  master$selected_harv[i] <- counter[3]

  master$mostSelected[i] <- temp$most_selected     # most_selected
  master$freqSelected[i] <- temp$freq_selected     # round(freq_selected*nreps,0)
  master$sdSelected[i] <- temp$sd_selected           #  round(wm2_fracharv2,3)
  master$veryRight[i] <- temp$correct_selected[1]    
  master$weaklyRight[i] <- temp$correct_selected[2]
  master$veryWrong[i] <- temp$incorrect_selected[1]
  master$weaklyWrong[i] <- temp$incorrect_selected[2]
  master$pw_null[i] <- temp$constK
  master$K_null[i] <- temp$Kflat1
  master$pw_habloss[i] <- temp$habloss
  master$hlrate_ricker[i] <- temp$Kslope1
  master$hlrate_ceiling[i] <- temp$Kslope2
  master$pw_harvest[i] <- temp$fracharv
  master$K_harvest[i] <- temp$Kflat2
  master$harvrate_frac[i] <- temp$harvrate2
  master$harvrate_const[i] <- temp$harvrate3
  master$variability[i]  <- temp$processerr
  master$growthrate_null[i] <- temp$growthrate1
  master$growthrate_hl_ricker[i] <- temp$growthrate2
  master$growthrate_hl_ceiling[i] <- temp$growthrate4
  master$growthrate_harv_frac[i] <- temp$growthrate3
  master$growthrate_harv_const[i]  <- temp$growthrate5
  master$growthrate_hl_ceiling2[i] <- temp$growthrate6
}
 

################# Write CSV file..


#colnames(temp10) <- c("filename","most selected model","freq selected","strongly select true model",
#                "weakly select true model","strongly select wrong model","weakly select wrong model",
#                "posterior weight, null model", "K, null model", "posterior weight, habitat loss",
#                "hl rate, ricker", "hl rate, ceiling", "posterior weight, exploitation", "K, exploitation",
#                 "harv rate, const frac", "harv rate, const absolute", "env. variability")  # "harv rate, const absolute", 

setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\TFP Results")
filename <- paste("master_",Sys.Date(),".csv",sep="")     # "_",Sys.time(),
write.table(master,file=filename,row.names=F,sep=",")



LatestMasterFile = filename
save(LatestMasterFile,file="LatestMasterFile.RData")

















