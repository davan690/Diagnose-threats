####################################
########  SCRIPT: MAKE TABLES TO EXPLORE DIAGNOSIS PERFORMANCE.




####################
########  READ IN DATA

setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\TFP Results")

load("LatestMasterFile.RData")    # read in name of latest master file

  # ls()

master <- read.csv(LatestMasterFile,header=T)

names(master)

#####################
#########  FIRST TABLE...

rownames <- c("null", "weak exploitation", "moderate exploitation", "severe exploitation",
                 "weak habitat loss", "moderate habitat loss", "severe habitat loss")

colnames <- c("Sample size","True Threat Process", "Correct", "Strongly correct", "Wrong", 
              "Highly misleading",   "None", "Exploitation", "Habitat loss")

    #### display only medium variability, 30 years duration, moderate DD, single threat process

yearfilter <- grep("30years",as.character(master$filename))
varfilter <- grep("medvar",as.character(master$filename))
modfilter <- which(master$threat%in%c(1,2,3))
DDfilter <- grep("moderateDD",as.character(master$filename))

                  ### might try running without the constant harvest scenario
harvfilter2 <- grep("Constharv",as.character(master$filename))
harvfilter <- c(1:nrow(master))[-harvfilter2]

filter <- sort(intersect(yearfilter,intersect(varfilter,intersect(modfilter,intersect(DDfilter,harvfilter)))) )

subset1 <- master[filter,]

names(subset1)

nrow(subset1)

rows_severities <- c(0,1,2,3,1,2,3)
rows_threats <- c(1,3,3,3,2,2,2)

newtable <- data.frame(temp=rep(NA,length(rownames)))
row.names(newtable) <- rownames
replicates=numeric(length(rownames))
for(r in 1:length(rownames)){
  ndx <- which((subset1$severity==rows_severities[r])&
               (subset1$threat==rows_threats[r]))
  denom <- (sum(subset1$veryRight[ndx])+sum(subset1$weaklyRight[ndx])+sum(subset1$veryWrong[ndx])+sum(subset1$weaklyWrong[ndx]))
  newtable$Correct[r] <- (sum(subset1$veryRight[ndx])+sum(subset1$weaklyRight[ndx]))/denom
                             
  newtable$sc[r] <- sum(subset1$veryRight[ndx]) / denom
  newtable$Wrong[r] <- (sum(subset1$veryWrong[ndx])+sum(subset1$weaklyWrong[ndx])) / denom
  newtable$hm[r] <- sum(subset1$veryWrong[ndx]) / denom

  newtable$aw_null[r] <- mean(subset1$pw_null[ndx])
  newtable$aw_harv[r] <- mean(subset1$pw_harvest[ndx])
  newtable$aw_hl[r] <- mean(subset1$pw_habloss[ndx])

  replicates[r] <- length(ndx)
}

newtable <- newtable[,-1]   # remove temporary column
 

names(newtable) <- colnames[-1]


newtable


######### Write to file

#filename=paste("PerformanceTable_moderatePower_",Sys.Date(),".csv",sep="")
filename=paste("PerformanceTable_highPower_",Sys.Date(),".csv",sep="")

setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\TFP Results")
write.table(newtable,sep=",",row.names=F,file=filename)



#####################
#########  FIRST TABLE...   REPEAT FOR CONSTANT HARVEST....

rownames <- c("null", "weak exploitation", "moderate exploitation", "severe exploitation",
                 "weak habitat loss", "moderate habitat loss", "severe habitat loss")

colnames <- c("Sample size","True Threat Process", "Correct", "Strongly correct", "Wrong", 
              "Highly misleading",   "None", "Exploitation", "Habitat loss")

    #### display only medium variability, 30 years duration, moderate DD, single threat process

yearfilter <- grep("45years",as.character(master$filename))
varfilter <- grep("lowvar",as.character(master$filename))
modfilter <- which(master$threat%in%c(1,2,3))
DDfilter <- grep("strongDD",as.character(master$filename))

                  ### might try running without the constant harvest scenario
harvfilter2 <- grep("Fracharv",as.character(master$filename))
harvfilter <- c(1:nrow(master))[-harvfilter2]

filter <- sort(intersect(yearfilter,intersect(varfilter,intersect(modfilter,intersect(DDfilter,harvfilter)))) )

subset1 <- master[filter,]

names(subset1)

nrow(subset1)

rows_severities <- c(0,1,2,3,1,2,3)
rows_threats <- c(1,3,3,3,2,2,2)

newtable <- data.frame(temp=rep(NA,length(rownames)))
row.names(newtable) <- rownames
replicates=numeric(length(rownames))
for(r in 1:length(rownames)){
  ndx <- which((subset1$severity==rows_severities[r])&
               (subset1$threat==rows_threats[r]))
  denom <- (sum(subset1$veryRight[ndx])+sum(subset1$weaklyRight[ndx])+sum(subset1$veryWrong[ndx])+sum(subset1$weaklyWrong[ndx]))
  newtable$Correct[r] <- (sum(subset1$veryRight[ndx])+sum(subset1$weaklyRight[ndx]))/denom
                             
  newtable$sc[r] <- sum(subset1$veryRight[ndx]) / denom
  newtable$Wrong[r] <- (sum(subset1$veryWrong[ndx])+sum(subset1$weaklyWrong[ndx])) / denom
  newtable$hm[r] <- sum(subset1$veryWrong[ndx]) / denom

  newtable$aw_null[r] <- mean(subset1$pw_null[ndx])
  newtable$aw_harv[r] <- mean(subset1$pw_harvest[ndx])
  newtable$aw_hl[r] <- mean(subset1$pw_habloss[ndx])

  replicates[r] <- length(ndx)
}

newtable <- newtable[,-1]   # remove temporary column
 

names(newtable) <- colnames[-1]


newtable


######### Write to file

#filename=paste("PerformanceTable_moderatePower_ConstHarv_",Sys.Date(),".csv",sep="")
filename=paste("PerformanceTable_highPower_ConstHarv_",Sys.Date(),".csv",sep="")

setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\TFP Results")
write.table(newtable,sep=",",row.names=F,file=filename)





#######################################
#####################
#########  SECOND TABLE...  explore multiple threats...

names(master)

  #comb_ndx <- grep("Combined",master$filename)
  #length(comb_ndx)  # all ran successfully

    ## to use for table: medium variance only???


rownames <- c("weak exploitation, weak habitat loss", "severe exploitation, weak habitat loss", 
                 "weak exploitation, severe habitat loss",  "severe exploitation, severe habitat loss")

colnames <- c("True Threat Process", "# scenarios","Exploitation (proportion)", "Habitat Loss (proportion)", "Exploitation (weight)", 
              "Habitat Loss (weight)")

    #### display only medium variability, 30 years duration, moderate DD, single threat process

yearfilter <- grep("45years",as.character(master$filename))
varfilter <- c(grep("lowvar",as.character(master$filename)),grep("medvar",as.character(master$filename)),grep("highvar",as.character(master$filename)))
modfilter <- which(master$threat%in%c(4))
DDfilter <- grep("moderateDD",as.character(master$filename))

                  ### might try running without the constant harvest scenario
harvfilter2 <- grep("Constharv",as.character(master$filename))
harvfilter <- c(1:nrow(master))[-harvfilter2]

filter <- sort(intersect(yearfilter,intersect(varfilter,intersect(modfilter,intersect(DDfilter,harvfilter)))) )

subset1 <- master[filter,]

names(subset1)

        # new columns: harvest severity and habitat loss severity
harvlevels <- unique(master$FracHarv)
hllevels <- unique(master$HabLoss)

subset1$severity_hl <- match(subset1$HabLoss,hllevels)-1
subset1$severity_harv <- match(subset1$FracHarv,harvlevels)-1


        #  for subsetting for each row...
rows_severities_hl <- c(1,1,3,3)
rows_severities_harv <- c(1,3,1,3)

newtable <- data.frame(temp=rep(NA,length(rownames)))
row.names(newtable) <- rownames
replicates=numeric(length(rownames))

   # r=1
for(r in 1:length(rownames)){
  ndx <- which((subset1$severity_hl==rows_severities_hl[r])&
               (subset1$severity_harv==rows_severities_harv[r]))
  replicates[r] <- length(ndx)

  newtable$reps[r] <- replicates[r]
  newtable$harv1[r] <- mean(subset1$selected_harv[ndx])
  newtable$hl1[r] <- mean(subset1$selected_habloss[ndx])
  newtable$harv2[r] <- mean(subset1$pw_harvest[ndx])
  newtable$hl2[r] <- mean(subset1$pw_habloss[ndx])
}

newtable <- newtable[,-1]   # remove temporary column
 

names(newtable) <- colnames[-1]


newtable


######### Write to file

filename=paste("PerformanceTable_CombinedThreats_",Sys.Date(),".csv",sep="")
setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\TFP Results")
write.table(newtable,sep=",",row.names=F,file=filename)




##########  END SCRIPT





