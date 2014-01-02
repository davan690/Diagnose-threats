############  Goal: to infer threat from pattern. 


####################################
####################################
#####  DIAGNOSE THREATS USING BAYESIAN ALGORITHM



####################################
####################################

library(R2WinBUGS)
library(snow)
library(snowfall)
  #?snowfall
library(foreach)
library(doSNOW)


#################
####    SET WORKING DIRECTORY

#setwd("C:/Users/Kevin/Documents/Academic/Employment/Stony Brook/Proposals/Threat from Pattern")
#setwd("C:/Users/Kevin/Desktop/Kevin/TFP_transfer")
#setwd("C:/Users/grads/Dropbox/TFP stuff")
#setwd("C:/Users/Kevin/Dropbox/TFP stuff")
#setwd("C:/Documents and Settings/JS/My Documents/Dropbox/TFP stuff")
setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\Simulated_TS")

############################3
#################################################################################


#######################################################
######################## READ IN SIMULATED DATA


#setwd("C:/Users/Kevin/Documents/Academic/Employment/Stony Brook/Proposals/Threat from Pattern/Time series4_PVA")
#setwd("C:/Users/grads/Desktop/Timeseries-set3")

 #setwd("C:\\Users\\Kevin\\Desktop\\TFP_DATA\\SimulatedTS")
setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\Simulated_TS")


temp <- read.csv("mp_filedata.txt",header=T,sep="")

##"filename"    "nyears"      "threat"      "severity"    "variability"     "DDStrength"  "DDType" 

nfiles <- length(temp[,1])  #18 #24
starttimes <- rep(1,times=nfiles)
lengths    <- temp$nyears
filenames  <- temp$filename
endtimes <- lengths
variabilities <- temp$variability
DDStrengths <- temp$DDStrength
DDTypes <- temp$DDType
 #severities <- temp$severity
truemodel <- temp$threat
threatname <- as.character(temp$threatname)

model <- as.character(filenames)


nmodels = nfiles #18

nreps = 50

ny = temp$nyears

gentime = c(rep(1,times=nfiles))

#  read all time series data into master array

ts_list <- list()

nyears <- array(0, dim=c(nmodels,nreps))
startabun <- array(0, dim=c(nmodels,nreps))
for(f in 1:nfiles){
  filename = paste("model",f,".csv",sep="")
  temp <- read.csv(filename,header=T)
  ts_list[[f]] <- list()
  for(r in 1:nreps){
    temp2 <- as.numeric(temp[r,1:ny[f]])[which(as.numeric(temp[r,1:ny[f]])>0)]
    nyears1 <- length(temp2)
    nyears[f,r] <- floor(length(temp2)/gentime[f])  # modify number of years to account for generation time
    startabun[f,r] <- mean(temp2[1:3])
    ts_list[[f]][[r]] <- array(0,dim=c(gentime[f],nyears[f,r]))
    for(g in 1:gentime[f]){
      sss <- seq(g,nyears1,by=gentime[f])[1:nyears[f,r]]
      ts_list[[f]][[r]][g,] <- temp2[sss]
    }
  }    
}


setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\BUGS files")

############################################################
###########################################################################################

  #### weed out inappropriate models that muddy the waters too much...
     # harvest with moderate or strong DD...     models 10 to 18
     # habitat loss with weak DD (doesn't pick up on the habitat loss signal)  models 38 to 46
     
      #modelset <- c(1:nmodels)[-c(10:18,38:46)]



##########################################################
############               INITIALIZE GLOBAL PARAMETERS


#  run each time series through BUGS and write results to a summary file

nmodels <- nfiles         #18   #24  # length(modelset)  # 
                          # define an object that is stored for each model
nreps <- 50  # 10 # 
ncandidates <- 3          #5


   #  current_m <- m
   # current_r <- r


truethreat = c("null (no threat)", "habitat Loss", "exploitation", "combined")

cand_models <- c("null (no threat)","habitat loss","exploitation")


#####################################################
###############   MAKE FUNCTION FOR PARALLELIZING PROCESS

 ### m=285;r=13

runSimModel <- function(m,nreps,ncandidates,truethreat,cand_models,truemodel,nfiles,
                     starttimes,lengths,filenames,endtimes,ts_list){

  
      #  instantiate data storage structure
  AllReps <- list()
  AllReps$details <- list()
  AllReps$summary <- list()    # summary for final table...

  AllReps$details$selected  <- array(0,dim=c(nreps))  
  AllReps$details$modselect_array <- array(0,dim=c(nreps,ncandidates))
  AllReps$details$beta_array <- array(0,dim=c(nreps,6))
  AllReps$details$slope_array <- array(0,dim=c(nreps,2))        
  AllReps$details$Kflat_array <- array(0,dim=c(nreps,2))
  AllReps$details$sd_array <- array(0,dim=c(nreps,ncandidates))
   #startyear_array <- array(0,dim=c(nreps,ncandidates))
   #length_array <- array(0,dim=c(nreps,ncandidates))
  AllReps$details$constharv_array <- array(0,dim=c(nreps))
  AllReps$details$fracharv_array <- array(0,dim=c(nreps))
    #fracharv2_array <- array(0,dim=c(nreps))

  AllReps$details$hltype_array <- array(0,dim=c(nreps))
  AllReps$details$harvtype_array <- array(0,dim=c(nreps))



  for(r in 1:nreps){   


    ##############
    #############

    ### m=288;r=14
 

    gentime = 1    #  approximate generation time, 
    realspecname <- truethreat[truemodel[m]] # "Simulated time series"
    #reallocation <- ifelse(truemodel[m]==3,paste("(",substr(filenames[m],3,3),"% per year)",sep=""),
         #                ifelse(truemodel[m]==2,paste("(K reduced by ",substr(filenames[m],3,3)," per year)",sep=""),""))  #""     # \n starting at year 15

    reallocation <- threatname[m]
    realsource <- filenames[m]    
    startyear = starttimes[m]
    endyear = nyears[m,r]   #lengths[m]    # 15  #  floor((nrow(realdata))/gentime)-1  # 
    yaxname = "Abundance index"
    realyears = seq(1,nyears[m,r],1)   #  realdata[,1]

    ts <- ts_list[[m]][[r]] *(100/max(ts_list[[m]][[r]]))  # rescale
    ncandidates <- 3
    StartTimes <- array(0,dim=nyears[m,r])
    StartTimes[startyear] <- 1
    Lengths <- array(0,dim=nyears[m,r])
    Lengths[(endyear-startyear)] <- 1
    realhabitat <- matrix(seq(1,nyears[m,r],1),nrow=gentime,ncol=nyears[m,r])
    realtime <- realhabitat
     #if(endyear==15) realhabitat[1,15:30] <- 15
     #if(startyear==15) realhabitat[1,1:30] <- c(rep(1,times=15),seq(2,16,1))

                           # plot out the time series
    PLOT = F # PLOT=T
    if(PLOT){
      par(mai=c(.75,1,.3,.3))
      plot(as.vector(ts),type="l",ylim=c(0,max(ts)*1.1),lwd=2,col=gray(.2),xlab="",
                  xaxt="n",ylab=yaxname)   # first plot out the original time series time (years), abundance
      axis(1,at=c(1:nyears[m,r]),labels=realyears,tick=F)
      text(endyear*gentime-5,95,realsource,cex=.8,col=gray(.3))
      text(nyears[m,r]/2,106,paste("Simulated threat scenario: ",realspecname," ",reallocation,sep=""),cex=1.2)
    }

    year <- seq(1,nyears[m,r],1)
    slope <- summary( lm( as.vector(ts)[startyear:(endyear*gentime)]~c(startyear:(endyear*gentime)) ) )$coefficients[2,1]
    hlhigh <- ifelse(slope<0,slope*1.5,-slope)
    hllow  <- ifelse(slope<0,slope/5,-slope/5)    # was divided by 5
    nullhigh <- (mean(ts[1:2])*0.1)/nyears[m,r]
    nulllow <- -1*(mean(ts[1:2])*0.1)/nyears[m,r]
    if (abs(hllow)<0.6) hllow <- -0.6 
    if (abs(hlhigh)<1) hlhigh <- -1
    stdev <- sd(log(as.vector(ts)))
     #sdlow <- stdev/100    # 0.01   #
     #sdhigh <- stdev       #   0.01  #

    betas <- log(as.vector(ts))[-1] - log(as.vector(ts))[-length(ts)] 
    maxbeta <- max(betas) 
    highbeta2 <- maxbeta*2    # for the habitat loss model with Ricker DD...  // was maxbeta*2
    highbeta <- maxbeta/2    # maxbeta/2 seems to work well?*** but remember that this is really a maximum growth rate, and in habitat loss scenarios is never reached even close
    medbeta <- maxbeta/4   # maxbeta/5 seems to work well?***
    lowbeta <- maxbeta/10
    minbeta <- min(betas[which(betas>0)])
    sdlow  <- sd(betas)/1000    # 1000 seems to work well?***
    sdhigh  <- sd(betas)*50    # 15 seems to work well?***


                     ######  generic BUGS code
    Data <- list( 
             y = ts,
             nyears = nyears[m,r],
             realhabitat = realhabitat, 
             realtime = realtime,
             gentime=gentime,
             startabun = mean(ts[1:2]),
              #endabun = startabun[m,r]-100,
             nmodels = ncandidates,
             startyear2 = as.vector(StartTimes),
             length2 = as.vector(Lengths),
             hlhigh = hlhigh,
             hllow = hllow,
             nullhigh = nullhigh,
             nulllow = nulllow,
	            #hlsdhigh = hlsdhigh,
                  #hlsdlow = hlsdlow,
             hrhigh = -1*hlhigh,
             hrlow  = 0.01, #-1*hllow,
                  #hrsdhigh = hlsdhigh,
                  #hrsdlow = hlsdlow,
                  #harvmax = slope*10,
                  #minbeta = minbeta,
             lowbeta = lowbeta,
             medbeta = medbeta,
             highbeta = highbeta,
             highbeta2 = highbeta2,
             sdlow = sdlow,
             sdhigh = sdhigh
    )

    Inits <- function() list(       #  initial values for all stochastic nodes.  
                mod = round(1,0),
                candbeta = c(3,3,3,3,3,3),
                K0 = c(mean(ts[1:2]),NA,mean(ts[1:2])*(17/15)),
                hlrate2 = hllow+(hlhigh-hllow)/2,
                hlrate3 = hllow+(hlhigh-hllow)/2,
                candDD = 1,
                fracharv2 = max(1-exp(mean(betas)),.01),
                constharv2 = -1*(hllow+(hlhigh-hllow)/2) ,  # abs(slope),
                candharv = 1,
                startyear = rep(startyear,3),
                length = rep(endyear-startyear,3),
                sd = rep((stdev/2),3)
    )

    Par <- c(
         "beta",
         #"hlmean",
         #"harvratemean",
         #"harvratesd",
         #"harvratemean2",
         #"harvratesd2",
         #"ceilinggrowth",
         #"hlsd", 
         "mod",
         "logN.new",
         "logN.new2",
         "K0",
         "sd", 
         "constharv2",
         "fracharv2",
         "hlrate2",
         "hlrate3",
         "nullslope2",
            #"overK",
         #"ceilinggrowth",
         #"rickergrowth",
            #"constharv3",
         "harvtype",
         "nulltype",
         "hltype"
           #"fracharv2",
           #"constharv0",
           #"fracharv0",
           #"startyear",
           #"length",
           #"threat_final"
    ) 

    BugFile <- ("tfp9_realdata8_mixture.bug")
    
    BugDir <- "C:\\Users\\Kevin\\Desktop\\Kevin\\WinBUGS14"  #   "C:\\Users\\Kevin\\Documents\\Employment\\ESF\\Bog Turtle\\DATA\\software\\BUGS\\WinBUGS14"   #   

    Mod <- bugs(data=Data, inits=Inits, parameters.to.save=Par,   #  run BUGS with 3 chains  
                   model.file=BugFile, n.chains=1, n.iter= 5000, #5000,   # 20000, 10000, 10          # use 5000 iterations, 2500 burnin
                   bugs.directory=BugDir,
                   n.burnin=2500,n.thin=2) #,debug=T)    # 2500 ,debug=T

      #    "C:/Users/Kevin/Desktop/Kevin/WinBUGS14",
       #    "C:/Users/grads/Desktop/winbugs14/WinBUGS14", 
         #  "C:/Users/Kevin/WinBUGS/winbugs14/WinBUGS14",
        # "C:/Users/Kevin/Documents/Academic/Bog Turtle/DATA/software/BUGS/WinBUGS14" #


    #############
    #############

    ####################################  write out the important information for each replicate: 
                                             
    AllReps$details$selected[r] <- as.numeric(names(which.max(table(Mod$sims.list$mod))))
    temp  <- table(Mod$sims.list$mod)
    temp2 <- as.numeric(dimnames(temp)[[1]])
    temp3 <- as.numeric(temp)
    for(c in 1:ncandidates){  
        if(c%in%temp2) AllReps$details$modselect_array[r,c] <- temp3[which(temp2==c)]  
    }
    AllReps$details$beta_array[r,1:6] <- c(mean(Mod$sims.list$beta[,1]),mean(Mod$sims.list$beta[,2]),mean(Mod$sims.list$beta[,3]),mean(Mod$sims.list$beta[,4]),mean(Mod$sims.list$beta[,5]),mean(Mod$sims.list$beta[,6]))   # slope_array, K_flat, K_final, sd_array, startyear_array, endyear_array, constharv_array, fracharv_array
    
    AllReps$details$slope_array[r,1:2] <- c(mean(Mod$sims.list$hlrate2),mean(Mod$sims.list$hlrate3))
    AllReps$details$Kflat_array[r,1:2] <- c(mean(Mod$sims.list$K0[,1]),mean(Mod$sims.list$K0[,2]))
      # allee_array[m,r] <- mean(Mod$sims.list$alleestr)
      # K_final[m,r] <- mean(Mod$sims.list$Kt)
    AllReps$details$sd_array[r,1:ncandidates] <- mean(Mod$sims.list$sd[,1:ncandidates])
      # startyear_array[m,r,1:ncandidates] <- as.numeric(names(which.max(table(Mod$sims.list$startyear[,1:ncandidates]))))
      # length_array[m,r,1:ncandidates]   <- as.numeric(names(which.max(table(Mod$sims.list$length[,1:ncandidates]))))
    AllReps$details$constharv_array[r] <- c(mean(Mod$sims.list$constharv2))  #,mean(Mod$sims.list$constharv0))
    AllReps$details$fracharv_array[r] <- c(mean(Mod$sims.list$fracharv2))  #,mean(Mod$sims.list$fracharv0) ) 
       # fracharv2_array[m,r] <- c(mean(Mod$sims.list$fracharv2))  #,mean(Mod$sims.list$fracharv0) )  

       #wm_Kslope[m,r]  <- sum(Mod$sims.list$hlrate[which(Mod$sims.list$mod==2)])/nsims
       #wm_fracharv[m,r]  <- sum(Mod$sims.list$fracharv[which(Mod$sims.list$mod==3)])/nsims

    AllReps$details$hltype_array[r] <- as.numeric(names(which.max(table(Mod$sims.list$hltype))))
    AllReps$details$harvtype_array[r] <- as.numeric(names(which.max(table(Mod$sims.list$harvtype))))

  ##############################################
  }   ## end loop through replicates


  #####################################
  ##############  Summarize the important information over all replicates

  nsims <- nrow(Mod$sims.list$beta)
  bm_threshold <- nsims * 0.5
  truemodel2 <- truemodel[m]


                                         ####  new storage variables
  AllReps$summary$correct_selected    <-    array(0,dim=c(2))            # (strong, weak)
  AllReps$summary$correct_selected    <-    array(0,dim=c(2))             # (strong, weak)
  AllReps$summary$incorrect_selected  <-    array(0,dim=c(2)) 
  AllReps$summary$incorrect_selected  <-    array(0,dim=c(2))
                                    
  tot_samples <- length(Mod$sims.list$mod)*nreps
  constK2 <- numeric(nreps)  
  habloss2 <- numeric(nreps)
  fracharv2 <- numeric(nreps)
  selected <- numeric(nreps)
   
    #r=1
  for(r in 1:nreps){
    constK2[r] = AllReps$details$modselect_array[r,1]
    habloss2[r] = AllReps$details$modselect_array[r,2]
    fracharv2[r] = AllReps$details$modselect_array[r,3]
    selected[r] <- which.max(AllReps$details$modselect_array[r,])
    if(truemodel2%in%c(1:3)){    # IF NOT COMBINED THREAT... 
      if( (sum(AllReps$details$modselect_array[r,truemodel2])>=bm_threshold)&
        (selected[r]%in%truemodel2) ) AllReps$summary$correct_selected[1] <- AllReps$summary$correct_selected[1] + 1 
      if( (sum(AllReps$details$modselect_array[r,truemodel2])<bm_threshold)&
        (selected[r]%in%truemodel2)) AllReps$summary$correct_selected[2] <- AllReps$summary$correct_selected[2] + 1 
      if( (AllReps$details$modselect_array[r,selected[r]]>=bm_threshold)&!(selected[r]%in%truemodel2)) AllReps$summary$incorrect_selected[1] <- AllReps$summary$incorrect_selected[1] + 1 
      if( (AllReps$details$modselect_array[r,selected[r]]<bm_threshold)&!(selected[r]%in%truemodel2)) AllReps$summary$incorrect_selected[2] <- AllReps$summary$incorrect_selected[2] + 1 
    }
  }
  AllReps$summary$model <- as.character(filenames[m])
  AllReps$summary$most_selected <- cand_models[as.numeric(names(which.max(table(selected))))]
  AllReps$summary$freq_selected <- as.numeric(table(selected)[which.max(table(selected))])/nreps

  sdtemp <- cbind(constK2/length(Mod$sims.list$mod), habloss2/length(Mod$sims.list$mod), fracharv2/length(Mod$sims.list$mod))
  AllReps$summary$sd_selected <- sd(sdtemp[,as.numeric(names(which.max(table(selected))))])
  AllReps$summary$constK <- sum(constK2)/tot_samples
  AllReps$summary$habloss <- sum(habloss2)/tot_samples
  AllReps$summary$fracharv <- sum(fracharv2)/tot_samples
   #fracharv_2[m] <- sum(fracharv3)/tot_samples
  AllReps$summary$Kslope1 <- mean(AllReps$details$slope_array[,1])
  AllReps$summary$Kslope2 <- mean(AllReps$details$slope_array[,2])
  AllReps$summary$harvrate2 <- mean(AllReps$details$fracharv_array[])
  AllReps$summary$harvrate3 <- mean(AllReps$details$constharv_array[])
  AllReps$summary$growthrate1 <- mean(AllReps$details$beta_array[,1])
  AllReps$summary$growthrate2 <- mean(AllReps$details$beta_array[,2])
  AllReps$summary$growthrate3 <- mean(AllReps$details$beta_array[,3])
  AllReps$summary$growthrate4 <- mean(AllReps$details$beta_array[,4])
  AllReps$summary$growthrate5 <- mean(AllReps$details$beta_array[,5])
  AllReps$summary$growthrate6 <- mean(AllReps$details$beta_array[,6])

  AllReps$summary$Kflat1    <- mean(AllReps$details$Kflat_array[,1])
  AllReps$summary$Kflat2    <- mean(AllReps$details$Kflat_array[,2])
  AllReps$summary$processerr  <- mean(AllReps$details$sd_array[,as.numeric(names(which.max(table(selected))))])

  filename <- paste(filenames[m],"_model",m,"_",Sys.Date(),".RData",sep="")
  #filename <- "temp"

  setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\DUMP")
  save(Mod,AllReps,file=filename)
  setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\BUGS files")

  return(AllReps)
}   



######################################
##############  LOOP THROUGH ALL SCENARIOS (PARALELLIZED)    use packages foreach and snow

##############   NOTE: this can take a long time!


cl <- makeCluster(4,"SOCK")
registerDoSNOW(cl)

bugs.runs <- foreach(x = 1:nmodels,     # nmodels
               .packages = c("R2WinBUGS"),
               .errorhandling=c("pass")
                                         ) %dopar% {   #    %dopar% {   #    nmodels
  runSimModel(x,nreps,ncandidates,truethreat,cand_models,truemodel,nfiles,   #nreps,nmodels
                     starttimes,lengths,filenames,endtimes,ts_list)
}

stopCluster(cl)

  #stack.out <- stack(bugs.runs)
  #names(stack.out)


######################################################
#############

 #nyears2 = nyears[1,1]
gentime2 = gentime[1]

  #save(list = ls(all=TRUE), file = "BugsOutput.RData")



#######################################
###############  determine if any model specifications did not run properly...

PREV_RUN <- T

if(PREV_RUN){
  setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\R Workspace")
  load("DiagnosisResults_Simulated.RData")
  setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\BUGS files")
}


badmodels <- numeric(0)
for(i in 1:nmodels){
  temp <- length(bugs.runs[[i]][[1]])
  if(temp<10) badmodels <- c(badmodels,i)
}

badmodels
filenames[badmodels]

length(badmodels)

##################################
################# re-run any bad models

cl <- makeCluster(4,"SOCK")
registerDoSNOW(cl)

bugs.runs2 <- foreach(x = 1:length(badmodels),          # nmodels
               .packages = c("R2WinBUGS"),
               .errorhandling=c("pass")
                                         ) %dopar% {   #    %dopar% {   #    nmodels
  runSimModel(m=badmodels[x],nreps,ncandidates,truethreat,cand_models,truemodel,nfiles,   #nreps,nmodels
                     starttimes,lengths,filenames,endtimes,ts_list)
}

stopCluster(cl)


#################################
###############  add bad models to the main bugs.runs object..

for(b in 1:length(badmodels)){
  bugs.runs[[badmodels[b]]] <- bugs.runs2[[b]]
}


##################################
################# run the mixture models only


mixmodels <- which(threatname=="combined") 


cl <- makeCluster(4,"SOCK")
registerDoSNOW(cl)

bugs.runs3 <- foreach(x = 1:length(mixmodels),          # nmodels
               .packages = c("R2WinBUGS"),
               .errorhandling=c("pass")
                                         ) %dopar% {   #    %dopar% {   #    nmodels
  runSimModel(m=mixmodels[x],nreps,ncandidates,truethreat,cand_models,truemodel,nfiles,   #nreps,nmodels
                     starttimes,lengths,filenames,endtimes,ts_list)
}

stopCluster(cl)


bugs.

#######################################################
###############    SAVE KEY INFORMATION FOR FURTHER ANALYSIS

        #### make sure all models specs are carried through

setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\Simulated_TS")
specs <- read.csv("mp_filedata.txt",header=T,sep="")


         ## not necessary- just read in from the latest MakeMP file
#specs <- data.frame(temp=numeric(length(filenames)))
#specs$filename <- filenames
#specs$DDStrength <- DDStrengths
#specs$tsLength <- endtimes
#specs$truethreat <- threatname
#specs$truemodelID <- truemodel
#specs$variabilities <- variabilities


   ## add a "Severity" column

strongndx <- c(grep("HL200",as.character(specs$filename)),grep("Constharv400",as.character(specs$filename)),grep("Fracharv6",as.character(specs$filename)))
moderatendx <- c(grep("HL125",as.character(specs$filename)),grep("Constharv275",as.character(specs$filename)),grep("Fracharv4.5",as.character(specs$filename)))
weakndx <- c(grep("HL50",as.character(specs$filename)),grep("Constharv150",as.character(specs$filename)),grep("Fracharv3",as.character(specs$filename)))
specs$severity <- 0
specs$severity[strongndx] <- 3
specs$severity[moderatendx] <- 2
specs$severity[weakndx] <- 1


true_models <- truethreat

setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\R Workspace")
save(cand_models,
     true_models,
     gentime2,
     nmodels,
     nreps,
     ts_list,  # all simulated time series
     specs,
     bugs.runs,
     file="DiagnosisResults_Simulated.RData"
)

##### NOTE: go to Rscript_processSimResults.r for further data processing...



##### 
setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\R Workspace")

save(ts_list,file="SimulatedTimeSeries.RData")





##################################################################################
##################################################################################
#   TEST WITH REAL TIME SERIES

# read in data


realhabitat_df <- read.table("realdata_skylark_habloss1.txt",sep="",header=F)


gentime2 = 1    #  approximate generation time, 
realspecname <- "European Skylark"  #  "Simulated dataset, exploitation, 6% per year"  #  ""Southern bluefin tuna"  #  "Red grouse"  #"Atlantic cod"  #  "Thorny skate"  #  
reallocation <- "in Great Britain" # "" # "" #  "in Langholm moor, Scotland"  # "in Greenland"  #  "in North Atlantic"  #   
realsource <- "BTO/JNCC/RSPB Breeding Bird Survey" #  "" #  "GPDD, Imperial College"  #  "(Thirgood et al. 2000)"  #  "NEFSC (www.nefsc.noaa.gov)"  #  
realfilename <- "realdata_skylark3.txt"  #  "realdata_SBT1.txt"  # "realdata_grouse1.txt" #  "realdata_cod2.txt"  # "realdata_skate1.txt" #  "realdata_cod2.txt"  #  "realdata_farmbirds1.txt"  #
realdata <- read.table(realfilename,header=F,sep="") 

 #realdata <- realdata[-which(realdata[,1]%in%c(1940:1948)),]     ### for grouse data, remove years 1940 to 1947
 
startyear = 1
endyear = floor((nrow(realdata))/gentime2)-1  # 15  #   
yaxname = "Abundance index"
realyears = realdata[,1]

realhabitat <- array(0,dim=c(gentime2,nrow(realdata)))
realhabitat[1,] <- realhabitat_df[,2]

ts <- array(0,dim=c(gentime2,floor(nrow(realdata)/gentime2)))  # c(1,30)
scaling_factor <- 100/max(realdata[,2])
for(i in 1:gentime2){
  ts[i,] <- realdata[seq(i,nrow(realdata),gentime2),2] * scaling_factor      # re-scale to approximately 1000
}

#ts <- ts_list[[34]][[11]] *(100/max(ts_list[[1]][[1]]))    # if real dataset
#endyear <- 30

nyears2 <- ncol(ts)
ncandidates <- 3
StartTimes <- array(0,dim=nyears2)
StartTimes[startyear] <- 1
Lengths <- array(0,dim=nyears2)
Lengths[(endyear-startyear)] = 1



# run BUGS code.
                           # plot out the real time series
#par(mai=c(.75,1,.3,.3))
#plot(as.vector(ts),type="l",ylim=c(0,max(ts)*1.1),lwd=2,col=gray(.2),xlab="",
#              xaxt="n",ylab=yaxname)   # first plot out the original time series time (years), abundance
#axis(1,at=c(1:nrow(realdata)),labels=realyears,tick=F)
#text(endyear*gentime2-5,95,realsource,cex=.8,col=gray(.3))
#text(nrow(realdata)/2,106,paste(realspecname,reallocation),cex=1.2)


year <- seq(1,nyears2,1)
#slope <- summary( lm( as.vector(ts)[startyear:(endyear*gentime2)]~c(startyear:(endyear*gentime2)) ) )$coefficients[2,1]
slope <- summary( lm( as.vector(ts)[startyear:(endyear*gentime2)]~as.vector(realhabitat)[startyear:(endyear*gentime2)] ) )$coefficients[2,1]
hlhigh <- slope*1.5
hllow  <- slope/5

stdev <- sd(log(as.vector(ts)))
 #sdlow <- stdev/100    # 0.01   #
 #sdhigh <- stdev       #   0.01  #

betas <- log(as.vector(ts))[-1] - log(as.vector(ts))[-length(ts)] 
maxbeta <- max(betas) 
highbeta <- maxbeta/2    # maxbeta/2 seems to work well?***
medbeta <- maxbeta/5   # maxbeta/5 seems to work well?***
lowbeta <- maxbeta/10
minbeta <- min(betas[which(betas>0)])
sdlow  <- sd(betas)/1000    # 1000 seems to work well?***
sdhigh  <- sd(betas)*50    # 15 seems to work well?***

#highbeta = 0.18   # debugging... 
#medbeta = 0.06
#lowbeta = 0.09
#sdlow = 0.02
#sdhigh = 0.2


                     ######  generic BUGS code
Data <- list( 
             y = ts,
             nyears = nyears2, 
             gentime=gentime2,
             startabun = mean(ts[1:2]),
              #endabun = startabun[m,r]-100,
             nmodels = ncandidates,
             startyear2 = as.vector(StartTimes),
             length2 = as.vector(Lengths),
             hlhigh = hlhigh,
             hllow = hllow,
             realhabitat = realhabitat,    #***
	            #hlsdhigh = hlsdhigh,
                  #hlsdlow = hlsdlow,
             hrhigh = -1*hlhigh,
             hrlow  = -1*hllow,
                  #hrsdhigh = hlsdhigh,
                  #hrsdlow = hlsdlow,
                  #harvmax = slope*10,
                  #minbeta = minbeta,
             lowbeta = lowbeta,
             medbeta = medbeta,
             highbeta = highbeta,
             sdlow = sdlow,
             sdhigh = sdhigh
)

Inits <- function() list(       #  initial values for all stochastic nodes.  
                mod = round(1,0),
                candbeta = c(3,3,3,3,3,3),
                K0 = c(mean(ts[1:2]),NA,mean(ts[1:2])*(17/15)),
                hlrate2 = hllow+(hlhigh-hllow)/2,
                hlrate3 = hllow+(hlhigh-hllow)/2,
                candDD = 1,
                fracharv2 = 1-exp(mean(betas)),
                constharv2 = abs(slope),
                candharv = 1,
                startyear = rep(startyear,3),
                length = rep(endyear-startyear,3),
                sd = rep((stdev/2),3)
)

Par <- c(
         "beta",
         #"hlmean",
         #"harvratemean",
         #"harvratesd",
         #"harvratemean2",
         #"harvratesd2",
         #"ceilinggrowth",
         #"hlsd", 
         "mod",
         "logN.new",
         "logN.new2",
         "K0",
         "sd", 
         "constharv2",
         "fracharv2",
         "hlrate2",
            #"overK",
         #"ceilinggrowth",
         #"rickergrowth",
            #"constharv3",
         "harvtype",
         "hltype"
           #"fracharv2",
           #"constharv0",
           #"fracharv0",
           #"startyear",
           #"length",
           #"threat_final"
) 

BugFile <- ("tfp9_realdata5_3d_hablossdata.bug")  #  ("tfp9_realdata5_3b.bug")   #3b has logN.new    

Mod <- bugs(data=Data, inits=Inits, parameters.to.save=Par,   #  run BUGS with 3 chains  
                   model.file=BugFile, n.chains=1, n.iter=2000,   # 20000, 10000, 10          # use 5000 iterations, 2500 burnin
                   bugs.directory="C:/Users/Kevin/WinBUGS/winbugs14/WinBUGS14",  #"C:/Users/Kevin/WinBUGS/winbugs14/WinBUGS14",
                   n.burnin=1000,n.thin=2,debug=T)    # ,debug=T




##############################













