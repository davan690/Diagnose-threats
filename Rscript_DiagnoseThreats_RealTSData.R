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



###################################
###################################
     # load previous results

setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\R Workspace")
load("DiagnosisResults_Real.RData")


##########################################################
############               INITIALIZE GLOBAL PARAMETERS


cand_models <- c("null (no threat)","habitat loss","exploitation")


##################################################################################
##################################################################################
#   TEST WITH REAL TIME SERIES

# read in data
setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\DATA")

realhabitat_df <- read.table("realdata_skylark_habloss1.txt",sep="",header=F)   # farmland habitat alteration index...


gentime2 = 1    #  approximate generation time, 
realspecname <- "Atlantic cod"  # "Southern bluefin tuna"  #  "Red grouse"  # "European Skylark"  #   "Thorny skate"  #  
reallocation <- "in Greenland"  #  "in South Pacific" # "in Langholm moor, Scotland"  # "in Great Britain" # "" #  "in North Atlantic"  #   
realsource <- "GPDD, Imperial College"  #  "(Thirgood et al. 2000)"  #  "BTO/JNCC/RSPB Breeding Bird Survey" #  "" #  "NEFSC (www.nefsc.noaa.gov)"  #  
realfilename <- "realdata_cod2.txt"  # "realdata_SBT1.txt"  # "realdata_grouse1.txt" #  "realdata_skylark3.txt"  #  "realdata_skate1.txt" #  "realdata_cod2.txt"  #  "realdata_farmbirds1.txt"  #
realdata <- read.table(realfilename,header=F,sep="") 

COD = T 

 #realdata <- realdata[-which(realdata[,1]%in%c(1940:1948)),]     ### for grouse data, remove years 1940 to 1947
 
startyear = 1
endyear = floor((nrow(realdata))/gentime2)-1  # 15  #   
yaxname = "Abundance index"
realyears = realdata[,1]

ts <- array(0,dim=c(gentime2,floor(nrow(realdata)/gentime2)))  # c(1,30)
scaling_factor <- 100/max(realdata[,2])
for(i in 1:gentime2){
  ts[i,] <- realdata[seq(i,nrow(realdata),gentime2),2] * scaling_factor      # re-scale to approximately 1000
}

nyears2 <- ncol(ts)
ncandidates <- 3
StartTimes <- array(0,dim=nyears2)
StartTimes[startyear] <- 1
Lengths <- array(0,dim=nyears2)
Lengths[(endyear-startyear)] = 1

if(COD){
  Lengths <- array(0,dim=nyears2)
  Lengths[17] = 1
}

realhabitat <- array(0,dim=c(gentime2,nrow(realdata)))
 #realhabitat[1,] <- realhabitat_df[,2]
realhabitat[1,] <- matrix(seq(1,nyears2,1),nrow=gentime2,ncol=nyears2)

realtime <- matrix(seq(1,nyears2,1),nrow=gentime2,ncol=nyears2)

# run BUGS code.
                           # plot out the real time series
par(mai=c(.75,1,.3,.3))
plot(as.vector(ts),type="l",ylim=c(0,max(ts)*1.1),lwd=2,col=gray(.2),xlab="",
              xaxt="n",ylab=yaxname)   # first plot out the original time series time (years), abundance
axis(1,at=c(1:nrow(realdata)),labels=realyears,tick=F)
text(endyear*gentime2-5,95,realsource,cex=.8,col=gray(.3))
text(nrow(realdata)/2,106,paste(realspecname,reallocation),cex=1.2)


year <- seq(1,nyears2,1)
slope <- summary( lm( as.vector(ts)[startyear:(endyear*gentime2)]~as.vector(realhabitat)[startyear:(endyear*gentime2)] ) )$coefficients[2,1]
hlhigh <- slope*1.5
hllow  <- slope/5
nullhigh <- (mean(ts[1:2])*0.1)/nyears2
nulllow <- -1*(mean(ts[1:2])*0.1)/nyears2

#if (abs(hllow)<0.6) hllow <- -0.6 
#if (abs(hlhigh)<1) hlhigh <- -1

stdev <- sd(log(as.vector(ts)))

betas <- log(as.vector(ts))[-1] - log(as.vector(ts))[-length(ts)] 
maxbeta <- max(betas) 
highbeta2 <- maxbeta*2    # for hl model with Ricker DD
highbeta <- maxbeta/2    # maxbeta/2 seems to work well?***
medbeta <- maxbeta/4   # maxbeta/5 seems to work well?***
lowbeta <- maxbeta/10
minbeta <- min(betas[which(betas>0)])
sdlow  <- sd(betas)/1000    # 1000 seems to work well?***
sdhigh  <- sd(betas)*50    # 15 seems to work well?***


############################3
#################################################################################


setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\BUGS files")


                     ######  generic BUGS code
Data <- list( 
             y = ts,
             nyears = nyears2,
             realhabitat = realhabitat,     
             realtime = realtime,
             gentime=gentime2,
             startabun = mean(ts[1:2]),
             nmodels = ncandidates,
             startyear2 = as.vector(StartTimes),
             length2 = as.vector(Lengths),
             hlhigh = hlhigh,
             hllow = hllow,
             nullhigh = nullhigh,
             nulllow = nulllow,
             hrhigh = -1*hlhigh,
             hrlow  = 0.01,  #  -1*hllow,  # 
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
                fracharv2 = 1-exp(mean(betas)),
                constharv2 = abs(slope),
                candharv = 1,
                startyear = rep(startyear,3),
                length = rep(endyear-startyear,3),
                sd = rep((stdev/2),3)
)

Par <- c(
         "beta", 
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
         "harvtype",
         "hltype"
) 

   
BugFile <- ("tfp9_realdata8_mixture.bug")  #BugFile <- ("tfp9_realdata5_3d_hablossdata.bug")  #  ("tfp9_realdata5_3b.bug")   #3b has logN.new 
    
BugDir <-  "C:\\Users\\Kevin\\Documents\\Employment\\ESF\\Bog Turtle\\DATA\\software\\BUGS\\WinBUGS14"   #   "C:\\Users\\Kevin\\Desktop\\Kevin\\WinBUGS14"  #   

Mod <- bugs(data=Data, inits=Inits, parameters.to.save=Par,   #  run BUGS with 3 chains  
                   model.file=BugFile, n.chains=1, n.iter= 5000, #5000,   # 20000, 10000, 10          # use 5000 iterations, 2500 burnin
                   bugs.directory=BugDir,
                   n.burnin=2500,n.thin=2) #,debug=T)    # 2500 ,debug=T



 # Mod_Skylark_AD <- Mod
 # Mod_Skylark_noAD <- Mod
 # Mod_Grouse <- Mod
 # Mod_Tuna <- Mod
 # Mod_Cod <- Mod
 # Mod_Cod2 <- Mod
 

table(Mod$sims.list$mod)


##############################





##################################################
############   SAVE THE BUGS OUTPUT FILES FOR PLOTTING 


setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\R Workspace")
save(Mod_Skylark_AD,
     Mod_Skylark_noAD,
     Mod_Grouse,
     Mod_Tuna,
     Mod_Cod,
     Mod_Cod2,   # with threat ending in 1975
     file="DiagnosisResults_Real.RData"
)

##### NOTE: go to plotting script to use these results....





#############################################
########################   END SCRIPT





