###### script for generating time series data for the discrimination algorithm


########################
#####    SET WD and LOAD PACKAGES


###  first load matt's mp.read function AND metapopversion AND fill.matrix
            #  AND mp.write

  #source([mp.read directory])  # source the scripts from dropbox???

source("C:/Users/Kevin/Dropbox/SACode_Sandbox/mp.write.r")
source("C:/Users/Kevin/Dropbox/SACode_Sandbox/mp.read.r")
source("C:/Users/Kevin/Dropbox/SACode_Sandbox/metapopversion.r")
source("C:/Users/Kevin/Dropbox/SACode_Sandbox/fill.matrix.df.r")


####################
########    SET WORKING DIRECTORY

 #setwd("C:\\Users\\Kevin\\Desktop\\TFP_DATA\\SimulatedTS")
setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\Simulated_TS")

#######################
##### READ IN TEMPLATE METAPOP FILE

               # then, specify directory for mp file
mpfilein <- "Template.mp"
 #mpfilein2 <- "Template_alt.mp"
 #mpfilein3 <- "Template_ceiling.mp"

template <- mp.read(mpfilein)
 #alt_template <- mp.read(mpfilein2)
 #ceil_template <- mp.read(mpfilein3)


# note: initial abundance and K are set at 10000

initAbund <- 10000
initK <- 10000
NReps <- 50

timeFrames <- c(15,30,45)
constHarvRates <- c(150,275,400)   #  c("75 ","150","275","400")
fracHarvRates <-  c(0.0300, 0.0450, 0.0600)   #  c("0.0150 ", "0.0300 ", "0.0450 ", "0.0600 ")
  #fracHarvRates_alt <- c(" 0.0300", " 0.0450", " 0.0600")  #  c(" 0.0150", " 0.0300", " 0.0450", " 0.0600")
habLossRates <- c(50,125,200)    #   c(50,125,200,300)

Variabilities <- c(0.03,0.05,0.07)
VarNames <- c("low","med","high")


 #### NOTE: variables to change
# newMP$mp.file$PopData_df$MaxR     # strength of density dependence: 1.0001 is none, 1.01 is weak, 1.08 is moderate, 1.15 is strong 
# newMP$mp.file$PopData_df$DensDep  # "CE" is ceiling, "LO" is scramble (for "local"??), "EX" is exponential
# newMP$mp.file$PopData_df$KchangeSt  # rate of K change or name of filename (.KCH) specifying carrying capacity over time 




#####################
###### first generate stable habitat models with no threat process and different levels of DD and variability  ...
### can we get away without a random walk model here? The random walk gets in the way frequently, since it can lead to severe declines
       # that have more or less of the signature of one or the other threat categories.  

nullDDLevels <- c(1.05,1.10,1.16)    # different DD levels for null model
nullDDLevelsNames <- c("weak","moderate","strong")

newName <- character(1000)
XYears <- numeric(1000)
XThreatScenario <- numeric(1000)
XThreatName <- character(1000)
XHabLoss <- numeric(1000)
XConstHarv <- numeric(1000)
XFracHarv <- numeric(1000)
XVariability <- numeric(1000)
XDDStrength <- numeric(1000)
XDDType <- character(1000)

counter = 1

i=1;t=1;v=1
for(i in 1:length(nullDDLevels)){
  for (t in 1:length(timeFrames)){
    for(v in 1:length(Variabilities)){
      newMP <- template
      newMP$mp.file$MaxDur <- timeFrames[t]

      newMP$mp.file$MaxRep <- NReps

      newMP$mp.file$PopData_df$InitAbund <- initAbund
      newMP$mp.file$PopList[[1]]$InitAbund <- initAbund

      newMP$mp.file$PopData_df$K <- initK
      newMP$mp.file$PopList[[1]]$K <- initK

      newMP$mp.file$PopData_df$MaxR <- nullDDLevels[i]
      newMP$mp.file$PopList[[1]]$MaxR <- nullDDLevels[i]

      newMP$mp.file$PopData_df$KchangeSt <- 0
      newMP$mp.file$PopList[[1]]$KchangeSt <- 0

      # newMP$mp.file$PopManageProp
      newMP$mp.file$PopManageProp$active <- 0
      newMP$mp.file$PopManageProp$number <- 0

      #substr(newMP$mp.file$PopManageProp,18,19) <- as.character(timeFrames[t])   # time frame of harvest
      #substr(newMP$mp.file$PopManageProp,25,28) <- "0 0 "   # specify fractional harvest... (for null model doesn't matter)
      #substr(newMP$mp.file$PopManageProp,29,31) <- "00 "     # specify constant harvest (shouldn't matter for null model)
      #substr(newMP$mp.file$PopManageProp,32,38) <- "0.0000 "  # specify rate of fractional harvest

      newMP$mp.file$SDMatr[[1]]$Matr[1,1] <- Variabilities[v]

      newName[counter] <- paste("NULL_",nullDDLevelsNames[i],"DD_",timeFrames[t],"years_",VarNames[v],"var.mp",sep="")
      XYears[counter] <- timeFrames[t]
      XThreatScenario[counter] <- 1           # means null
      XThreatName[counter] <- "null"
      XHabLoss[counter] <- 0
      XConstHarv[counter] <- 0
      XFracHarv[counter] <- 0
      XVariability[counter] <- Variabilities[v]  #c(1:length(Variabilities))[v]
      XDDStrength[counter] <- nullDDLevels[i]   # c(0:(length(nullDDLevels)-1))[i]
      XDDType[counter] <- "Ricker"   # ifelse(i==1,"Exponential","Ricker") 

      version <- template$version
      mpfileout <- newName[counter]

      mp.write(newMP$mp.file,version,mpfileout)
      counter = counter + 1
    }
  }
}


####################
###### generate habitat loss models

DDTypes <- c("LO")  #  c("LO","CE")
DDNames <- c("Ricker")  # c("Ricker","Ceiling")
HLDDLevels <- c(1.05,1.10,1.16)  # only use two DD levels for HL model?   # lower Rmax and HL is no longer a major threat due to lack of population response to carrying capacity under the Ricker model  
HLDDLevelsNames <- c("weak","moderate","strong")         # this parameter can be considered strength of response to carrying capacity...

t=1;i=1;d=2;f=1;v=1
for(t in 1:length(timeFrames)){
  for(i in 1:length(habLossRates)){
    ndx <- 0
          # develop a KCH file for this habitat loss rate and this time frame ...
    KCH_name <- paste("HL",habLossRates[i],"_",timeFrames[t],"years.KCH",sep="") 
    KCHvec <- numeric(timeFrames[t])
    for(j in 1:timeFrames[t]){
      KCHvec[j] <- 10000 - habLossRates[i]*(j-1)
      if(KCHvec[j]<600) ndx <- 1
      if(KCHvec[j]<600) break
    }
    if(ndx==1) break
    write.table(data.frame(var=KCHvec),file=KCH_name,col.names=F,row.names=F)

    for(d in 1:length(DDTypes)){
      for(f in 1:length(HLDDLevels)){
        for(v in 1:length(Variabilities)){
          newMP <- template
            #if(DDNames[d]=="Ceiling") newMP <- ceil_template

           # set the time frame

          newMP$mp.file$PopData_df$InitAbund <- initAbund
          newMP$mp.file$PopList[[1]]$InitAbund <- initAbund

          newMP$mp.file$PopData_df$K <- initK
          newMP$mp.file$PopList[[1]]$K <- initK

          newMP$mp.file$MaxDur <- timeFrames[t]

          newMP$mp.file$MaxRep <- NReps
 
           # set RMax  
          if(DDNames[d]=="Ricker") newMP$mp.file$PopData_df$MaxR <- HLDDLevels[f]      
          if(DDNames[d]=="Ricker") newMP$mp.file$PopList[[1]]$MaxR <- HLDDLevels[f]

          if(DDNames[d]=="Ceiling") newMP$mp.file$StMatr[[1]]$Matr[1,1] <- exp(log(HLDDLevels)/2)[f]        

           # add the name of the KCH file
          newMP$mp.file$PopData_df$KchangeSt <- KCH_name
          newMP$mp.file$PopList[[1]]$KchangeSt <- KCH_name

          # make sure no harvest is occurring
          # newMP$mp.file$PopManageProp
          newMP$mp.file$PopManageProp$active <- 0
          newMP$mp.file$PopManageProp$number <- 0

          #substr(newMP$mp.file$PopManageProp,18,19) <- as.character(timeFrames[t])   # time frame of harvest
          #substr(newMP$mp.file$PopManageProp,25,28) <- "0 0 "   # specify fractional harvest... (for null model doesn't matter)
          #substr(newMP$mp.file$PopManageProp,29,31) <- "00 "     # specify constant harvest (shouldn't matter for null model)
          #substr(newMP$mp.file$PopManageProp,32,38) <- "0.0000 "  # specify rate of fractional harvest

          # set the DD Type to the appropriate value
          newMP$mp.file$PopData_df$DensDep  <- DDTypes[d]
          newMP$mp.file$PopList[[1]]$DensDep  <- DDTypes[d]

            #if(DDNames[d]=="Ceiling") newMP$mp.file$StMatr[[1]]$Matr[1,1] <- HLDDLevels[f]

          # set the variability
          newMP$mp.file$SDMatr[[1]]$Matr[1,1] <- Variabilities[v]

          newName[counter] <- paste("HL",habLossRates[i],"_",DDNames[d],"_",HLDDLevelsNames[f],"DD_",timeFrames[t],"years_",VarNames[v],"var.mp",sep="")
          XYears[counter] <- timeFrames[t]
          XThreatScenario[counter] <- 2
          XThreatName[counter] <- "habitat loss"
          XHabLoss[counter] <- habLossRates[i]      #c(1:length(habLossRates))[i]
	    XConstHarv[counter] <- 0
          XFracHarv[counter] <- 0
          XVariability[counter] <- Variabilities[v]   #c(1:length(Variabilities))[v]
          XDDStrength[counter] <- HLDDLevels[f]   
          XDDType[counter] <- DDNames[d] 

          version <- template$version
          mpfileout <- newName[counter]

          mp.write(newMP$mp.file,version,mpfileout)
          counter = counter + 1 
        }
      }   
    }
  }
}



####################
###### generate harvest models

 #HarvTypes <- c("0 0 ","0 1 ")
HarvNames <- c("Const","Frac")

harvDDLevels <- c(1.05,1.10,1.16)    # 2 different DD levels for harvest model (higher DD and harvest is no longer a threat due to strong internal compensation)
harvDDLevelsNames <- c("weak","moderate","strong")

t=1;i=1;v=1;d=1
for(t in 1:length(timeFrames)){

        # first loop through constant harvest scenarios
  for(i in 1:length(constHarvRates)){
    ndx <- 0

    for(v in 1:length(Variabilities)){
      for(d in 1:length(harvDDLevels)){
                  # determine if pop likely to go negative ...
        harvvec <- numeric(timeFrames[t])
          #for(j in 1:timeFrames[t]){
          #  harvvec[j] <- 10000 - as.numeric(constHarvRates[i])*(j-1)
          #  if(harvvec[j]<-100&d==1) ndx <- 1
          #  if(harvvec[j]<-100&d==1) break
          #}
          #if(ndx==1) break

        newMP <- template
         #if(nchar(as.numeric(constHarvRates[i]))==3) newMP <- alt_template

           # set the time frame
        newMP$mp.file$MaxDur <- timeFrames[t]

        newMP$mp.file$MaxRep <- NReps

        newMP$mp.file$PopData_df$InitAbund <- initAbund
        newMP$mp.file$PopList[[1]]$InitAbund <- initAbund

        newMP$mp.file$PopData_df$K <- initK
        newMP$mp.file$PopList[[1]]$K <- initK

 
           # set RMax  
        newMP$mp.file$PopData_df$MaxR <- harvDDLevels[d]      
        newMP$mp.file$PopList[[1]]$MaxR <- harvDDLevels[d]  

           # add the name of the KCH file
        newMP$mp.file$PopData_df$KchangeSt <- 0
        newMP$mp.file$PopList[[1]]$KchangeSt <- 0

          # newMP$mp.file$PopManageProp
        newMP$mp.file$PopManageProp$active <- 1
        newMP$mp.file$PopManageProp$end.time <- timeFrames[t]
        newMP$mp.file$PopManageProp$num.or.prop <- 0
        newMP$mp.file$PopManageProp$number <- constHarvRates[i]
        newMP$mp.file$PopManageProp$proportion <- 0
        
          # set harvest model
        #substr(newMP$mp.file$PopManageProp,18,19) <- as.character(timeFrames[t])   # time frame of harvest
        #substr(newMP$mp.file$PopManageProp,25,28) <- "0 0 "   # specify constant harvest...
        #substr(newMP$mp.file$PopManageProp,29,31) <- constHarvRates[i]     # specify constant harvest rate
        #substr(newMP$mp.file$PopManageProp,32,38) <- "0.0000 "  # specify rate of fractional harvest
        #if(nchar(as.numeric(constHarvRates[i]))==3) substr(newMP$mp.file$PopManageProp,32,38) <- " 0.0000"

          # set the DD Type to the appropriate value
        newMP$mp.file$PopData_df$DensDep  <- "LO"
        newMP$mp.file$PopList[[1]]$DensDep  <- "LO"

          # set the variability
        newMP$mp.file$SDMatr[[1]]$Matr[1,1] <- Variabilities[v]

        newName[counter] <- paste("Constharv",as.numeric(constHarvRates[i]),"_",harvDDLevelsNames[d],"DD_",timeFrames[t],"years_",VarNames[v],"var.mp",sep="")
        XYears[counter] <- timeFrames[t]
        XThreatScenario[counter] <- 3
        XThreatName[counter] <- "harvest"
        XHabLoss[counter] <- 0
        XConstHarv[counter] <- constHarvRates[i]   #c(1:length(constHarvRates))[i]
        XFracHarv[counter] <- 0
        XVariability[counter] <- Variabilities[v]   #  c(1:length(Variabilities))[v]
        XDDStrength[counter] <- nullDDLevels[d]  #  c(0:(length(nullDDLevels)-1))[d]   # set at high DD
        XDDType[counter] <- "Ricker"   #  ifelse(d==1,"Exponential","Ricker")

        version <- template$version
        mpfileout <- newName[counter]

        mp.write(newMP$mp.file,version,mpfileout)
        counter = counter + 1  
      }  
    }
  }

        # then loop through fractional harvest scenarios
  for(i in 1:length(fracHarvRates)){
    ndx <- 0

    for(v in 1:length(Variabilities)){
      for(d in 1:length(harvDDLevels)){
        newMP <- template

           # set the time frame
        newMP$mp.file$MaxDur <- timeFrames[t]
        newMP$mp.file$MaxRep <- NReps

        newMP$mp.file$PopData_df$InitAbund <- initAbund
        newMP$mp.file$PopList[[1]]$InitAbund <- initAbund

        newMP$mp.file$PopData_df$K <- initK
        newMP$mp.file$PopList[[1]]$K <- initK
 
           # set RMax  
        newMP$mp.file$PopData_df$MaxR <- harvDDLevels[d]      
        newMP$mp.file$PopList[[1]]$MaxR <- harvDDLevels[d]  

           # remove KCH file
        newMP$mp.file$PopData_df$KchangeSt <- 0
        newMP$mp.file$PopList[[1]]$KchangeSt <- 0

          # set harvest model

          # newMP$mp.file$PopManageProp
        newMP$mp.file$PopManageProp$active <- 1
        newMP$mp.file$PopManageProp$end.time <- timeFrames[t]
        newMP$mp.file$PopManageProp$num.or.prop <- 1
        newMP$mp.file$PopManageProp$number <- 0
        newMP$mp.file$PopManageProp$proportion <- fracHarvRates[i]

        #substr(newMP$mp.file$PopManageProp,18,19) <- as.character(timeFrames[t])   # time frame of harvest
        #substr(newMP$mp.file$PopManageProp,25,28) <- "0 1 "   # specify fractional harvest...
        #substr(newMP$mp.file$PopManageProp,29,31) <- "0  "     # specify constant harvest rate
        #substr(newMP$mp.file$PopManageProp,32,38) <- fracHarvRates[i]  # specify rate of fractional harvest

          # set the DD Type to the appropriate value
        newMP$mp.file$PopData_df$DensDep  <- "LO"
        newMP$mp.file$PopList[[1]]$DensDep  <- "LO"

          # set the variability
        newMP$mp.file$SDMatr[[1]]$Matr[1,1] <- Variabilities[v]

        newName[counter] <- paste("Fracharv",as.numeric(fracHarvRates[i])*100,"_",harvDDLevelsNames[d],"DD_",timeFrames[t],"years_",VarNames[v],"var.mp",sep="")
        XYears[counter] <- timeFrames[t]
        XThreatScenario[counter] <- 3
        XThreatName[counter] <- "harvest"
        XHabLoss[counter] <- 0
        XConstHarv[counter] <- 0
        XFracHarv[counter] <- fracHarvRates[i]  #  c(1:length(fracHarvRates))[i]
        XVariability[counter] <- Variabilities[v]   #  c(1:length(Variabilities))[v]
        XDDStrength[counter] <- harvDDLevels[d]   #  c(0:(length(harvDDLevels)-1))[d]   # set at high DD
        XDDType[counter] <- "Ricker"  # ifelse(d==1,"Exponential","Ricker")

        version <- template$version
        mpfileout <- newName[counter]

        mp.write(newMP$mp.file,version,mpfileout)
        counter = counter + 1  
      }  
    }
  }
}



############################
####  MAKE COMBINED THREAT SCENARIOS


####################
###### generate habitat loss models

DDTypes <- c("LO")  #  c("LO","CE")
DDNames <- c("Ricker")  # c("Ricker","Ceiling")
HLDDLevels <- c(1.10)  # only use two DD levels for HL model?   # lower Rmax and HL is no longer a major threat due to lack of population response to carrying capacity under the Ricker model  
HLDDLevelsNames <- c("moderate")         # this parameter can be considered strength of response to carrying capacity...
timeFrames <- c(45)
habLossRates <- c(50,200)
   #Variabilities <- c(0.03,0.07)


t=1;i=1;d=2;f=1;h=1;v=1
for(t in 1:length(timeFrames)){   
  for(i in 1:length(habLossRates)){
    ndx <- 0
          # develop a KCH file for this habitat loss rate and this time frame ...
    KCH_name <- paste("HLComb",habLossRates[i],"_",timeFrames[t],"years.KCH",sep="") 
    KCHvec <- numeric(timeFrames[t])
    for(j in 1:timeFrames[t]){
      KCHvec[j] <- 10000 - habLossRates[i]*(j-1)
      if(KCHvec[j]<600) ndx <- 1
      if(KCHvec[j]<600) break
    }
    if(ndx==1) break
    write.table(data.frame(var=KCHvec),file=KCH_name,col.names=F,row.names=F)

    for(d in 1:length(DDTypes)){
      for(f in 1:length(HLDDLevels)){
        for(h in 1:length(fracHarvRates)){
          for(v in 1:length(Variabilities)){
            newMP <- template
              #if(DDNames[d]=="Ceiling") newMP <- ceil_template

             # set the time frame

            newMP$mp.file$PopData_df$InitAbund <- initAbund
            newMP$mp.file$PopList[[1]]$InitAbund <- initAbund

            newMP$mp.file$PopData_df$K <- initK
            newMP$mp.file$PopList[[1]]$K <- initK

            newMP$mp.file$MaxDur <- timeFrames[t]

            newMP$mp.file$MaxRep <- NReps

 
             # set RMax  
            if(DDNames[d]=="Ricker") newMP$mp.file$PopData_df$MaxR <- HLDDLevels[f]      
            if(DDNames[d]=="Ricker") newMP$mp.file$PopList[[1]]$MaxR <- HLDDLevels[f]
  
            if(DDNames[d]=="Ceiling") newMP$mp.file$StMatr[[1]]$Matr[1,1] <- exp(log(HLDDLevels)/2)[f]        

             # add the name of the KCH file
            newMP$mp.file$PopData_df$KchangeSt <- KCH_name
            newMP$mp.file$PopList[[1]]$KchangeSt <- KCH_name

                # newMP$mp.file$PopManageProp
            newMP$mp.file$PopManageProp$active <- 1
            newMP$mp.file$PopManageProp$end.time <- timeFrames[t]
            newMP$mp.file$PopManageProp$num.or.prop <- 1
            newMP$mp.file$PopManageProp$number <- 0
            newMP$mp.file$PopManageProp$proportion <- fracHarvRates[h]

            #substr(newMP$mp.file$PopManageProp,18,19) <- as.character(timeFrames[t])   # time frame of harvest
            #substr(newMP$mp.file$PopManageProp,25,28) <- "0 0 "   # specify fractional harvest... (for null model doesn't matter)
            #substr(newMP$mp.file$PopManageProp,29,31) <- "00 "     # specify constant harvest (shouldn't matter for null model)
            #substr(newMP$mp.file$PopManageProp,32,38) <- "0.0000 "  # specify rate of fractional harvest

            # set the DD Type to the appropriate value
            newMP$mp.file$PopData_df$DensDep  <- DDTypes[d]
            newMP$mp.file$PopList[[1]]$DensDep  <- DDTypes[d]

            #if(DDNames[d]=="Ceiling") newMP$mp.file$StMatr[[1]]$Matr[1,1] <- HLDDLevels[f]

            # set the variability
            newMP$mp.file$SDMatr[[1]]$Matr[1,1] <- Variabilities[v]

            newName[counter] <- paste("Combined_HL",habLossRates[i],"_","Fracharv",as.numeric(fracHarvRates[h])*100,"_",DDNames[d],"_",HLDDLevelsNames[f],"DD_",timeFrames[t],"years_",VarNames[v],"var.mp",sep="")
            XYears[counter] <- timeFrames[t]
            XThreatScenario[counter] <- 4    # combined
            XThreatName[counter] <- "combined"
            XHabLoss[counter] <- habLossRates[i]  #  c(1:length(habLossRates))[i]
            XConstHarv[counter] <- 0
            XFracHarv[counter] <- fracHarvRates[h]
            XVariability[counter] <- Variabilities[v]   # c(1:length(Variabilities))[v]
            XDDStrength[counter] <- HLDDLevels[f]   
            XDDType[counter] <- DDNames[d] 

            version <- template$version
            mpfileout <- newName[counter]

            mp.write(newMP$mp.file,version,mpfileout)
            counter = counter + 1 
          }
        }
      }   
    }
  }
}



#############################
####  WRITE OUT OF FILE NAMES and other data

   # remove null rows 
newName <- newName[1:(counter-1)]
XYears <- XYears[1:(counter-1)]
XThreatScenario <- XThreatScenario[1:(counter-1)]
XThreatName <- XThreatName[1:(counter-1)]
#XSeverity <- XSeverity[1:(counter-1)]
XHabLoss <- XHabLoss[1:(counter-1)]
XConstHarv <- XConstHarv[1:(counter-1)]
XFracHarv <- XFracHarv[1:(counter-1)]
XVariability <- XVariability[1:(counter-1)]
XDDStrength <- XDDStrength[1:(counter-1)]
XDDType <- XDDType[1:(counter-1)]

   ### generate a data frame with all key information for each MP file: number of years, DD type, variability,  

tempdf <- data.frame(filename=newName,nyears=XYears,threat=XThreatScenario,threatname=XThreatName,
               HabLoss=XHabLoss,ConstHarv=XConstHarv,FracHarv=XFracHarv,
                variability=XVariability,DDStrength=XDDStrength,DDType=XDDType)

write.table(tempdf,"mp_filedata.txt",row.names=F,col.names=T)


#############################
####   WRITE OUT BATCH FILE FOR RUNNING MODELS

   ##  structure of each line:   call runraw Lw15-b01

textlines <- character(length(newName))
for(i in 1:length(newName)){
  textlines[i] <- paste("call runraw ",substr(newName[i],1,(nchar(newName[i])-3)),sep="")
}

write.table(data.frame(col1=textlines),file="run.bat",row.names=F,col.names=F,quote=F)


############################

