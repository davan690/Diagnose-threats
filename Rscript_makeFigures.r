
########################################################
########################################################
#####    R Script for plotting results from the TFP project
### K Shoemaker, October 2013




#############
#       READ IN DATA

   #load("BugsOutput.RData")

setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\TFP Results")

load("LatestMasterFile.RData")    # read in name of latest master file

  # ls()

master <- read.csv(LatestMasterFile,header=T)

names(master)


   # read in simulated time series...

#setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\R Workspace")
#load("SimulatedTimeSeries.RData")

#ls()

setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\R Workspace")
load("DiagnosisResults_Simulated.RData")


   #### load real time series data..

load("DiagnosisResults_Real.RData")


#############
#       LOAD PACKAGES

library(Hmisc)
library(ggplot2)
library(lattice)
library(party)
library(RColorBrewer)


#########
#     DEVELOP a "CORRECT" VARIABLE

nmodels = nrow(master)
master$Correct<-0
master$Wrong<-0
for(m in 1:nmodels){
  if(master$threatname[m]!="combined"){
    denom <- (sum(master$veryRight[m])+sum(master$weaklyRight[m])+sum(master$veryWrong[m])+sum(master$weaklyWrong[m]))
    master$Correct[m] <- (sum(master$veryRight[m])+sum(master$weaklyRight[m]))/denom
                             
    master$Wrong[m] <- (sum(master$veryWrong[m])+sum(master$weaklyWrong[m])) / denom
  } else{
    master$Correct[m] <- NA
    master$Wrong[m] <- NA
  }  
}

master$Correct
master$Wrong


#### DEVELOP A "VARCLASS" VARIABLE

lowvarndx <- grep("lowvar",master$filename)
medvarndx <- grep("medvar",master$filename)
highvarndx <- grep("highvar",master$filename)

master$varclass <- 0
master$varclass[lowvarndx] <- 1
master$varclass[medvarndx] <- 2
master$varclass[highvarndx] <- 3


### DEVELOP A "DDCLASS" variable

fac <- factor(master$DDStrength)  #,levels=c("1","2","3"))
levels(fac) = c("1","2","3")
master$DDclass <- as.numeric(fac)

#######################
#        FINALIZE ORGANIZATION

names(master)
predictors <- c("nyears","threatname","varclass","DDclass","severity") 
npredictors <- length(predictors)
predictorNames <- c("Length of Time Series","Type of Threat","Stochastic Process Noise",
                    "Strength of Density Dependence","Threat Severity")

cbind(predictors,predictorNames)


##############################
###########  MAKE DATASET WITH REPLICATE TIME SERIES AS ESSENTIAL UNIT...

### first, determine the percentage decline of each time series...

nreps=50

    # make sure original time series are properly linked to the master data set.
ndx<-match(master$filename,specs$filename)

   # m=1; r=1
perDecline <- array(0,dim=c(nrow(master),nreps))
Correct <- array(0,dim=c(nrow(master),nreps))
trueModel <- array(0,dim=c(nrow(master),nreps))
nYears <- array(0,dim=c(nrow(master),nreps))
ddClass <- array(0,dim=c(nrow(master),nreps))
Severity <- array(0,dim=c(nrow(master),nreps))  
EnvVar <- array(0,dim=c(nrow(master),nreps))    
EnvVar2 <- array(0,dim=c(nrow(master),nreps))
ConstHarv <- array(0,dim=c(nrow(master),nreps))
RelEnvVar <- array(0,dim=c(nrow(master),nreps)) # env var relative to decline...

    ##### START HERE ********************************
       ### do we use threat severity or actual decline?
       ### do we use variability as set, or actual noisiness?
names(master)

   # m=100; r=1
for(m in 1:nmodels){
  for(r in 1:nreps){
    trueModel[m,r] <- specs$threat[ndx[m]]
    timeSeries <- as.numeric(ts_list[[ndx[m]]][[r]])
    xvalues <- c(1:length(timeSeries)) 
    nYears[m,r] <- length(timeSeries)
    model1 <- lm(timeSeries~xvalues)
    model2 <- loess(timeSeries~xvalues)
       # plot(timeSeries)
       # abline(model1)
       # points(fitted(model2),type="l")
    year1 <- as.numeric(predict(model1,newdata=data.frame(xvalues=1))) 
    finalyear <- as.numeric(predict(model1,newdata=data.frame(xvalues=nYears[m,r])))
    perDecline[m,r] <-  round((1-(finalyear/year1))*100,3) 
    Correct[m,r] <- ifelse(bugs.runs[[ndx[m]]]$details$selected[r]==trueModel[m,r],1,0)
    ddClass[m,r] <- master$DDclass[m]
    Severity[m,r] <- master$severity[m]
    EnvVar[m,r] <- sqrt(mean(resid(model2)^2))   # RMSE from locally weighted model...
    EnvVar2[m,r] <- specs$variability[ndx[m]]
    ConstHarv[m,r] <- ifelse(specs$ConstHarv[ndx[m]]>1,1,0)
    denom <- year1-finalyear
    RelEnvVar[m,r] <- ifelse(denom!=0,EnvVar[m,r]/denom,NA)
  }
}

####################
 #####  make new data frame with observations at the replicate level...

df_ind <- data.frame(
  perDecline = as.vector(perDecline), 
  Correct = as.vector(Correct), 
  trueModel = as.vector(trueModel),
  nYears  = as.vector(nYears),
  ddClass = as.vector(ddClass),
  Severity = as.vector(Severity),
  EnvVar = as.vector(EnvVar),
  EnvVar2 = as.vector(EnvVar2),
  ConstHarv = as.vector(ConstHarv),
  RelEnvVar = as.vector(RelEnvVar)
)

head(df_ind)

   ## remove combined threat scenarios
ndx1 <- which(df_ind$trueModel==4)
ndx2 <- which(df_ind$ConstHarv==1)  ## remove constant harvest scenarios
ndx <- sort(union(ndx1,ndx2))
df_ind2 <- df_ind[-ndx,]

head(df_ind2)  
 
df_ind <- df_ind2     # df_ind$trueModel

#######################
#        FINALIZE ORGANIZATION of per-trajectory data frame

names(df_ind)
predictors_ind <- c("nYears","trueModel","EnvVar","ddClass","perDecline") 
npredictors_ind <- length(predictors)
predictorNames_ind <- c("Length of Time Series","Type of Threat","Variability of Realized Time Series",
                    "Strength of Density Dependence","Linear Decline of Realized Time Series")

cbind(predictors_ind,predictorNames_ind)



#############################
###########  NEW FIGURE: plot diagnostic performance against the degree of decline of the real time series...



### get mean percent decline

ndx1 <- which(master$threat==4)    # remove combined threat scenarios...
ndx2 <- grep("Constharv",master$filename)     # remove constant harvest scenarios
ndx <- union(ndx1,ndx2)   # all rows to be removed
master <- master[-ndx,]    # master$filename

head(master)

meanPerDecline <- apply(perDecline[-ndx,],1,mean)

perCorrect <- round(master$Correct*100,3)
threatName <- master$threatname
ny <- master$nyears
dd <- master$DDclass

df1 <- data.frame(pd = meanPerDecline, pc = perCorrect, tn = threatName, 
                   ny = ny, dd = dd )


########## make figure
setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\TFP Results\\Figures")
graphics.off()

filename = "PerDecline_vs_perRight.svg"

svg(filename = filename,
    width = 8, height = 4, pointsize = 10,
    onefile = TRUE, family = "sans", bg = "white")

par(mfrow=c(1,3))

  # y=1
pch = c(1,3,20)
timeframes <- unique(df1$ny)
threats <- as.character(unique(df1$tn))

  # y=1
for(y in 1:3){
  sub <- subset(df1,ny==timeframes[y])
  sub_hl <- subset(sub,tn==threats[2])
  sub_harv <- subset(sub,tn==threats[3])
  sub_null <- subset(sub,tn==threats[1])

  plot(sub_hl$pc~sub_hl$pd,main=paste(timeframes[y],"Year Time Frame",sep=""),
          xlab="% Decline", ylab="% Correct", pch=pch[1], ylim=c(min(sub$pc),100),
          xlim=c(0,max(sub$pd)) )
  points(sub_harv$pd,sub_harv$pc,pch=pch[2])
  points(sub_null$pd,sub_null$pc,pch=pch[3])

  legend(max(sub$pd)-20,15,pch=pch,legend=threats[c(2,3,1)],bty="n")
}

dev.off()



##########################
######    NEW FIGURE: plot out with trajectory as the essential observation unit...

head(df_ind)


########## make figure


setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\TFP Results\\Figures")
graphics.off()

filename = "PerDecline_vs_perRight_NEW.svg"       # this time at the individual trajectory level
svg(filename = filename,
    width = 8, height = 4, pointsize = 10,
    onefile = TRUE, family = "sans", bg = "white")

par(mfrow=c(1,3))

  # y=1
pch = c(1,3,20)
timeframes <- unique(df1$ny)
threats <- as.character(unique(master$threatname))
lty = c(1,2,3)

legx <- c(40,55,55)
legy <- c(.9,.4,.5)


  # y=1
for(y in 1:3){
  sub <- subset(df_ind,nYears==timeframes[y])
  sub_hl <- subset(sub,trueModel==2)
  sub_harv <- subset(sub,trueModel==3)
  sub_null <- subset(sub,trueModel==1)

     ## generate bins for percent decline...
  declineBins <- c(min(sub$perDecline),seq(0,100,by=10))
 
  pc2 <- array(0,dim=c(11,3))    # percent correct for each bin
   # b=1
  for(b in 1:11){
    pc2[b,1] <- sum(sub_null$Correct[which((sub_null$perDecline>=declineBins[b])&(sub_null$perDecline<=declineBins[(b+1)]))])/
                    length(which((sub_null$perDecline>=declineBins[b])&(sub_null$perDecline<=declineBins[(b+1)]))) 
    pc2[b,2] <- sum(sub_hl$Correct[which((sub_hl$perDecline>=declineBins[b])&(sub_hl$perDecline<=declineBins[(b+1)]))])/
                    length(which((sub_hl$perDecline>=declineBins[b])&(sub_hl$perDecline<=declineBins[(b+1)]))) 

    pc2[b,3] <- sum(sub_harv$Correct[which((sub_harv$perDecline>=declineBins[b])&(sub_harv$perDecline<=declineBins[(b+1)]))])/
                    length(which((sub_harv$perDecline>=declineBins[b])&(sub_harv$perDecline<=declineBins[(b+1)]))) 

  }  
  xvals <- declineBins[2:12]-5
  plot(pc2[,2]~xvals,main=paste(timeframes[y]," Year Time Frame",sep=""),
          xlab="% Decline", ylab="Proportion Correct", pch=pch[1], ylim=c(0,1),
          xlim=c(0,100) )
  points(xvals,pc2[,3],pch=pch[2])
  points(xvals,pc2[,1],pch=pch[3])

      #### lines: based on logistic regression... 

  test <- glm(sub_hl$Correct ~ sub_hl$perDecline, family = binomial)
  o <- order(sub_hl$perDecline)
  lines(sub_hl$perDecline[o], fitted(test)[o],lty=lty[1])

  test <- glm(sub_harv$Correct ~ sub_harv$perDecline, family = binomial)
  o <- order(sub_harv$perDecline)
  lines(sub_harv$perDecline[o], fitted(test)[o],lty=lty[2]) 

  test <- glm(sub_null$Correct ~ sub_null$perDecline, family = binomial)
  o <- order(sub_null$perDecline)
  lines(sub_null$perDecline[o], fitted(test)[o],lty=lty[3])     

  legend(legx[y],legy[y],pch=pch,lty=lty,legend=threats[c(2,3,1)],bty="n")
}

dev.off()



#############################
##############  NEW FIGURE: ASSESS DIAGNOSTIC SUCCESS AS A FUNCTION OF TIME SERIES "NOISE"

  ### note, to make this figure make sense, only time series with decline >=25% are included here...

names(df_ind)

df_ind$EnvVar3 <- (df_ind$EnvVar/10000)*100   
 # df_ind$EnvVar3

########## make figure


setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\TFP Results\\Figures")
graphics.off()

filename = "PerDecline_vs_Noise_NEW.svg"       # this time at the individual trajectory level
svg(filename = filename,
    width = 8, height = 4, pointsize = 10,
    onefile = TRUE, family = "sans", bg = "white")

par(mfrow=c(1,3))

    # y=1
pch = c(1,3,20)
timeframes <- unique(df1$ny)
threats <- as.character(unique(master$threatname))
lty = c(1,2,3)

legx <- c(40,55,10)   # change...
legy <- c(.6,.8,.3)

  # y=1
for(y in 1:3){
  sub <- subset(df_ind,(nYears==timeframes[y])&(perDecline>=25))
  sub_hl <- subset(sub,trueModel==2)
  sub_harv <- subset(sub,trueModel==3)
  sub_null <- subset(sub,trueModel==1)

     ## generate bins for noise as percent of initial abundance...
  EnvVarBins <- c(min(sub$EnvVar3)-.1,seq(1,5,length=10),max(sub$EnvVar3)+1)
 
  pc2 <- array(0,dim=c(11,3))    # percent correct for each bin
   # b=1
  for(b in 1:11){
    pc2[b,1] <- sum(sub_null$Correct[which((sub_null$EnvVar3>=EnvVarBins[b])&(sub_null$EnvVar3<=EnvVarBins[(b+1)]))])/
                    length(which((sub_null$EnvVar3>=EnvVarBins[b])&(sub_null$EnvVar3<=EnvVarBins[(b+1)]))) 
    pc2[b,2] <- sum(sub_hl$Correct[which((sub_hl$EnvVar3>=EnvVarBins[b])&(sub_hl$EnvVar3<=EnvVarBins[(b+1)]))])/
                    length(which((sub_hl$EnvVar3>=EnvVarBins[b])&(sub_hl$EnvVar3<=EnvVarBins[(b+1)]))) 

    pc2[b,3] <- sum(sub_harv$Correct[which((sub_harv$EnvVar3>=EnvVarBins[b])&(sub_harv$EnvVar3<=EnvVarBins[(b+1)]))])/
                    length(which((sub_harv$EnvVar3>=EnvVarBins[b])&(sub_harv$EnvVar3<=EnvVarBins[(b+1)]))) 

  }  
  xvals <- EnvVarBins[2:12]-diff(EnvVarBins)/2
  plot(pc2[,2]~xvals,main=paste(timeframes[y]," Year Time Frame",sep=""),
          xlab="Variability of Realized Time Series", ylab="Proportion Correct", pch=pch[1], ylim=c(0,1),
          xlim=c(0,10) )
  points(xvals,pc2[,3],pch=pch[2])
   #points(xvals,pc2[,1],pch=pch[3])

      #### lines: based on logistic regression... 

  test <- glm(sub_hl$Correct ~ sub_hl$EnvVar3, family = binomial)
  o <- order(sub_hl$EnvVar3)
  lines(sub_hl$EnvVar3[o], fitted(test)[o],lty=lty[1])

  test <- glm(sub_harv$Correct ~ sub_harv$EnvVar3, family = binomial)
  o <- order(sub_harv$EnvVar3)
  lines(sub_harv$EnvVar3[o], fitted(test)[o],lty=lty[2]) 

  #test <- glm(sub_null$Correct ~ sub_null$EnvVar3, family = binomial)
  #o <- order(sub_null$EnvVar3)
  #lines(sub_null$EnvVar3[o], fitted(test)[o],lty=lty[3])     

  legend(legx[y],legy[y],pch=pch[1:2],lty=lty[1:2],legend=threats[c(2,3)],bty="n")
}

dev.off()





##############################
#################  DEVELOP TREE FIGURES 

  #predndx <- #   
response <- "Correct"    
plotname =   "Figure_DiagnosisSuccessTree"  #  "Figure_isOscTree"  #  "Figure_MinPDTree" #    
             #"Figure_ExtTreeAll"  # "Figure_ExtTreePD"  # "Figure_ExtTreeOsc"  # 
             #"Figure_emaTreePD"  #     
SVG=T
PDF=F


predictors2 <-  paste(predictors,collapse="+")

formula <- eval(parse(text=paste(response,"~",predictors2,sep="")))

ndx <- c(1:nrow(master))

subset1 <- master

con <- ctree_control(maxdepth = 4)
treeobj <- ctree(formula, control=con,data=subset1)


DEV=!(SVG|PDF)


############  MAKE PLOT

setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\TFP Results\\Figures")
graphics.off()

  #title <- predictorNames[which(predictors==response)]

if(SVG){
  filename <- paste(plotname,".svg",sep="")
   #dev.off()
  par(mfrow=c(1,1))
  svg(filename = filename,
      width = 7.3, height = 6, pointsize = 10,
      onefile = TRUE, family = "sans", bg = "white")
}

if(PDF){
  filename <- paste(plotname,".pdf",sep="")
  pdf(file = filename,
     width = 7.3, height = 6, pointsize = 500,
     onefile = TRUE, family = "sans", bg = "white")
}


plot(treeobj,main="Prob of Correct Diagnosis",type="simple" , inner_panel=node_inner(treeobj, digits = 2, abbreviate = FALSE, 
  fill = "white", pval = FALSE, id = FALSE),terminal_panel=node_terminal(treeobj, digits = 2, abbreviate = FALSE, 
  fill = c("lightgray", "white"), id = FALSE))    #  #   # type="extended"

dev.off()
graphics.off()


####### IMPORTANCE VALUES FROM CFOREST   ###########

cforest.controls <- cforest_unbiased(ntree=1000,mtry=3)

cforest.controls@fraction=0.4
cforest.controls@splitctrl@minprob=0.05

cforest.fit <- cforest(formula=formula, controls=cforest.controls, data=subset1)

#replace spaces with _
#names(oc)<-gsub(" ", "_",names(oc))

print(cforest.fit)
data.cforest.varimp<-varimp(cforest.fit, conditional= FALSE, OOB=TRUE) #try conditional = TRUE
print(data.cforest.varimp)

#Variables can be considered informative and important if their variable importance value is above the absolute value of the lowest negative-scoring variable

filename2 <- paste(plotname,"ImpVal.svg",sep="")

if(SVG){
  graphics.off()
  svg(filename = filename2,
      width = 7.3, height = 6, pointsize = 15,
      onefile = TRUE, family = "sans", bg = "white")
}

lengthndx <- length(data.cforest.varimp)
par(mai=c(1,3.5,0.2,0.2))
#par(mai=c(1.4,3.4,0.6,0.9))
col <- gray(seq(.7,.2,length=lengthndx))      # rep(brewer.pal(6,"Blues"),each=2)
barplot(height=data.cforest.varimp[order(data.cforest.varimp,decreasing = FALSE)],
        horiz=T,las=1,main="",
        xlab="Index of overall importance",col=col,           
        #names.arg=names(data.cforest.varimp)[order(data.cforest.varimp,decreasing = FALSE)])
        names.arg=predictorNames[match(names(data.cforest.varimp),predictors)][order(data.cforest.varimp,decreasing = FALSE)])

dev.off()


#############################  RANDOM FOREST
#############################  Display univariate plots
  

graphics.off()

if(SVG){
  graphics.off()
  svg(filename = "Figure_ThreatDiag_UnivPlots.svg",
      width = 7.3, height = 6, pointsize = 15,
      onefile = TRUE, family = "sans", bg = "white")
}

RF_UnivariatePlots(object=cforest.fit, varimp=data.cforest.varimp, data=subset1, 
                   predictors=predictors,allpredictors=predictors,labels=predictorNames,plot.layout=c(2,3))

dev.off()  
graphics.off()


####################################
#######################   RANDOM FOREST FIND AND PLOT INTERACTIONS

          # NOTE: this one can take a very long time- maybe up to 10 minutes...
rf_findint <- RF_FindInteractions(object=cforest.fit,data=subset1,predictors=predictors)

    # display and plot out interactions...
rf_findint$interactions1

rf_findint$interactions2

rf_findint$rank.list1

rf_findint$rank.list2

  ### plot interaction strength
graphics.off()
lengthndx <- min(9,nrow(rf_findint$rank.list1))
par(mai=c(0.95,3.1,0.6,0.4))
barplot(height=(rf_findint$rank.list1[c(1:min(9,nrow(rf_findint$rank.list1))),5][c(lengthndx:1)]),
            horiz=T,las=1,main=paste(response, sep=""),
            xlab="Index of interaction strength",col=brewer.pal(lengthndx,"Blues"),           
            names.arg=paste("",predictorNames[match(rf_findint$rank.list1[,2][c(lengthndx:1)],predictors)],"\n",predictorNames[match(rf_findint$rank.list1[,4][c(lengthndx:1)],predictors)],sep="") )


rf_findint$rank.list1


graphics.off()

filename3 = paste(plotname,"IntPlot2.svg",sep="")
svg(filename = filename3,
    width = 6, height = 5, pointsize = 15,
    onefile = TRUE, family = "sans", bg = "white")

#par(mai=c(0.9,1.4,0,1))


# graphics.off()
fam <- "gaussian"  #  "bernoulli"  #
RF_InteractionPlots(1,3,object=cforest.fit,data=subset1,predictors=predictors,family=fam) # zlim=c(0,0.5) 

dev.off()
graphics.off()







#########################################################
#############################################
#   NEXT SERIES OF PLOTS: RANDOM FOREST AT THE INDIVIDUAL TRAJECTORY LEVEL

 ##### NOTE: THIS SECTION IS UNFINISHED AND MAY NOT BE NECESSARY... ****************

names(df_ind)
 ## df_ind$Correct

  #predndx <- #   
response <- "Correct"    
plotname =   "Figure_DiagnosisSuccessTree_binaryResponse"       
SVG=T
PDF=F


predictors2 <-  paste(predictors,collapse="+")

formula <- eval(parse(text=paste(response,"~",predictors2,sep="")))

ndx <- c(1:nrow(master))

subset1 <- master

con <- ctree_control(maxdepth = 4)
treeobj <- ctree(formula, control=con,data=subset1)


DEV=!(SVG|PDF)


############  MAKE PLOT

setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\TFP Results\\Figures")
graphics.off()

  #title <- predictorNames[which(predictors==response)]

if(SVG){
  filename <- paste(plotname,".svg",sep="")
   #dev.off()
  par(mfrow=c(1,1))
  svg(filename = filename,
      width = 7.3, height = 6, pointsize = 10,
      onefile = TRUE, family = "sans", bg = "white")
}

if(PDF){
  filename <- paste(plotname,".pdf",sep="")
  pdf(file = filename,
     width = 7.3, height = 6, pointsize = 500,
     onefile = TRUE, family = "sans", bg = "white")
}


plot(treeobj,main="Prob of Correct Diagnosis",type="simple" , inner_panel=node_inner(treeobj, digits = 2, abbreviate = FALSE, 
  fill = "white", pval = FALSE, id = FALSE),terminal_panel=node_terminal(treeobj, digits = 2, abbreviate = FALSE, 
  fill = c("lightgray", "white"), id = FALSE))    #  #   # type="extended"

dev.off()
graphics.off()


####### IMPORTANCE VALUES FROM CFOREST   ###########

cforest.controls <- cforest_unbiased(ntree=1000,mtry=3)

cforest.controls@fraction=0.4
cforest.controls@splitctrl@minprob=0.05

cforest.fit <- cforest(formula=formula, controls=cforest.controls, data=subset1)

#replace spaces with _
#names(oc)<-gsub(" ", "_",names(oc))

print(cforest.fit)
data.cforest.varimp<-varimp(cforest.fit, conditional= FALSE, OOB=TRUE) #try conditional = TRUE
print(data.cforest.varimp)

#Variables can be considered informative and important if their variable importance value is above the absolute value of the lowest negative-scoring variable

filename2 <- paste(plotname,"ImpVal.svg",sep="")

if(SVG){
  graphics.off()
  svg(filename = filename2,
      width = 7.3, height = 6, pointsize = 15,
      onefile = TRUE, family = "sans", bg = "white")
}

lengthndx <- length(data.cforest.varimp)
par(mai=c(1,3.5,0.2,0.2))
#par(mai=c(1.4,3.4,0.6,0.9))
col <- gray(seq(.7,.2,length=lengthndx))      # rep(brewer.pal(6,"Blues"),each=2)
barplot(height=data.cforest.varimp[order(data.cforest.varimp,decreasing = FALSE)],
        horiz=T,las=1,main="",
        xlab="Index of overall importance",col=col,           
        #names.arg=names(data.cforest.varimp)[order(data.cforest.varimp,decreasing = FALSE)])
        names.arg=predictorNames[match(names(data.cforest.varimp),predictors)][order(data.cforest.varimp,decreasing = FALSE)])

dev.off()


#############################  RANDOM FOREST
#############################  Display univariate plots
  
setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\TFP Results\\Figures")
graphics.off()

if(SVG){
  graphics.off()
  svg(filename = "Figure_ThreatDiag_UnivPlots.svg",
      width = 7.3, height = 6, pointsize = 15,
      onefile = TRUE, family = "sans", bg = "white")
}

RF_UnivariatePlots(object=cforest.fit, varimp=data.cforest.varimp, data=subset1, 
                   predictors=predictors,allpredictors=predictors,labels=predictorNames,plot.layout=c(2,3))

dev.off()  
graphics.off()


####################################
#######################   RANDOM FOREST FIND AND PLOT INTERACTIONS

          # NOTE: this one can take a very long time- maybe up to 10 minutes...
rf_findint <- RF_FindInteractions(object=cforest.fit,data=subset1,predictors=predictors)

    # display and plot out interactions...
rf_findint$interactions1

rf_findint$interactions2

rf_findint$rank.list1

rf_findint$rank.list2

  ### plot interaction strength
graphics.off()
lengthndx <- min(9,nrow(rf_findint$rank.list1))
par(mai=c(0.95,3.1,0.6,0.4))
barplot(height=(rf_findint$rank.list1[c(1:min(9,nrow(rf_findint$rank.list1))),5][c(lengthndx:1)]),
            horiz=T,las=1,main=paste(response, sep=""),
            xlab="Index of interaction strength",col=brewer.pal(lengthndx,"Blues"),           
            names.arg=paste("",predictorNames[match(rf_findint$rank.list1[,2][c(lengthndx:1)],predictors)],"\n",predictorNames[match(rf_findint$rank.list1[,4][c(lengthndx:1)],predictors)],sep="") )


rf_findint$rank.list1


graphics.off()

filename3 = paste(plotname,"IntPlot2.svg",sep="")
svg(filename = filename3,
    width = 6, height = 5, pointsize = 15,
    onefile = TRUE, family = "sans", bg = "white")

#par(mai=c(0.9,1.4,0,1))


# graphics.off()
fam <- "gaussian"  #  "bernoulli"  #
RF_InteractionPlots(4,3,object=cforest.fit,data=subset1,predictors=predictors,family=fam) # zlim=c(0,0.5) 

dev.off()
graphics.off()







###################################################################
################################################
#############
#  NEXT PLOT: examine the effect of Density dependence on the rate of correct diagnosis...

names(subset1)
predictors

setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\TFP Results\\Figures")
graphics.off()

if(SVG){
  graphics.off()
  svg(filename = "Figure_Univar_SensAnal.svg",
      width = 8, height = 8, pointsize = 15,
      onefile = TRUE, family = "sans", bg = "white")
}


par(mfrow=c(3,3))
ndx <- which(subset1$threat==1)
plot(Correct~as.factor(DDclass),data=subset1[ndx,],xlab="DD Strength",ylab="Proportion correct",
                 main="No Threat")
ndx <- which(subset1$threat==2)
plot(Correct~as.factor(DDclass),data=subset1[ndx,],xlab="DD Strength",ylab="Proportion correct",
                 main="Habitat Loss")
ndx <- which(subset1$threat==3)
plot(Correct~as.factor(DDclass),data=subset1[ndx,],xlab="DD Strength",ylab="Proportion correct",
                 main="Exploitation")



################################################
#############
#  NEXT PLOT: examine the effect of variability on the rate of correct diagnosis...

names(subset1)
predictors


  #par(mfrow=c(2,2))
ndx <- which(subset1$threat==1)
plot(Correct~as.factor(varclass),data=subset1[ndx,],xlab="Env. Process Noise",ylab="Proportion correct",
                 main="No Threat")
ndx <- which(subset1$threat==2)
plot(Correct~as.factor(varclass),data=subset1[ndx,],xlab="Env. Process Noise",ylab="Proportion correct",
                 main="Habitat Loss")
ndx <- which(subset1$threat==3)
plot(Correct~as.factor(varclass),data=subset1[ndx,],xlab="Env. Process Noise",ylab="Proportion correct",
                 main="Exploitation")



################################################
#############
#  NEXT PLOT: examine the effect of time series length on the rate of correct diagnosis...

names(subset1)
predictors


 #par(mfrow=c(2,2))
ndx <- which(subset1$threat==1)
plot(Correct~as.factor(nyears),data=subset1[ndx,],xlab="Time Series Duration",ylab="Proportion correct",
                 main="No Threat")
ndx <- which(subset1$threat==2)
plot(Correct~as.factor(nyears),data=subset1[ndx,],xlab="Time Series Duration",ylab="Proportion correct",
                 main="Habitat Loss")
ndx <- which(subset1$threat==3)
plot(Correct~as.factor(nyears),data=subset1[ndx,],xlab="Time Series Duration",ylab="Proportion correct",
                 main="Exploitation")


dev.off()


graphics.off()











###########################################################################################################################
########################################################################
########################################################################
####################################################    PLOTS TO ASSESS GOODNESS OF FIT AND MODEL PERFORMANCE
# test the models    check model fit

REAL = T

   ### select a model and load the information...
  #names(specs)

filenames <- as.character(specs$filename)

setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\DUMP")


   ### HABLOSS EXAMPLE
# target <- 77            ## 77: moderate habitat loss, med variability, 30 years
# filename <- list.files()[grep(as.character(specs$filename[target]),list.files())][1]

   ### HARVEST EXAMPLE
# target <- 203            ## 203: moderate frac harvest, med variability, 30 years
# filename <- list.files()[grep(as.character(specs$filename[target]),list.files())][1]

   ### NULL EXAMPLE
# target <- 14            ## 14: no threat, med variability, 30 years
# filename <- list.files()[grep(as.character(specs$filename[target]),list.files())][1]

if(!REAL){
  nyears2 <- specs$nyears[target]
  gentime2 <- 1

  load(filename)
  ts <- ts_list[[target]][[nreps]] *(100/max(ts_list[[target]][[nreps]]))
}

##################################
#########  Alternatively, load results from a real time series run...

if(REAL){
  Mod <- Mod_Cod2 # Mod_Cod #  Mod_Tuna # Mod_Grouse # Mod_Skylark_AD # Mod_Skylark_noAD # 

  setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\DATA")
  gentime2 = 1    #  approximate generation time, 
  realfilename <- "realdata_cod2.txt"  #  "realdata_SBT1.txt"  #  "realdata_grouse1.txt" #  "realdata_skylark3.txt"  #   "realdata_skate1.txt" #   "realdata_farmbirds1.txt"  #
  realdata <- read.table(realfilename,header=F,sep="") 

  ts <- array(0,dim=c(gentime2,floor(nrow(realdata)/gentime2)))  # c(1,30)
  scaling_factor <- 100/max(realdata[,2])
  for(i in 1:gentime2){
    ts[i,] <- realdata[seq(i,nrow(realdata),gentime2),2] * scaling_factor      # re-scale to approximately 100
  }

  nyears2 <- ncol(ts)
}

setwd("C:\\Users\\Kevin\\Dropbox\\TFP stuff\\TFP Results\\Figures")
graphics.off()

svg(filename = "SimPlot_Cod2.svg", # "SimPlot_Cod.svg", # "SimPlot_Tuna.svg", # "SimPlot_Grouse.svg", # "SimPlot_Skylark_AD.svg", # "SimPlot_Skylark_noAD.svg",   # "SimPlot_nullExample.svg", # "SimPlot_harvestExample.svg",# "SimPlot_hablossExample.svg", #       
    width = 3, height = 4, pointsize = 15,
    onefile = TRUE, family = "sans", bg = "white")



nsims = length(Mod$sims.list$mod)
plottype = "pop size"    #    "growth rate"   #              

cand_models <- c("Null (no threat)","Habitat loss","Exploitation")

if(plottype=="pop size") par(mfcol = c(3,1))   # was mfrow=c(3,1)  *******

harvtype = as.numeric(names(which.max(table(Mod$sims.list$harvtype))))
hltype = as.numeric(names(which.max(table(Mod$sims.list$hltype))))
gentime = c(rep(1,times=length(specs$filename)))

 ## plot=1
for(plot in 1:3){

plotmodel = plot   #  1      #  which mechanism should we plot???

      ##########  plot out models accounting for generation time structure.

tempnames <- as.numeric(names(table(Mod$sims.list$mod)))
postweights <- numeric(3)
for(p2 in 1:3){
  if(p2%in%tempnames) postweights[p2] <- round((as.numeric(table(Mod$sims.list$mod))[which(tempnames==p2)]/nsims)*100,3)  # round((as.numeric(table(Mod$sims.list$mod))/nsims)[p2]*100,3)
}
  
if(plottype == "pop size"){                                                   #  initialize vectors for storing the time series at each year
  ts2 <- numeric(nyears2*gentime2) 
  m1 <- array(0,dim=c(nsims,nyears2*gentime2))      # predicted time series, model 1, null
  m2 <- array(0,dim=c(nsims,nyears2*gentime2))       # predicted time series, model 2, habitat loss
  m3 <- array(0,dim=c(nsims,nyears2*gentime2))        # predicted time series, model 3, harvest
}

if(plottype == "growth rate"){
  ts2 <- numeric((nyears2-1)*gentime2)         # initialize vectors for storing the growth rate at each year
  m1 <- array(0,dim=c(nsims,(nyears2-1)*gentime2))
  m2 <- array(0,dim=c(nsims,(nyears2-1)*gentime2))
  m3 <- array(0,dim=c(nsims,(nyears2-1)*gentime2))
}

if(plottype == "growth rate"){
  counter <- 0
  for(t in 1:nyears2){
    for(g in 1:gentime2){
      if(t>1){ 
        counter <- counter+1
        ts2[counter] <- ts[g,t]/ts[g,t-1] 
      }
    }
  }
}

if(plottype == "pop size"){
  counter <- 0
  for(t in 1:nyears2){
    for(g in 1:gentime2){
      counter = counter+1
      ts2[counter] <- ts[g,t]   #ts_list[[target]][[nreps]][g,t]
    }
  }
}

if(plot==1&plottype=="pop size") par(mai=c(.1,.4,.4,.05))
if(plot==2&plottype=="pop size") par(mai=c(.25,.4,.25,.05))
if(plot==3&plottype=="pop size") par(mai=c(.4,.4,.1,.05))

if(plottype=="growth rate"&plot==1) par(mai=c(.7,.7,.05,.05))

                                           # plot the real time series...
if(plottype=="pop size") par(cex.axis=1.2)
if(plottype=="pop size") plot(ts2,type="b",ylim=c(0,max(ts2)*2.3),yaxt="n",xlab="",ylab="",lwd=2,col=gray(.5))   # first plot out the original time series time (years), abundance

if(!REAL) axis(2,at=c(0,25,50,75,100,125,200),labels=c("0",NA,"5K",NA,"10K",NA,"20K"))
if(REAL) axis(2,at=c(0,25,50,75,100,125,200),labels=c("0",NA,"50",NA,"100",NA,"200"))

#if(plottype=="growth rate"&plot==1)  plot(ts2~c((gentime2+1):(nyears2*gentime2)),pch=1,ylab="",xlab="",ylim=c(.8,1.2),cex=1.5)
if(plottype=="growth rate"&plot==1)  plot(ts2~c((gentime2+1):(nyears2*gentime2)),pch=1,ylab="",xlab="",ylim=c(.8,1.2),cex=1.5)

  # i=1; t=1; g=1
for(i in 1:nsims){
  counter <- 0
   #counter2 <- 0
  for(t in 1:nyears2){
    for(g in 1:gentime2){
      if(plottype=="pop size"){
        counter = counter+1
        if(t==1) m1[i,counter] <- ts[g,t]      # year 1 is always the same...
        if(t==1) m2[i,counter] <- ts[g,t]
        if(t==1) m3[i,counter] <- ts[g,t]

        if(t>1) m2[i,counter] <- exp(Mod$sims.list$logN.new[i,g,t-1,2]) 
        if(t>1) m3[i,counter] <- exp(Mod$sims.list$logN.new[i,g,t-1,3])
        if(t>1) m1[i,counter] <- exp(Mod$sims.list$logN.new[i,g,t-1,1])
      }

      if(plottype=="growth rate"){
        if(t>1){
          counter <- counter+1
          m1[i,counter] <- exp(Mod$sims.list$logN.new2[i,g,t-1,1])/ts[g,t-1]
          m2[i,counter] <- exp(Mod$sims.list$logN.new2[i,g,t-1,2])/ts[g,t-1]
          m3[i,counter] <- exp(Mod$sims.list$logN.new2[i,g,t-1,3])/ts[g,t-1]
        }
      }
    }
  }
  if(plottype == "pop size"){
    if(i%%5 == 0 & plotmodel==2) lines(c(1:counter),m2[i,],col=gray(.8),lty=1,lwd=1.2)  # plot the expected time series under linear trend
    if(i%%5 == 0 & plotmodel==3) lines(c(1:counter),m3[i,],col=gray(.8),lty=1,lwd=1.2)  # plot the expected time series under exponential trend 
    if(i%%5 == 0 & plotmodel==1) lines(c(1:counter),m1[i,],col=gray(.8),lty=1,lwd=1.2)  # plot the expected time series under stable pop.  
  }
  if(plottype == "growth rate"){
    #if(i%%5 == 0 & plotmodel==2) points(c(2:(counter+1)),m2[i,],col=gray(.8),pch=20,cex=.5)  # plot the expected time series under linear trend
    #if(i%%5 == 0 & plotmodel==3) points(c(2:(counter+1)),m3[i,],col=gray(.8),pch=20,cex=.5)  # plot the expected time series under exponential trend 
    #if(i%%5 == 0 & plotmodel==1) points(c(2:(counter+1)),m1[i,],col=gray(.8),pch=20,cex=.5)  # plot the expected time series under stable pop.  
  }
}


            ##################### graphical parameters
     #selected = which.max(modselect_array[m,r,])
linetype <- c(1,3,5)
linewd <- c(2.5,1.5,1.5)
density <- c(10,10,10)
color <- c(gray(.4),gray(.4),gray(.4))
angle <- c(90,90,90)

table(Mod$sims.list$mod)
                         ###############################  MECHANISMS
                                                                     #############
                                                                        # MODEL 1: CONSTANT K
#if(plotmodel==1 & plottype=="pop size"){
#  abline(h=mean(Mod$sims.list$K0[,1]),col=gray(.1),lwd=linewd[1],lty=linetype[1])     # startabun[m,r]
#}

                                                                     #############
                                                                        # MODEL 2: RANDOM WALK (don't show for now)                                                                     
     

                                                                     #############
                                                                        # MODEL 2: HABITAT LOSS

#if((plotmodel == 2 & plottype=="pop size") & hltype == 0){
#  temp1 <- starttimes[m] #as.numeric(names(which.max(table(Mod$sims.list$startyear[,2])))) # threshold year1
#  temp2 <- lengths[m] #temp1+as.numeric(names(which.max(table(Mod$sims.list$length[,2]))))   # threshold year2
#  if(temp2>=29) temp2 <- 60
#  model2func <- function(x){
#    new <- startabun[m,r]
#    K <- startabun[m,r]
#    if(x>0){
#      for (t in 1:x){
#        temp3 <- ifelse((t<=temp2&t>temp1),1,0)
#        lognew <- log(new) + mean(Mod$sims.list$beta[,2])*((K-new)/K)
#        K <- K + mean(Mod$sims.list$hlrate2)*temp3
#        new <- exp(lognew)
#      }
#    }
#    r <- c(new,K)
#    return(r)
#  }
#  newx <- numeric((nyears2*gentime2)+10)
#  newy <- numeric((nyears2*gentime2)+10)
#  newK <- numeric((nyears2*gentime2)+10)
#  for (x in 1:(nyears2*gentime2+10)){
#    newy[x] <- model2func(x-1)[1]
#    newK[x] <- model2func(x-1)[2]
#    newx[x] <- x-1
#  }
#  points(newx,newy,type="l",col="black",lwd=linewd[3],lty=linetype[3])
#  points(newx,newK,type="l",col=gray(.1),lwd=linewd[1],lty=linetype[1])
#}

#table(Mod$sims.list$mod)
                                                                     #############
                                                                        # MODEL 4: CONSTANT HARVEST
#temp1 <- as.numeric(names(which.max(table(Mod$sims.list$startyear[,3]))))-1 # threshold year1
#temp2 <- temp1 + as.numeric(names(which.max(table(Mod$sims.list$length[,3]))))   # threshold year2
#model4func <- function(x){
#  new <- startabun[m,r]
#  if(x>0){
#    for (t in 1:x){
#      temp3 <- ifelse((t<=temp2&t>temp1),1,0)
#      lognew <- log(new-(0/gentime2)-temp3*(mean(Mod$sims.list$constharv)/gentime2)) + mean(Mod$sims.list$beta[,3])*((mean(Mod$sims.list$K0[,3])-new)/mean(Mod$sims.list$K0[,3]))   # mean(Mod$sims.list$constharv0)
#      new <- exp(lognew)
#    }
#  }
#  return(new)
#}
#newx <- numeric((nyears2*gentime2)+6)
#newy <- numeric((nyears2*gentime2)+6)
#for (x in 1:(nyears2*gentime2+6)){
#  newy[x] <- model4func(x-1)
#  newx[x] <- x-1
#}
#points(newx,newy,type="l",col="black",lwd=linewd[4],lty=linetype[4])

                                                                     #############
                                                                       # MODEL 3: CONSTANT FRACTION HARVEST

#if(plotmodel == 3 & plottype=="pop size"){
#  temp1 <- as.numeric(names(which.max(table(Mod$sims.list$startyear[,3]))))-1 # threshold year1
#  temp2 <- temp1+as.numeric(names(which.max(table(Mod$sims.list$length[,3]))))   # threshold year2
#  if(temp2>=29) temp2 <- 60
#  model3func <- function(x){
#    new <- startabun[m,r]
#    if(x>0){
#      for (t in 1:x){
#        temp3 <- ifelse((+t<=temp2&t>temp1),1,0)
#        lognew <- log(new) + log(1-0) + temp3*log(1-as.numeric(quantile(Mod$sims.list$fracharv, 0.5))) + 
#                    mean(Mod$sims.list$beta[,3])*((mean(Mod$sims.list$K0[,2])-new)/mean(Mod$sims.list$K0[,2])) -
#                            ifelse(new<=mean(Mod$sims.list$K0[,2])*.5,mean(Mod$sims.list$alleestr)*((mean(Mod$sims.list$K0[,2])*.5-new)/mean(Mod$sims.list$K0[,2])*.5),0 )      # as.numeric(quantile(Mod$sims.list$fracharv0, 0.5),0)
#        new <- exp(lognew)
#      }
#    }
#    return(new)
#  }
#  newx <- numeric((nyears2*gentime2)+10)
#  newy <- numeric((nyears2*gentime2)+10)
#  for (x in 1:(nyears2*gentime2+10)){
#    newy[x] <- model3func(x-1)
#    newx[x] <- x-1
#  }
#  points(newx,newy,type="l",col="black",lwd=linewd[3],lty=linetype[3])
#  abline(h=mean(Mod$sims.list$K0[,2]),col=gray(.1),lwd=linewd[1],lty=linetype[1])
#}

                                                                     #############
                                                                        # MODEL 5: CONSTANT FRACTION HARVEST, exponential growth

#temp1 <- as.numeric(names(which.max(table(Mod$sims.list$startyear[,5]))))-1 # threshold year1
#temp2 <- temp1+as.numeric(names(which.max(table(Mod$sims.list$length[,5]))))   # threat length
#model6func <- function(x){
#  new <- startabun[m,r]
#  if(x>0){
#    for (t in 1:x){
#      temp3 <- ifelse((t<=temp2&t>temp1),1,0)
#      lognew <- log(new) + log(1-0)+ temp3*log(1-as.numeric(quantile(Mod$sims.list$fracharv2, 0.5))) + mean(Mod$sims.list$beta[,5])*((mean(Mod$sims.list$K0[,4])-new)/mean(Mod$sims.list$K0[,4]))
#      new <- exp(lognew)
#    }
#  }
#  return(new)
#}
#newx <- numeric((nyears2*gentime2)+6)
#newy <- numeric((nyears2*gentime2)+6)
#for (x in 1:(nyears2*gentime2+6)){
#  newy[x] <- model6func(x-1)
#  newx[x] <- x-1
#}
#points(newx,newy,type="l",col="black",lwd=linewd[6],lty=linetype[6])


##########################################
#############            PLOT RANGE OF EXPECTED TRAJECTORIES  (95% posterior predictive interval)

if(plottype=="growth rate"){
  if(plotmodel == 1){
    lines(c(2:(counter+1)),predict(loess(apply(m1,2, mean)~c(2:(counter+1)),span=.4)),lty=1, lwd=2,col=gray(.5) )
    #lines(c(2:(counter+1)),predict(loess(apply(m1,2,function(t) quantile(t,.25))~c(2:(counter+1)),span=.5)),lwd=1,col=gray(.5) )
    #lines(c(2:(counter+1)),predict(loess(apply(m1,2,function(t) quantile(t,.75))~c(2:(counter+1)),span=.5)),lwd=1,col=gray(.5) ) 
     #lines(c(1:counter),predict(loess(apply(m1,2,function(t) quantile(t,.025))~c(1:counter),span=1)),lwd=2,col=gray(.99) ) 
    #polygon(c(c(1:counter),c(counter:1)),c(predict(loess(apply(m1,2,function(t) quantile(t,.975))~c(1:counter),span=.5)), 
    #       predict(loess(apply(m1,2,function(t) quantile(t,.025))~c(1:counter),span=1))[counter:1]),density=density[1],
    #       angle=angle[1],border=NA,col=gray(.6) )
  }

  if(plotmodel == 2){
    lines(c(2:(counter+1)),predict(loess(apply(m2,2, mean)~c(2:(counter+1)),span=.4)),lty=2,lwd=2,col=gray(.5) )
    #lines(c(2:(counter+1)),predict(loess(apply(m2,2,function(t) quantile(t,.25))~c(2:(counter+1)),span=.5)),lwd=1,col=gray(.5) )
    #lines(c(2:(counter+1)),predict(loess(apply(m2,2,function(t) quantile(t,.75))~c(2:(counter+1)),span=.5)),lwd=1,col=gray(.5) ) 
     #lines(c(1:counter),predict(loess(apply(m2,2,function(t) quantile(t,.025))~c(1:counter),span=1)),lwd=2,col=gray(.99) ) 
    #polygon(c(c(1:counter),c(counter:1)),c(predict(loess(apply(m2,2,function(t) quantile(t,.975))~c(1:counter),span=.5)), 
    #         predict(loess(apply(m2,2,function(t) quantile(t,.025))~c(1:counter),span=1))[counter:1]),density=density[2],
    #         angle=angle[1],border=NA,col=gray(.6) )
  }
  
  if(plotmodel == 3){
    lines(c(2:(counter+1)),predict(loess(apply(m3,2, mean)~c(2:(counter+1)),span=.4)),lty=3,lwd=2,col=gray(.5) )
    #lines(c(2:(counter+1)),predict(loess(apply(m3,2,function(t) quantile(t,.25))~c(2:(counter+1)),span=.5)),lwd=1,col=gray(.5) )
    #lines(c(2:(counter+1)),predict(loess(apply(m3,2,function(t) quantile(t,.75))~c(2:(counter+1)),span=.5)),lwd=1,col=gray(.5) ) 
     #lines(c(1:counter),predict(loess(apply(m3,2,function(t) quantile(t,.025))~c(1:counter),span=.5)),lwd=2,col=gray(.99) ) 
    #polygon(c(c(1:counter),c(counter:1)),c(predict(loess(apply(m3,2,function(t) quantile(t,.975))~c(1:counter),span=.5)), 
    #         predict(loess(apply(m3,2,function(t) quantile(t,.025))~c(1:counter),span=1))[counter:1]),density=density[3],
    #         angle=angle[1],border=NA, col=gray(.6) )  # angle[3]
  }
}
                                           # plot the time series again...
if(plottype=="pop size")  points(c(1:counter),ts2,type="b",lwd=2,col=gray(.4))   # first plot out the original time series time (years), abundance

if(plottype=="growth rate")  points(c((gentime2+1):(nyears2*gentime2)),ts2,type="p",col=gray(.4),cex=1.5)   # first plot out the original time series time (years), abundance

if(plottype=="growth rate"&plot==3){
  #text(locator(1),paste("Most selected mechanism: \n", most_selected[m]))
  legend(locator(1),lwd=c(2,2,2),lty=c(1,2,3),
       legend= c(paste(cand_models[1]," ",postweights[1],"%",sep=""),paste(cand_models[2]," ",postweights[2],"%",sep=""),paste(cand_models[3]," ",postweights[3],"%",sep="")),bty="n",
        col=c(gray(.5),gray(.5),gray(.5)))
}

#if(plot==3) mtext("Time (years)", side=1, padj = 3)
#if(plottype=="pop size"&plot==2) mtext("Abundance", side=2, padj = -3.2)
if(plottype=="growth rate"&plot==1) mtext("Population growth rate", side=2, padj = -3.2)
#if(plottype=="pop size") mtext(paste(cand_models[plot],", posterior weight = ",round(postweights[plot],1),"%",sep=""),side=3,padj=2,cex=.8)
if(plottype=="pop size") mtext(paste(cand_models[plot],": ",round(postweights[plot],1),"%",sep=""),side=3,padj=2,cex=.8,adj=0.1)

}




dev.off()




graphics.off()






#############################
################  END SCRIPT


