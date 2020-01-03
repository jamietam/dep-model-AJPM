# Set directory where results will be saved
mainDir <- "~/GitHub/dep-model-AJPM/"
setwd(file.path(mainDir))
library(openxlsx)
library(splines)
load("depprevs_2005-2017.rda") # load NSDUH data for calibration

# If running the fully calibrated model (no need for additional estimation of parameters), use this line of code:
allparamsF = read.xlsx(paste0(mainDir, "parameters.xlsx"),sheet=paste0("dep_females_nocalib"),rowNames=TRUE,colNames=TRUE) # Adjust parameters in this excel file with estimates to be used in the model
allparamsM = read.xlsx(paste0(mainDir, "parameters.xlsx"),sheet=paste0("dep_males_splines"),rowNames=TRUE,colNames=TRUE) 

# If performing calibration of underreporting and incidence parameters use this line of code:
# allparamsF = read.xlsx(paste0(mainDir, "parameters.xlsx"),sheet=paste0("dep_females_splines"),rowNames=TRUE,colNames=TRUE) 
# allparamsM = read.xlsx(paste0(mainDir, "parameters.xlsx"),sheet=paste0("dep_males_splines"),rowNames=TRUE,colNames=TRUE) 

# If calibrationg increased incidence probabilities starting in year_SF (for 18-25 year olds), use this line of code:
# allparamsF = read.xlsx(paste0(mainDir, "parameters.xlsx"),sheet=paste0("dep_females_incSF"),rowNames=TRUE,colNames=TRUE) 
# allparamsM = read.xlsx(paste0(mainDir, "parameters.xlsx"),sheet=paste0("dep_males_incSF"),rowNames=TRUE,colNames=TRUE) 

## FEMALES
paramsF = as.vector(subset(allparamsF,bhat==1)[['estimate']]) # Parameters where bhat = 1 can be estimated by bhat package in 'calibration.R'
paramsnamesF = rownames(subset(allparamsF,bhat==1))  # names of parameters
lowervectorF = as.vector(subset(allparamsF,bhat==1)[['lower']]) # lower bounds for bhat parameters
uppervectorF = as.vector(subset(allparamsF,bhat==1)[['upper']]) # upper bounds for bhat parameters

## MALES
paramsM = as.vector(subset(allparamsM,bhat==1)[['estimate']]) 
paramsnamesM = rownames(subset(allparamsM,bhat==1))
lowervectorM = as.vector(subset(allparamsM,bhat==1)[['lower']]) 
uppervectorM = as.vector(subset(allparamsM,bhat==1)[['upper']]) 

startyear = 1900 # burn-in period starting point
endyear = 2017 # 2017 
startage = 0
endage = 99
Ny= endyear - startyear + 1 
Na= endage - startage + 1
emptycompartment <- matrix(0, nrow = Na, ncol = Ny, dimnames=list(c(startage:endage),c(startyear:endyear))) # create matrix of zeroes for compartments
agerownames<-c("18to25", "26to34", "35to49", "50to64",  "65plus", "total")
agegroupstart <- c(18,26,35,50,65,18)
agegroupend <- c(25,34,49,64,99,99)

# Get model prevs by age group --------------------------------------------
getmodelprevs <- function(numerator,denominator){
  numerator = as.data.frame(numerator)
  denominator = as.data.frame(denominator)
  prevs = NULL
  for (a in c(1:length(agegroupstart))) {
    prevs <- rbind(prevs,colSums(numerator[c((agegroupstart[a]+1):(agegroupend[a]+1)), ],na.rm=TRUE)/
                     colSums(denominator[c((agegroupstart[a]+1):(agegroupend[a]+1)), ],na.rm=TRUE))
  }
  row.names(prevs)<-agerownames 
  return(prevs)
}

# Main model --------------------------------------------------------------
main <- function(getmodelprevs, whichgender, allparamsF, paramsF,paramsnamesF, year_SF){
  
  # model parameters are retrieved from the parameters.xlsx file, if they are bhat=1, then they are read in as a params vector for bhat estimation
  RRdepr_death =  c(rep(1,17),rep(ifelse(allparamsF["RRdepr_death","bhat"]==0,allparamsF["RRdepr_death","estimate"],paramsF[match("RRdepr_death",paramsnamesF)]),Na-18))
  death_rate = read.xlsx(paste0(mainDir, "hmd_death_rates.xlsx"),sheet=paste0("hmd_",whichgender),rowNames=TRUE, colNames=TRUE, check.names=FALSE) 
  pop = read.xlsx(paste0(mainDir, "/usproj2000-2050.xlsx"),sheet=whichgender,rowNames=TRUE, colNames=TRUE, check.names=FALSE)
  dep1inc = read.xlsx(paste0(mainDir, "incidence_eaton.xlsx"),sheet=whichgender,rowNames=TRUE, colNames=FALSE, check.names=FALSE)
  
  # scale all recovery rates by the scaling factor deprecov_SF
  dep1recov_rate =  c(rep(ifelse(allparamsF["dep1recov_rate","bhat"]==0,allparamsF["dep1recov_rate","estimate"],paramsF[match("dep1recov_rate",paramsnamesF)]),Na-1))
  deprecov_SF = matrix(1,99,1)
  deprecov_rate = dep1recov_rate*deprecov_SF
  
  # recovery rates for depression are the same regardless of whether 1st or subsequent depressive episode
  
  deprinc = rep(ifelse(allparamsF["depr_inc","bhat"]==0,allparamsF["depr_inc","estimate"],paramsF[match("depr_inc",paramsnamesF)]),Na-1)
  deprinc_SF = matrix(1,99,1)
  depr_inc = deprinc * deprinc_SF
  
  depinc1 = ifelse(allparamsF["depinc1","bhat"]==0,allparamsF["depinc1","estimate"],paramsF[match("depinc1",paramsnamesF)])
  depinc2 = ifelse(allparamsF["depinc2","bhat"]==0,allparamsF["depinc2","estimate"],paramsF[match("depinc2",paramsnamesF)])
  depinc3 = ifelse(allparamsF["depinc3","bhat"]==0,allparamsF["depinc3","estimate"],paramsF[match("depinc3",paramsnamesF)])
  
  # scale all incidence rates by scaling factor starting in year_SF
  inc_SF = ifelse(allparamsF["inc_SF","bhat"]==0,allparamsF["inc_SF","estimate"],paramsF[match("inc_SF",paramsnamesF)])
  
  if (whichgender=="females"){
    MP=ns(0:21,knots=c(13,18)) ### Matrix of X's
    Rps=predict(MP,21)[1,] ## Predicts y value given a set of X's for age 22  # 0.0051285304 anchor at age 22
    y=c()
    for (j in 0:21){
      y=c(y,0.0051285304*exp(sum((MP[j,]-Rps)*c(depinc1,depinc2,depinc3)))) ## multiply by coefficients and sum , exp makes it positive for incidence
    }
    dep1_inc=dep1inc
    dep1_inc[0:22,]<-y
    dep1_inc[0:12,]<-rep(0,12) # assumes no 1st MDE before age 12
    scaleddep1_inc <- dep1_inc*inc_SF # multiply all incidence probabilities by scaling factor
  }
  
  if (whichgender=="males"){
    MP=ns(0:28,knots=c(13,18))
    Rps=predict(MP,28)[1,] ## Predicts y value given a set of X's for age 22  # 0.0019072600 anchor at age 29
    y=c()
    for (j in 0:28){
      y=c(y,0.0019072600*exp(sum((MP[j,]-Rps)*c(depinc1,depinc2,depinc3)))) ## multiply by coefficients and sum , exp makes it positive for incidence
    }
    dep1_inc=dep1inc
    dep1_inc[0:29,]<-y
    dep1_inc[0:12,]<-rep(0,12) # assumes no 1st MDE before age 12
    scaleddep1_inc <- dep1_inc*inc_SF # multiply all incidence probabilities by scaling factor
  }
  
  # Age-group categorical underreporting probabilities 
  forget1 = ifelse(allparamsF["forget1","bhat"]==0,allparamsF["forget1","estimate"],paramsF[match("forget1",paramsnamesF)])
  forget2 = ifelse(allparamsF["forget2","bhat"]==0,allparamsF["forget2","estimate"],paramsF[match("forget2",paramsnamesF)])
  forget3 = ifelse(allparamsF["forget3","bhat"]==0,allparamsF["forget3","estimate"],paramsF[match("forget3",paramsnamesF)])
  forget4 = ifelse(allparamsF["forget4","bhat"]==0,allparamsF["forget4","estimate"],paramsF[match("forget4",paramsnamesF)])
  forget5 = ifelse(allparamsF["forget5","bhat"]==0,allparamsF["forget5","estimate"],paramsF[match("forget5",paramsnamesF)])
  forget = matrix(NA,99,1)
  forget[1:17] = rep(0,17)
  forget[18:25] = forget1
  forget[26:34] = forget2
  forget[35:49] = forget3
  forget[50:64] = forget4
  forget[65:99] = forget5
  
  # Compartments / state variables ------------------------------------------
  
  # Initialize population - Model compartments are organized with age 0-99 as rows, and year as columns
  matrix.names<-c('nevdeptrue', 'dep1','hdepr','depr','forgot')
  for (name in matrix.names) assign(name,emptycompartment)
  
  nevdeptrue[paste(startage),1:Ny] <- as.matrix(pop[paste(startage),1:Ny])   # Takes empty compartment and populates the top row of the matrix with the number of 0-yrolds
  
  for (y in c((startyear+1):(year_SF-1))){
    py = paste(y - 1)
    nevdeptrue[2:Na,paste(y)] <- nevdeptrue[1:Na-1,py]*(1-dep1_inc[(startage+1):(endage),])*(1-death_rate[(startage+1):endage,py])
    dep1[2:Na,paste(y)] <- dep1[1:Na-1,py]*(1-deprecov_rate)*(1-RRdepr_death*death_rate[(startage+1):endage,py]) + (dep1_inc[(startage+1):(endage),])*nevdeptrue[1:Na-1,py]
    hdepr[2:Na,paste(y)] <- hdepr[1:Na-1,py]*(1-depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_rate[(startage+1):endage,py])*(1-forget[(startage+1):(endage)])  + (deprecov_rate)*dep1[1:Na-1,py] + (deprecov_rate)*depr[1:Na-1,py]
    depr[2:Na,paste(y)] <- depr[1:Na-1,py]*(1-deprecov_rate)*(1-RRdepr_death*death_rate[(startage+1):endage,py])  + (depr_inc[(startage+1):(endage),])*hdepr[1:Na-1,py] + (depr_inc[(startage+1):(endage),])*forgot[1:Na-1,py]
    forgot[2:Na,paste(y)] <- forgot[1:Na-1,py]*(1-depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_rate[(startage+1):endage,py])+(forget[(startage+1):(endage)])*hdepr[1:Na-1,py]
  }
  for (y2 in c(year_SF:endyear)){
    py2 = paste(y2 - 1)
    nevdeptrue[2:Na,paste(y2)] <- nevdeptrue[1:Na-1,py2]*(1-scaleddep1_inc[(startage+1):(endage),])*(1-death_rate[(startage+1):endage,py2])
    dep1[2:Na,paste(y2)] <- dep1[1:Na-1,py2]*(1-deprecov_rate)*(1-RRdepr_death*death_rate[(startage+1):endage,py2]) + (scaleddep1_inc[(startage+1):(endage),])*nevdeptrue[1:Na-1,py2]
    hdepr[2:Na,paste(y2)] <- hdepr[1:Na-1,py2]*(1-depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_rate[(startage+1):endage,py2])*(1-forget[(startage+1):(endage)])  + (deprecov_rate)*dep1[1:Na-1,py2] + (deprecov_rate)*depr[1:Na-1,py2]
    depr[2:Na,paste(y2)] <- depr[1:Na-1,py2]*(1-deprecov_rate)*(1-RRdepr_death*death_rate[(startage+1):endage,py2])  + (depr_inc[(startage+1):(endage),])*hdepr[1:Na-1,py2] + (depr_inc[(startage+1):(endage),])*forgot[1:Na-1,py2]
    forgot[2:Na,paste(y2)] <- forgot[1:Na-1,py2]*(1-depr_inc[(startage+1):(endage),])*(1-RRdepr_death*death_rate[(startage+1):endage,py2])+(forget[(startage+1):(endage)])*hdepr[1:Na-1,py2]
  }
  
  # Get population counts for each subpopulation ----------------------------
  nevdep = nevdeptrue + forgot
  totalpop = nevdep+dep1 + hdepr + depr
  deppop = dep1+depr
  everdeppop = dep1 + hdepr + depr # does not account for recall error
  everdeppoptrue = dep1 + hdepr + depr + forgot # adjusts for recall error
  
  d1 <- getmodelprevs(nevdep,totalpop) # never MDE prevalence
  d2 <- getmodelprevs(deppop,totalpop) # current MDE prevalence
  d3 <- getmodelprevs(hdepr,totalpop) # former MDE prevalence 
  e1 <- getmodelprevs(everdeppop,totalpop) # lifetime MDE prevalence (no recall error adjustment)
  
  f1 <- getmodelprevs(forgot,totalpop) # prevalence of people with past MDE who 'forget' and transition to nevdep
  f2 <- getmodelprevs(forgot,nevdep) 
  
  n1 <- getmodelprevs(nevdeptrue,totalpop)
  n2 <- getmodelprevs(nevdeptrue,nevdep)
  
  d4 <- getmodelprevs(dep1, totalpop)
  d5 <- getmodelprevs(depr,totalpop)
  
  incidence <- data.frame(cbind(dep1inc[-1,],dep1_inc[-1,],scaleddep1_inc[-1,],deprinc,deprinc_SF,data.frame(depr_inc),c(1:99)))
  colnames(incidence) <- c("dep1inc","splines","scaleddep1inc_yr", "deprinc","SF_depr", "scaledrates_depr", "age")
  if (whichgender=="females"){
    incidence$dep1inc[0:21]<- NA
    incidence$dep1inc[68:99]<- NA
  }
  if (whichgender=="males"){
    incidence$dep1inc[0:28]<- NA
    incidence$dep1inc[59:99]<- NA
  }
  recovery <- data.frame(cbind(dep1recov_rate, deprecov_SF,deprecov_rate,c(1:99)))
  colnames(recovery) <- c("dep1recov","SF","scaledrates","age")
  forget <- data.frame(cbind(forget,c(1:99)))
  colnames(forget) <- c("forget_prob","age")
  modeldepprevdata <- list(d1,d2,d3,e1,d4,d5,incidence,recovery,forget,f1,f2,n1,n2, everdeppoptrue, totalpop)
  
  return(modeldepprevdata)
}

# Run the model
outF = main(getmodelprevs, "females", allparamsF, paramsF,paramsnamesF, 2016) # runs model using parameters specified in excel sheet OR using bhat estimates
outM = main(getmodelprevs, "males", allparamsM, paramsM,paramsnamesM, 2016)
