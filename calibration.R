library(Bhat)
load("depprevs_2005-2017.rda") # load NSDUH data for calibration

# If performing calibration of underreporting and incidence parameters, use these lines of code in 'main.R' before running "source('~/GitHub/dep-model-AJPM/main.r')"
# allparamsF = read.xlsx(paste0(mainDir, "parameters.xlsx"),sheet=paste0("dep_females_splines"),rowNames=TRUE,colNames=TRUE) 
# allparamsM = read.xlsx(paste0(mainDir, "parameters.xlsx"),sheet=paste0("dep_males_splines"),rowNames=TRUE,colNames=TRUE) 

# If calibrationg increased incidence probabilities starting in year_SF (for 18-25 year olds), use these lines of code in 'main.R' before running "source('~/GitHub/dep-model-AJPM/main.r')"
# allparamsF = read.xlsx(paste0(mainDir, "parameters.xlsx"),sheet=paste0("dep_females_incSF"),rowNames=TRUE,colNames=TRUE) 
# allparamsM = read.xlsx(paste0(mainDir, "parameters.xlsx"),sheet=paste0("dep_males_incSF"),rowNames=TRUE,colNames=TRUE) 

source('~/GitHub/dep-model-AJPM/main.r')

xF <- list(label=paramsnamesF, est=paramsF,low=lowervectorF,upp=uppervectorF) # est = parameter starting values
xM <- list(label=paramsnamesM, est=paramsM,low=lowervectorM,upp=uppervectorM) # est = parameter starting values

# Get NSDUH prevs by age group 
getnsduhprevs <- function(depsmkprevs_by_year,assignedsex,numpop,denompop){
  nsduhdata = NULL
  nsduhdata <- melt(depsmkprevs_by_year,id.vars=c("group","survey_year","age","gender","subpopulation","status"),measure.vars = c("prev") )
  nsduhdata <- cast(nsduhdata, group~survey_year,mean,subset=gender==assignedsex&subpopulation==denompop&status==numpop)
  nsduhdata <- nsduhdata[c(-1)] # remove group name column
  row.names(nsduhdata)<-agerownames
  return(nsduhdata)
}

# Get model sum of squares all ages 
getsumdiffs = function(out,whichgender){
  years<- c("X2005","X2006","X2007","X2008","X2009","X2010","X2011","X2012","X2013","X2014","X2015","X2016","X2017") # only look at output for years where NSDUH data are available
  years2 <- c("2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017")
  depdiffs <-
    rowSums(subset(as.data.frame(out[1]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,whichgender,"nevdep","totalpop"),select=years2))^2 +
    rowSums(subset(as.data.frame(out[2]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,whichgender,"dep","totalpop"),select=years2))^2  +
    rowSums(subset(as.data.frame(out[3]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,whichgender,"notdep","totalpop"),select=years2))^2 +
    rowSums(subset(as.data.frame(out[4]),select=years) - subset(getnsduhprevs(depsmkprevs_by_year,whichgender,"everdep","totalpop"),select=years2))^2
  return(sum(depdiffs))
}

# Get model sum of squares for 18-25 year olds 2016-2017 
getsumdiffs18to25 = function(out,whichgender){
  years<- c("X2016","X2017") # only fit output for years 2016 and 2017 for ages 18-25 group
  years2 <- c("2016","2017") 
  depdiffs <- 
    rowSums(subset(as.data.frame(out[1]),select=years)["18to25",] - subset(getnsduhprevs(depsmkprevs_by_year,whichgender,"nevdep","totalpop"),select=years2)["18to25",])^2 +
    rowSums(subset(as.data.frame(out[2]),select=years)["18to25",] - subset(getnsduhprevs(depsmkprevs_by_year,whichgender,"dep","totalpop"),select=years2)["18to25",])^2  +
    rowSums(subset(as.data.frame(out[3]),select=years)["18to25",] - subset(getnsduhprevs(depsmkprevs_by_year,whichgender,"notdep","totalpop"),select=years2)["18to25",])^2 +
    rowSums(subset(as.data.frame(out[4]),select=years)["18to25",] - subset(getnsduhprevs(depsmkprevs_by_year,whichgender,"everdep","totalpop"),select=years2)["18to25",])^2
  return(sum(depdiffs))
}

# Perform bhat estimation for females
ML_bhatF=function(paramsF){
  out = main(getmodelprevs,"females",allparamsF,paramsF,paramsnamesF, 2016)
  LL = sum(getsumdiffs(out,"females"))   # Least squares for all ages and all years
  # LL = sum(getsumdiffs18to25(out,"females"))   # Least squares for 18-25 year olds in 2016 & 2017
  cat(LL,paramsF,'\n')
  return(LL)
}
resbhatF=dfp(xF,ML_bhatF) 

# Perform bhat estimation for males
ML_bhatM=function(paramsM){
  out = main(getmodelprevs,"males",allparamsM,paramsM,paramsnamesM, 2016)
  LL = sum(getsumdiffs(out,"males"))   # Least squares for all ages and all years
  # LL = sum(getsumdiffs18to25(out,"males"))   # Least squares for 18-25 year olds in 2016 & 2017
  cat(LL,paramsM,'\n')
  return(LL)
}
resbhatM=dfp(xM,ML_bhatM) 
