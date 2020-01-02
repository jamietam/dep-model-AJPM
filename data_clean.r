# Set directory where data will be saved
mainDir <- "C:/Users/JT936/Documents/GitHub/dep-model-AJPM/"
setwd(file.path(mainDir))

# Load NSDUH 2002-2017 data (may take some time as this is a large file)
load("NSDUH_2002_2017.RData")
nsduh <- CONCATPUF_0217_031919

## Harmonize variables and reduce size of dataset
nsduh <- nsduh[c("year", "CATAG6","irsex","amdelt","amdeyr", "vestr","verep","ANALWC1","ajamdelt","ajamdeyr")] # Keep only the variables needed 
nsduh <- subset(nsduh, year>=2005 & CATAG6 >1) # only keep data after 2005 for adults ages 18+

# Due to questionnaire changes in 2008, the variables AMDELT and AMDEYR are not comparable pre-2008 vs. 2008-2017. 
table(nsduh$ajamdeyr, nsduh$year) 
table(nsduh$amdeyr, nsduh$year)
nsduh$amdeyr <- ifelse(nsduh$year<2008,nsduh$ajamdeyr,nsduh$amdeyr) # Adjusted variables AJAMDELT and AJAMDEYR were developed to allow for comparisons of adult MDE data from 2005-2017
table(nsduh$amdeyr, nsduh$year)

table(nsduh$ajamdelt, nsduh$year)
table(nsduh$amdelt, nsduh$year)
nsduh$amdelt <- ifelse(nsduh$year<2008,nsduh$ajamdelt,nsduh$amdelt) # Adjusted variables AJAMDELT and AJAMDEYR were developed to allow for comparisons of adult MDE data from 2005-2017
table(nsduh$amdelt, nsduh$year)

## Assign MD status based on past year and lifetime reports of MDE

# Current MD = past year MDE
nsduh$dep[nsduh$amdeyr==1] <- 1
nsduh$dep[nsduh$amdeyr==2] <- 0

# Former MD = no past year MDE, but lifetime MDE
nsduh$notdep[nsduh$amdeyr==2 & nsduh$amdelt==1] <- 1
nsduh$notdep[nsduh$amdeyr==2 & nsduh$amdelt==2] <- 0
nsduh$notdep[nsduh$amdeyr==1] <- 0

# Never MD = no reported lifetime MDE (but may be subject to recall error)
nsduh$nevdep[nsduh$amdeyr==2 & nsduh$amdelt==2] <- 1 
nsduh$nevdep[nsduh$amdeyr==1& nsduh$amdelt==1] <- 0
nsduh$nevdep[nsduh$amdeyr==2 & nsduh$amdelt==1] <- 0

# Ever MD = reported lifetime MDE
nsduh$everdep[nsduh$nevdep==1] <- 0 
nsduh$everdep[nsduh$nevdep==0] <- 1

# Load survey packages
library(survey)
options(survey.lonely.psu="adjust")
nsduh <- nsduh[order(nsduh$vestr, nsduh$verep),] # re-order survey design variables

## Generate age group-specific prevalences and confidence intervals across adult population
getprevsbyage <- function(groupvar, subpop){
  byagegroup = NULL
  alladults = NULL
  for (y in 2005:2017){
    svy <-svydesign(id=~verep, strata=~vestr, nest=TRUE, weights=~ANALWC1, data=subset(subpop,year==y))

    prev <-svymean(as.formula(paste("~",groupvar)),design=svy,na.rm=TRUE) 
    alladults  <- rbind(alladults, data.frame(y,"total",groupvar, deparse(substitute(subpop)), prev[1], confint(prev)[1,1], confint(prev)[1,2]))          
    
    agegroupnames <- c("18to25", "26to34","35to49", "50to64","65plus")
    for (k in 2:6){
      prev <-svymean(as.formula(paste("~",groupvar)),design=subset(svy,CATAG6==k),na.rm=TRUE) # 
      byagegroup <- rbind(byagegroup, data.frame(y,agegroupnames[k-1],groupvar, deparse(substitute(subpop)), prev[1],confint(prev)[1,1], confint(prev)[1,2]))          
    }
  }
  names(alladults) <- names(byagegroup)
  byagegroup <- rbind(alladults, byagegroup)
  colnames(byagegroup) <- c("survey_year","age","status","subpopulation","prev","prev_lowCI","prev_highCI")
  return(byagegroup)
}

## Create new dataframe with NSDUH prevalences for model fitting
depprevs_by_year <- NULL
gender <- c("males", "females")
for (x in 1:2){# irsex: 1 = males, 2 = females
  totalpop <- subset(nsduh, nsduh$irsex==x)
  
  depprevs_by_year <- rbind(depprevs_by_year, cbind(gender[x], getprevsbyage("dep",totalpop)))
  depprevs_by_year <- rbind(depprevs_by_year, cbind(gender[x], getprevsbyage("notdep",totalpop)))
  depprevs_by_year <- rbind(depprevs_by_year, cbind(gender[x], getprevsbyage("nevdep",totalpop)))
  depprevs_by_year <- rbind(depprevs_by_year, cbind(gender[x], getprevsbyage("everdep",totalpop)))
} 

colnames(depprevs_by_year)[1] <- "gender"
depprevs_by_year$prev=as.numeric(depprevs_by_year$prev)
depprevs_by_year$prev_highCI=as.numeric(depprevs_by_year$prev_highCI)
depprevs_by_year$prev_lowCI=as.numeric(depprevs_by_year$prev_lowCI)
depprevs_by_year$group <- paste(depprevs_by_year$gender, depprevs_by_year$status, depprevs_by_year$age, sep="_")

save(depprevs_by_year, file="depprevs_2005-2017.rda")