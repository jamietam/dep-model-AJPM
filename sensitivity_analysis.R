# Sensitivity Analyses using the 'pse' package
# Parameter Space Exploration of deterministic models - Tutorial here: https://cran.r-project.org/web/packages/pse/vignettes/pse_tutorial.pdf
library(pse)
library(Bhat)
mainDir <- "~/GitHub/dep-model-AJPM/"
setwd(file.path(mainDir))
source('~/GitHub/dep-model-AJPM/main.r')

# initial values
recovery= 0.173
recurrence = 0.058
RRmd = 1.71

# Function that runs the model for a single combination of parameter values
singlerun <- function(recovery,recurrence,RRmd){
  X <- array()
  
  allparamsF["depr_inc","estimate"] = recurrence
  allparamsF["dep1recov_rate","estimate"] = recovery
  allparamsF["RRdepr_death","estimate"] = RRmd 
  allparamsM["depr_inc","estimate"] = recurrence
  allparamsM["dep1recov_rate","estimate"] = recovery
  allparamsM["RRdepr_death","estimate"] = RRmd
  
  xF <- list(label=paramsnamesF, est=paramsF,low=lowervectorF,upp=uppervectorF) # est = parameter starting values
  resbhatF=dfp(xF,ML_bhatF)
  paramsFbhat=resbhatF$est
  xM <- list(label=paramsnamesM, est=paramsM,low=lowervectorM,upp=uppervectorM) # est = parameter starting values
  resbhatM=dfp(xM,ML_bhatM)
  paramsMbhat=resbhatM$est
  
  outF = main(getmodelprevs, "females", allparamsF, paramsFbhat,paramsnamesF, 2016) # runs model using parameters specified in excel sheet OR using bhat estimates
  d1F <- data.frame(outF[1],check.names = FALSE) # never depressed prevalence among totalpop
  f1F <- data.frame(outF[10],check.names=FALSE) # forgot prevalence among total pop
  d0F <- d1F-f1F # nevdeptrue for totalpop
  e0F <- 1-d0F # everdeptrue for totalpop
  
  outM = main(getmodelprevs, "males", allparamsM, paramsMbhat,paramsnamesM, 2016)
  d1M <- data.frame(outM[1],check.names = FALSE) # never depressed prevalence among totalpop
  f1M <- data.frame(outM[10],check.names=FALSE) # forgot prevalence among total pop
  d0M <- d1M-f1M # nevdeptrue for totalpop
  e0M <- 1-d0M # everdeptrue for totalpop
  
  # Figure S3 use "2017"
  # Figure S5 use "2015"
  X[1]<- e0F["total","2017"] # ever depressed prevalence - females
  X[2]<- e0M["total","2017"] # ever depressed prevalence - males
  X[3]<- f1F["total","2017"] # forgot prevalence - females
  X[4]<- f1M["total","2017"] # forgot prevalence - males
  
  return(X)
}


# Use LHS function to generate a hypercube for your model

modelrun <- function (my.data) { # Function that represents the model and runs the model for multiple parameter combinations
  return(mapply(singlerun, my.data[,1], my.data[,2], my.data[,3]))
}

factors <- c("recovery", "recurrence", "RRmd")  # factors: an array with the parameter names

N = 200 # number of parameter combinations to be generated

q <- c("qunif","qunif","qnorm") # names of the Probability density functions to generate the parameter values

# To get the SD, 1.71 + (1.96)(0.09693878) = 1.90 and 1.71 - (1.96)(0.08673469) = 1.54. (0.09693878+0.08673469)/2 = 0.09183674
# List containing the lists with all the parameters to the density functions
q.arg <- list( list(min=0.0865, max=0.346),list(min=0.029,max=0.116), list(mean=1.71,sd=0.09183674)) # a  list  with  the  arguments  of  each  pdf

res.names <- c("Lifetime MDE - Females","Lifetime MDE - Males","Forgot MDE - Females","Forgot MDE - Males")

myLHS <- LHS(modelrun, factors, N, q, q.arg, res.names, nboot=50)
# save(myLHS,file=paste0("LHS_200N_50B.rda"))

# Generate uncertainty analysis figures for appendix
library(ggplot2)
library(grid)
library(gridExtra)
library(plyr)
mdelhs <- as.data.frame(cbind(get.results(myLHS)[,1],"Females","Lifetime MDE"))
mdelhs <- rbind(mdelhs,cbind(get.results(myLHS)[,2],"Males","Lifetime MDE"))
colnames(mdelhs) <-c("prev","gender","status")
mdelhs$prev <- as.numeric(as.character(mdelhs$prev))
avg <- ddply(mdelhs, "gender", summarise, prev.mean=mean(prev))

mdelhs_plot<-ggplot(data=mdelhs, aes(x=prev*100, fill=gender)) + geom_density(alpha=0.5) +
  labs(title="A) Lifetime MD")+
  geom_vline(data=avg, aes(xintercept=prev.mean*100,  colour=gender), linetype="dashed", size=1)+
  scale_x_continuous(name="Prevalence (%)",limits=c(min(mdelhs$prev*100-1),max(mdelhs$prev*100+1)),breaks=seq(0,max(mdelhs$prev*100+1),1)) +
  theme(legend.title = element_blank(),plot.title = element_text(hjust = -0.14))

forgotlhs <- as.data.frame(cbind(get.results(myLHS)[,3],"Females","Forgot MDE"))
forgotlhs <- rbind(forgotlhs,cbind(get.results(myLHS)[,4],"Males","Forgot MDE"))
colnames(forgotlhs) <-c("prev","gender","status")
forgotlhs$prev <- as.numeric(as.character(forgotlhs$prev))
avg2 <- ddply(forgotlhs, "gender", summarise, prev.mean=mean(prev))

forgotlhs_plot<-ggplot(data=forgotlhs, aes(x=prev*100, fill=gender)) + geom_density(alpha=0.5) +
  labs(title="B) Recall error")+
  geom_vline(data=avg2, aes(xintercept=prev.mean*100,  colour=gender), linetype="dashed", size=1)+
  scale_x_continuous(name="Prevalence (%)",limits=c(0,max(forgotlhs$prev*100+1)),breaks=seq(0,max(forgotlhs$prev*100+1),1)) +
  theme(legend.title = element_blank(),plot.title = element_text(hjust = -0.14))

grid_arrange_shared_legend <- function(plots,columns,titletext) {
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(arrangeGrob(grobs= lapply(plots, function(x)
    x + theme(legend.position="none", plot.title = element_text(size = rel(1.0)))),ncol=columns),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight),
    top=textGrob(titletext,just="top", vjust=1,check.overlap=TRUE,gp=gpar(fontsize=9, fontface="bold"))
  )
}

# Figure S3 or Figure S5
jpeg(filename = paste0("Figure_S3_uncertainty_dist.jpg"),width=10, height=6, units ="in", res=1000)
grid_arrange_shared_legend(list(mdelhs_plot,forgotlhs_plot),2,"Uncertainty distributions")
dev.off()

# Numbers for manuscript text 
round(quantile(subset(forgotlhs,gender=="Males")$prev,probs=c(.025,.975))*100,1)
round(quantile(subset(forgotlhs,gender=="Females")$prev,probs=c(.025,.975))*100,1)
round(quantile(subset(mdelhs,gender=="Males")$prev,probs=c(.025,.975))*100,1)
round(quantile(subset(mdelhs,gender=="Females")$prev,probs=c(.025,.975))*100,1)

