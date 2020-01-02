library(ggplot2)
library(reshape)
library(grid)
library(gridBase)
library(gridExtra)
library(plyr)
library(Hmisc)
options(scipen=999)
mainDir <- "C:/Users/JT936/Documents/GitHub/dep-model-AJPM/"
setwd(file.path(mainDir)) # Set directory where results will be saved
source('~/GitHub/dep-model-AJPM/main.r')

# Females results
d1F <- data.frame(outF[1],check.names = FALSE) # never depressed prevalence among totalpop
d2F <- data.frame(outF[2],check.names = FALSE) # depressed prevalence among totalpop
d3F <- data.frame(outF[3],check.names = FALSE) # notdep prev among totalpop
e1F <- data.frame(outF[4],check.names = FALSE) # everdeppop,totalpop
d4F <- data.frame(outF[5],check.names = FALSE) # dep1, totalpop
d5F <- data.frame(outF[6],check.names = FALSE) # depr, totalpop
incidenceF = data.frame(outF[7])
recoveryF = data.frame(outF[8])
forgetF = data.frame(outF[9])
f1F <- data.frame(outF[10],check.names=FALSE) #totalpop
d0F <- d1F-f1F #nevdeptrue for totalpop
e0F <- 1-d0F #everdeptrue for totalpop
f2F <- data.frame(outF[11],check.names=FALSE) #nevpop
n1F <- data.frame(outF[12],check.names=FALSE) #totalpop
n2F <- data.frame(outF[13],check.names=FALSE) #nevpop
everdeppoptrueF <- data.frame(outF[14],check.names=FALSE)
totalpopF <- data.frame(outF[15],check.names=FALSE)

# Males results
d1M <- data.frame(outM[1],check.names = FALSE) # never depressed prevalence among totalpop
d2M <- data.frame(outM[2],check.names = FALSE) # depressed prevalence among totalpop
d3M <- data.frame(outM[3],check.names = FALSE) # notdep prev among totalpop
e1M <- data.frame(outM[4],check.names = FALSE) # everdeppop,totalpop
d4M <- data.frame(outM[5],check.names = FALSE) # dep1, totalpop
d5M <- data.frame(outM[6],check.names = FALSE) # depr, totalpop
incidenceM = data.frame(outM[7])
recoveryM = data.frame(outM[8])
forgetM = data.frame(outM[9])
f1M <- data.frame(outM[10],check.names=FALSE) #totalpop
d0M <- d1M-f1M #nevdeptrue for totalpop
e0M <- 1-d0M #everdeptrue for totalpop
f2M <- data.frame(outM[11],check.names=FALSE) #nevpop
n1M <- data.frame(outM[12],check.names=FALSE) #totalpop
n2M <- data.frame(outM[13],check.names=FALSE) #nevpop
everdeppoptrueM <- data.frame(outM[14],check.names=FALSE)
totalpopM <- data.frame(outM[15],check.names=FALSE)

# Combine plots with shared legend
grid_arrange_shared_legend <- function(plots,columns,titletext) {
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(arrangeGrob(grobs= lapply(plots, function(x)
    x + theme(legend.position="none", plot.title = element_text(size = rel(1.2)))),ncol=columns),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight),
    top=textGrob(titletext,just="top", vjust=1,check.overlap=TRUE,gp=gpar(fontsize=9, fontface="bold"))
  )
}

# Get NSDUH prevs and confidence intervals by age group 
getnsduhprevs <- function(depsmkprevs_by_year,assignedsex,numpop,denompop){
  nsduhdata = NULL
  nsduhdata <- melt(depsmkprevs_by_year,id.vars=c("group","survey_year","age","gender","subpopulation","status"),measure.vars = c("prev") )
  nsduhdata <- cast(nsduhdata, group~survey_year,mean,subset=gender==assignedsex&subpopulation==denompop&status==numpop)
  nsduhdata <- nsduhdata[c(-1)] # remove group name column
  row.names(nsduhdata)<-agerownames
  return(nsduhdata)
}
getnsduhprevsCI <- function(depsmkprevs_by_year,assignedsex,numpop,denompop){
  nsduhdatalow = NULL
  nsduhdatalow <- melt(depsmkprevs_by_year,id.vars=c("group","survey_year","age","gender","subpopulation","status"),measure.vars = c("prev_lowCI") )
  nsduhdatalow <- cast(nsduhdatalow, group~survey_year,mean,subset=gender==assignedsex&subpopulation==denompop&status==numpop)
  nsduhdatalow <- nsduhdatalow[c(-1)] # remove group name column
  row.names(nsduhdatalow)<-agerownames
  
  nsduhdatahigh = NULL
  nsduhdatahigh <- melt(depsmkprevs_by_year,id.vars=c("group","survey_year","age","gender","subpopulation","status"),measure.vars = c("prev_highCI") )
  nsduhdatahigh <- cast(nsduhdatahigh, group~survey_year,mean,subset=gender==assignedsex&subpopulation==denompop&status==numpop)
  nsduhdatahigh <- nsduhdatahigh[c(-1)] # remove group name column
  row.names(nsduhdatahigh)<-agerownames
  return(list(nsduhdatalow, nsduhdatahigh))
}

# Figure 2
incidenceplot <- ggplot(incidenceF, aes(age,dep1inc)) + 
  geom_point(aes(x=age, y=incidenceF$dep1inc, color="Eaton, 1997 - Females")) +
  geom_line(aes(x=age, y= incidenceF$splines, color="Calibrated estimates - Females" ))+
  
  geom_point(aes(x=age, y=incidenceM$dep1inc, color="Eaton, 1997 - Males"),shape=21) +
  geom_line(aes(x=age, y= incidenceM$splines, color="Calibrated estimates - Males" ),linetype="dashed")+
  
  labs(title=paste0(""))+
  xlab("Age") +
  scale_x_continuous(limits=c(0,100),breaks=seq(0,100,10)) +
  scale_y_continuous(name="Annual probability",limits=c(0,0.07000),breaks=seq(0,0.07,0.005)) +
  theme(legend.position = c(1, 1),legend.justification = c(1, 1), legend.background = element_rect(fill = "transparent"))+
  scale_colour_manual(name="",values = c("black", "black","black","black"), guide = guide_legend(override.aes = list(
    linetype = c("solid", "dashed","blank","blank"), shape = c(NA,NA,16,21))))

jpeg(filename = paste0(mainDir, "Figure_2_incidence.jpg"),width=5, height=6, units ="in", res=1000)
incidenceplot
dev.off()

# Figure 3
underreportingplot <- ggplot(forgetF, aes(age,forget_prob)) +
  geom_line(aes(x=age,y=forgetF$forget_prob,color="Females")) +
  geom_line(aes(x=age,y=forgetM$forget_prob,color="Males"),linetype="dashed")+
  labs(title=paste0(""))+
  xlab("Age")+
  theme(legend.position = c(1, 0),legend.justification = c(1, 0), legend.background = element_rect(fill = "transparent"))+
  scale_x_continuous(limits=c(0,100),breaks=seq(0,100,10)) +
  scale_y_continuous(name="Annual probability",limits=c(0,1),breaks=seq(0,1,0.10))+
  scale_color_manual(name="",values=c("black","black"),
                     guide=guide_legend(override.aes=list(linetype=c("solid","dashed"))))

jpeg(filename = paste0(mainDir, "Figure_3_underreport.jpg"),width=5, height=6, units ="in", res=1000)
underreportingplot
dev.off()

# Figure 4
lifetimeMDE <- function(d2,d3,f1,d0,e0,whichgender, whichyear){
  d2 = d2*100
  d3 = d3*100
  f1 = f1*100
  d0 = d0*100
  e0 = e0*100
  d2$age<-c("18-25","26-34","35-49","50-64","65+","Total")
  d2$status<-"Current MD"
  d3$age<-c("18-25","26-34","35-49","50-64","65+","Total")
  d3$status<-"Former MD"
  f1$age<-c("18-25","26-34","35-49","50-64","65+","Total")
  f1$status<-"Recall error"
  f1 <- melt(f1[c(whichyear,"status","age")],by=c("status","age"))
  d0$age<-c("18-25","26-34","35-49","50-64","65+","Total")
  d0$status <- "Never MD"
  e0$age<-c("18-25","26-34","35-49","50-64","65+","Total")
  e0$status <- "Lifetime MD"
  e0 <- melt(e0[c(whichyear,"status","age")],by=c("status","age"))
  
  comparedata <- rbind(melt(d2[c(whichyear,"status","age")],by=c("status","age")),
                       melt(d3[c(whichyear,"status","age")],by=c("status","age")),
                       f1,
                       melt(d0[c(whichyear,"status","age")],by=c("status","age")))
  
  comparedata$status <- factor(comparedata$status, levels = rev(unique(comparedata$status))) # re-orders the status variable
  f1plot <- ggplot() + geom_bar(aes(y = value, x = age, fill = status), colour="black", data = comparedata, stat="identity") +
    geom_text(data=e0, aes(x = age, y = value,
                           label = paste0(sprintf("%.1f", round(value,3)),"%\n",
                                          "(",sprintf("%.1f",round(f1$value,3)),"%)")  ,vjust=-1), size=3) +
    # geom_text(data=f1, aes(x = age, y= value,
    #                        label = paste0("(f=",sprintf("%.1f", round(value,3)),"%)"),vjust=-10), size=3)+
    ggtitle(paste0(whichgender)) +
    theme(legend.position="bottom", legend.direction="horizontal",
          legend.title = element_blank(), plot.title = element_text(hjust = -0.14)) +
    scale_y_continuous(name="Prevalence (%)",limits=c(0,102),breaks=seq(0,100,10)) +
    scale_fill_manual(values=c("white","white","dark gray","black"))+
    labs(x="Age group")
  return(f1plot)
}
f1plotF <- lifetimeMDE(d2F, d3F, f1F, d0F, e0F,"A) Females", "2017") # Used photoshop to add cross-hatch to 'forgot' bars
f1plotM <- lifetimeMDE(d2M, d3M, f1M, d0M, e0M,"B) Males", "2017") 
jpeg(filename = paste0(mainDir, "Figure_4_lifetimeprev.jpg"),width=10, height=6, units ="in", res=1000)
grid_arrange_shared_legend(list(f1plotF,f1plotM),2,"")
dev.off()

## Supplement Figures
plotprevbyage <- function(df,status,nsduhstatus,subpopulation,lim,brks,whichgender){
  thissubset <- depsmkprevs_by_year[depsmkprevs_by_year$status==nsduhstatus & depsmkprevs_by_year$gender==whichgender & depsmkprevs_by_year$subpopulation==subpopulation & depsmkprevs_by_year$age!="total",]
  nsduhdata <- melt(thissubset, id.vars = c("age","survey_year"),measure.vars = c("prev"))
  nsduhdatalow <- melt(thissubset, id.vars = c("age","survey_year"),measure.vars = c("prev_lowCI"))
  nsduhdatahigh <- melt(thissubset, id.vars = c("age","survey_year"),measure.vars = c("prev_highCI"))
  colnames(nsduhdata)[4] <- "nsduh_prev"
  nsduhdata$nsduh_low <- nsduhdatalow$value
  nsduhdata$nsduh_high <- nsduhdatahigh$value
  
  df$age<-rownames(df)
  df <- melt(subset(df, df$age!="total"), id.vars=c("age"))
  colnames(df)[2] <- "survey_year"
  colnames(df)[3] <- "model_prev"
  comparedata <- merge(nsduhdata, df, by=c("age","survey_year"))
  comparedata$survey_year <- as.numeric(comparedata$survey_year)
  g <- ggplot(comparedata, aes(x=survey_year, y=model_prev*100, group=age)) + 
    geom_line(aes(x=survey_year, y= model_prev*100,colour=age ))+
    geom_pointrange(aes(x = survey_year, y = nsduh_prev*100, ymin=nsduh_low*100, ymax = nsduh_high*100, colour=age))  +
    scale_y_continuous(name="Prevalence (%)",limits=lim,breaks=brks) +
    scale_x_continuous(name="Year",limits=c(2005,endyear),breaks=seq(2005,2017,1))  +
    theme(axis.text.x=element_text(angle=45, hjust=1),plot.title = element_text(hjust = -0.14)) +
    labs(title=paste0(status," prev - ",subpopulation," ", whichgender),color="age") +
    ggtitle(bquote(atop(.(paste0(status," prevalence - ", capitalize(whichgender)))))) 
  return(g)
}

# Figure S1
e1plotF <- plotprevbyage(e1F, "everdep","everdep","totalpop",c(0,35),seq(0,35,5),"females") + ggtitle(paste0("A) Females")) + labs(color='Age group') 
e1plotM <- plotprevbyage(e1M, "everdep","everdep","totalpop",c(0,35),seq(0,35,5),"males") + ggtitle(paste0("B) Males")) + labs(color='Age group') 
jpeg(filename = paste0(mainDir, "Figure_S1_NOrecallerror.jpg"),width=10, height=6, units ="in", res=1000)
grid_arrange_shared_legend(list(e1plotF,e1plotM),2,"")
dev.off()

# Figure S2
e0plotF <- plotprevbyage(e0F, "everdep","everdep","totalpop",c(0,35),seq(0,35,5),"females")+ ggtitle(paste0("A) Females")) + labs(color='Age group')
e0plotM <-plotprevbyage(e0M, "everdep","everdep","totalpop",c(0,35),seq(0,35,5),"males")+ ggtitle(paste0("B) Males")) + labs(color='Age group')
jpeg(filename = paste0(mainDir, "Figure_S2_recallerror.jpg"),width=10, height=6, units ="in", res=1000)
grid_arrange_shared_legend(list(e0plotF,e0plotM),2,"")
dev.off()
# Figure S4
Sf1plotF <- lifetimeMDE(d2F, d3F, f1F, d0F, e0F,"A) Females", "2015") 
Sf1plotM <- lifetimeMDE(d2M, d3M, f1M, d0M, e0M,"B) Males", "2015") 
jpeg(filename = paste0(mainDir, "Figure_S4_lifetimeprev.jpg"),width=10, height=6, units ="in", res=1000)
grid_arrange_shared_legend(list(Sf1plotF,Sf1plotM),2,"")
dev.off()

# Numbers for manuscript maintext 
test = subset(depsmkprevs_by_year, gender=="females" & status=="everdep" & subpopulation=="totalpop" & age!="total" & age!="65plus")
min(test$prev)
max(test$prev)

test = subset(depsmkprevs_by_year, gender=="males" & status=="everdep" & subpopulation=="totalpop" & age!="total" & age!="65plus")
min(test$prev)
max(test$prev)

test = subset(depsmkprevs_by_year, gender=="females" & status=="everdep" & subpopulation=="totalpop" & age=="65plus")
min(test$prev)
max(test$prev)


e1F["2017"]
e1M["2017"]

d2F["2017"]
d2M["2017"]

e0F["2017"]
e0M["2017"]

sum(everdeppoptrueM[c(19:100),"2017" ] + everdeppoptrueF[c(19:100),"2017" ]) / sum(totalpopM[c(19:100),"2017" ] + totalpopF[c(19:100),"2017" ])