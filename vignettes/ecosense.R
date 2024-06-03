## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(Rpath); library(data.table); library(viridis)
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())

## ----load the unbalanced models, balance the models, setup rsim scenario objects----
# load the unbalanced models
load("data/Ecosense.EBS.rda")
load("data/Ecosense.ECS.rda")
load("data/Ecosense.GOA.rda")
# balance the models
EBS_bal <- rpath(Ecosense.EBS)
ECS_bal <- rpath(Ecosense.ECS)
GOA_bal <- rpath(Ecosense.GOA)
# create rsim scenario objects
EBS_scene <- rsim.scenario(EBS_bal, Ecosense.EBS, years=1:100)
ECS_scene <- rsim.scenario(ECS_bal, Ecosense.ECS, years=1:100)
GOA_scene <- rsim.scenario(GOA_bal, Ecosense.GOA, years=1:100)
# source ecosense.R
source("R/ecosense.R", local = knitr::knit_global())
ls()

## ----echo=TRUE----------------------------------------------------------------
# One set of Ecosim parameters for the EBS model
rsim.sense(EBS_scene,Ecosense.EBS,Vvary = c(0,0), Dvary = c(0,0))

## ----echo=TRUE----------------------------------------------------------------
# Setting the burn-in period in the EBS scenario object to 50 years.
EBS_scene$params$BURN_YEARS <- 50

## ----echo=TRUE----------------------------------------------------------------
NUM_RUNS <- 1000 # how many ecosystem parameter sets to generate
parlist<-as.list(rep(NA,NUM_RUNS)) # create lists to store the generated parameters
kept<-rep(NA,NUM_RUNS) # object to keep track of kept systems
set.seed(666) # Optional, set seed so output can be replicated

## ----generator loop, echo=TRUE, results='hide'--------------------------------
for (i in 1:NUM_RUNS){
  EBSsense <- EBS_scene # scenario object
  # INSERT SENSE ROUTINE BELOW
  parlist[[i]]<- EBS_scene$params 		# Base ecosim params
  parlist[[i]]<- rsim.sense(EBS_scene,Ecosense.EBS,Vvary = c(-4.5,4.5), Dvary = c(0,0))	# Replace the base params with Ecosense params
  EBSsense$start_state$Biomass <- parlist[[i]]$B_BaseRef # Apply the Ecosense starting biomass
  parlist[[i]]$BURN_YEARS <- 50			# Set Burn Years to 50
  EBSsense$params <- parlist[[i]] # replace base params with the Ecosense generated params
  EBStest <- rsim.run(EBSsense, method="AB") # Run rsim with the generated system
  failList <- which(is.na(EBStest$end_state$Biomass))
  {if (length(failList)>0)
  {cat(i,": fail in year ",EBStest$crash_year,": ",EBStest$params$spname[failList],"\n"); kept[i]<-F; flush.console()}
    else 
    {cat(i,": success!\n"); kept[i]<-T;  flush.console()}} # output for the console
  parlist[[i]]$BURN_YEARS <- 1
}

## ----echo=TRUE----------------------------------------------------------------
KEPT <- which(kept==T); KEPT # the number associated with the kept system
nkept <- length(KEPT); nkept # how many were kept
1-(nkept/NUM_RUNS) # rejection rate

## ----echo=TRUE, results='hide'------------------------------------------------
ecos <- as.list(rep(NA,length(KEPT))) # lists for simulated ecosystems
k <- 0  # counter for simulated ecosystems
for (i in KEPT) {
  EBS_scene$start_state$Biomass <- parlist[[i]]$B_BaseRef # set the starting Biomass to the generated values
  EBSsense <- EBS_scene # set up the scenario object
  EBSsense$params <- parlist[[i]] # set the params in the scenario object equal to the generated params.
  EBSsense$BURN_YEARS <- -1 # no burn-in period
  k <- k + 1 # set the number for the simulated ecosystem
  ecos[[k]] <- rsim.run(EBSsense,method='AB') # run rsim.run on the generated system
  print(c("Ecosystem no.",k,"out of",nkept)) # progress output to console
}

## ----echo=TRUE----------------------------------------------------------------
relB_ecos <- as.list(rep(NA,length(KEPT))) # list to output relative biomass
k <- 0
for (i in 1:nkept) {
  spname <- colnames(ecos[[i]]$out_Biomass[,2:ncol(ecos[[i]]$out_Biomass)])
  biomass <- ecos[[i]]$out_Biomass[, spname]
  n <- ncol(biomass)
  start.bio <- biomass[1, ] # the drawn starting biomass 
  start.bio[which(start.bio == 0)] <- 1
  rel.bio <- matrix(NA, dim(biomass)[1], dim(biomass)[2])
  for(isp in 1:n) rel.bio[, isp] <- biomass[, isp] / start.bio[isp] # biomass relative to biomass at t=1
  colnames(rel.bio) <- spname
  k <- k + 1
  relB_ecos[[k]] <- rel.bio
}

## ----echo=TRUE----------------------------------------------------------------
this_species <- "Walleye pollock"
plot_mat <- matrix(nrow=1200, ncol=nkept) # matrix of pollock trajectories from all generated systems
for (i in 1:nkept) {
  plot_mat[,i] <- relB_ecos[[i]][,this_species]
}

## ----echo=TRUE----------------------------------------------------------------
plot_col <- viridis(nkept)
layout(matrix(c(1,1,1,1,1,1,2,2), nrow = 1, ncol = 8, byrow = TRUE))
plot(1:1200, relB_ecos[[1]][,this_species], type='n', xlab="Months",
     ylab="Relative biomass", ylim=c(min(plot_mat),max(plot_mat)), main=this_species)
# one line for pollock in each of the generated systems
for (i in 1:nkept) {
  lines(1:1200, relB_ecos[[i]][,this_species], lwd=2, col=plot_col[i])
}
# distribution of pollock relative biomasses
boxplot(plot_mat[1200,], ylim=c(min(plot_mat),max(plot_mat)), yaxt='n')
axis(side=2, at=c(seq(0,150,50)), tick=TRUE, labels=F)

## ----echo=TRUE, results='hide'------------------------------------------------
ecos_sp <- as.list(rep(NA,length(KEPT))) # lists for simulated ecosystems
k <- 0  # counter for simulated ecosystems
for (i in KEPT) {
  EBS_scene$start_state$Biomass <- parlist[[i]]$B_BaseRef # set the starting Biomass to the generated values
  EBSsense <- EBS_scene # set up the scenario object
  EBSsense <- adjust.fishing(EBSsense, "ForcedFRate", group=this_species, sim.year=51:100, value=2) # perturb pollock FRate
  EBSsense$params <- parlist[[i]] # set the params in the scenario object equal to the generated params.
  EBSsense$BURN_YEARS <- -1 # no burn-in period
  k <- k + 1 # set the number for the simulated ecosystem
  ecos_sp[[k]] <- rsim.run(EBSsense,method='AB') # run rsim.run on the generated system
  print(c("Ecosystem no.", k, "out of", nkept)) # progress output to console
}

## ----echo=TRUE----------------------------------------------------------------
relB_ecos_sp <- as.list(rep(NA,length(KEPT)))
k <- 0
for (i in 1:nkept) {
  spname <- colnames(ecos_sp[[i]]$out_Biomass[,2:ncol(ecos_sp[[i]]$out_Biomass)])
  biomass <- ecos_sp[[i]]$out_Biomass[, spname]
  n <- ncol(biomass)
  start.bio <- biomass[600, ] # end of burn-in
  start.bio[which(start.bio == 0)] <- 1
  rel.bio <- matrix(NA, dim(biomass)[1], dim(biomass)[2])
  for(isp in 1:n) rel.bio[, isp] <- biomass[, isp] / start.bio[isp] # biomass relative to the end of burn-in biomass
  colnames(rel.bio) <- spname
  k <- k + 1
  relB_ecos_sp[[k]] <- rel.bio
}

## ----echo=TRUE----------------------------------------------------------------
plot_mat_sp <- matrix(nrow=1200, ncol=nkept)
for (i in 1:nkept) {
  plot_mat_sp[,i] <- relB_ecos_sp[[i]][,this_species]
}
layout(matrix(c(1,1,1,1,1,1,2,2), nrow = 1, ncol = 8, byrow = TRUE))
plot(601:1200, relB_ecos_sp[[1]][601:1200,this_species], type='n', xlab="Months",
     ylab="Relative biomass", ylim=c(min(plot_mat_sp[601:1200,]),max(plot_mat_sp[601:1200,])), main=this_species)
for (i in 1:nkept) {
  lines(601:1200, relB_ecos_sp[[i]][601:1200,this_species], lwd=2, col=plot_col[i])
}
boxplot(plot_mat_sp[1200,], ylim=c(min(plot_mat_sp[601:1200,]),max(plot_mat_sp[601:1200,])))

## ----echo=TRUE----------------------------------------------------------------
sp_perturb_out <- matrix(nrow=nkept, ncol=(EBS_bal$NUM_LIVING+EBS_bal$NUM_DEAD))
for(i in 1:nkept){
  sp_perturb_out[i,] <- relB_ecos_sp[[i]][1200,]
}
colnames(sp_perturb_out) <- colnames(relB_ecos_sp[[1]])
par(mfrow=c(1,1), mar=c(9,3,0.5,0.5))
boxplot(sp_perturb_out, outline=FALSE, las=2, cex.axis=0.6)
abline(h=1, lty=2)

