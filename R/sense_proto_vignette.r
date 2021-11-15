# Ecosense vignette

library('Rpath')
library('viridis')

data("sense_vignette")
source("R/ecosense.R")
# NOTE: The depricated sense code is still included in ecosense.R (lines 211:401) for testing purposes.

# One generated system (params)
# one_ecosys <- rsim.sense(EBS_scene,EBS_unbal,Vvary = c(-4.5,4.5), Dvary = c(0,0))

# Set up generator loop
EBS_scene$params$BURN_YEARS <- 50 # burn-in period
NUM_RUNS <- 2000 # how many systems to generate
parlist<-as.list(rep(NA,NUM_RUNS)) # create lists to store generated parameters
kept<-rep(NA,NUM_RUNS) # object to keep track of kept systems
set.seed(666) # Optional, set seed so output can be replicated

for (i in 1:NUM_RUNS){
  EBSsense <- EBS_scene # scenario object
  # INSERT SENSE ROUTINE BELOW
  parlist[[i]]<- EBS_scene$params 		# Base ecosim params
  parlist[[i]]<- rsim.sense(EBS_scene,EBS_unbal,Vvary = c(-4.5,4.5), Dvary = c(0,0))	# Replace the base params with Ecosense params
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

# KEPT tells you which ecosystems were kept
KEPT <- which(kept==T) # the number associated with the kept system
nkept <- length(KEPT); nkept # how many were kept
1-(nkept/NUM_RUNS) # rejection rate

# Running rsim.run with the kept systems (no perturbation)
EBS_ecos <- as.list(rep(NA,length(KEPT))) # lists for simulated ecosystems
k <- 0  # counter for simulated ecosystems
for (i in KEPT) {
  EBS_scene$start_state$Biomass <- parlist[[i]]$B_BaseRef # set the starting Biomass to the generated values
  EBSsense <- EBS_scene # set up the scenario object
  EBSsense$params <- parlist[[i]] # set the params in the scenario object equal to the generated params.
  EBSsense$BURN_YEARS <- -1 # no burn-in period
  k <- k + 1 # set the number for the simulated ecosystem
  EBS_ecos[[k]] <- rsim.run(EBSsense,method='AB') # run rsim.run on the generated system
  print(c("Ecosystem no.", k, "out of", nkept)) # progress output to console
}

# Calculate relative biomass at new equilibrium (month 1200) versus the drawn value (t=1)
relB_EBS_ecos <- as.list(rep(NA,length(KEPT))) # list to output relative biomass
k <- 0
for (i in 1:nkept) {
  spname <- colnames(EBS_ecos[[i]]$out_Biomass[,2:ncol(EBS_ecos[[i]]$out_Biomass)])
  biomass <- EBS_ecos[[i]]$out_Biomass[, spname]
  n <- ncol(biomass)
  start.bio <- biomass[1, ] # the drawn starting biomass 
  start.bio[which(start.bio == 0)] <- 1
  rel.bio <- matrix(NA, dim(biomass)[1], dim(biomass)[2])
  for(isp in 1:n) rel.bio[, isp] <- biomass[, isp] / start.bio[isp] # biomass relative to biomass at t=1
  colnames(rel.bio) <- spname
  k <- k + 1
  relB_EBS_ecos[[k]] <- rel.bio
}

# Example plot of pollock
plot_mat <- matrix(nrow=1200, ncol=nkept) # matrix of pollock trajectories from all generated systems
for (i in 1:nkept) {
  plot_mat[,i] <- relB_EBS_ecos[[i]][,"Walleye pollock"]
}
plot_col <- viridis(nkept)
# par(mfrow=c(1,2))
layout(matrix(c(1,1,1,1,1,1,1,2), nrow = 1, ncol = 8, byrow = TRUE))
plot(1:1200, relB_EBS_ecos[[1]][,"Walleye pollock"], type='n', xlab="Months",
     ylab="Relative biomass", ylim=c(min(plot_mat),max(plot_mat)), main="Walleye pollock")
# one line for pollock in each of the generated systems
for (i in 1:nkept) {
  lines(1:1200, relB_EBS_ecos[[i]][,"Walleye pollock"], lwd=2, col=plot_col[i])
}
# distribution of pollock relative biomasses
boxplot(plot_mat[1200,], ylim=c(min(plot_mat),max(plot_mat)), yaxt='n')
axis(side=2, at=c(seq(0,150,50)), tick=TRUE, labels=F)


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Now with a perturbation
# Increase FRate on pollock
# Running rsim.run with the kept systems (no perturbation)
EBS_ecos_pol <- as.list(rep(NA,length(KEPT))) # lists for simulated ecosystems
k <- 0  # counter for simulated ecosystems
for (i in KEPT) {
  EBS_scene$start_state$Biomass <- parlist[[i]]$B_BaseRef # set the starting Biomass to the generated values
  EBSsense <- EBS_scene # set up the scenario object
  EBSsense <- adjust.fishing(EBSsense, "ForcedFRate", group="Walleye pollock", sim.year=51:100, value=2) # perturb pollock FRate
  EBSsense$params <- parlist[[i]] # set the params in the scenario object equal to the generated params.
  EBSsense$BURN_YEARS <- -1 # no burn-in period
  k <- k + 1 # set the number for the simulated ecosystem
  EBS_ecos_pol[[k]] <- rsim.run(EBSsense,method='AB') # run rsim.run on the generated system
  print(c("Ecosystem no.", k, "out of", nkept)) # progress output to console
}

# Calculate relative biomass at month 1200 relative to the end of the burn-in (t=600)
relB_EBS_ecos_pol <- as.list(rep(NA,length(KEPT)))
k <- 0
for (i in 1:nkept) {
  spname <- colnames(EBS_ecos_pol[[i]]$out_Biomass[,2:ncol(EBS_ecos_pol[[i]]$out_Biomass)])
  biomass <- EBS_ecos_pol[[i]]$out_Biomass[, spname]
  n <- ncol(biomass)
  start.bio <- biomass[600, ] # end of burn-in
  start.bio[which(start.bio == 0)] <- 1
  rel.bio <- matrix(NA, dim(biomass)[1], dim(biomass)[2])
  for(isp in 1:n) rel.bio[, isp] <- biomass[, isp] / start.bio[isp] # biomass relative to the end of burn-in biomass
  colnames(rel.bio) <- spname
  k <- k + 1
  relB_EBS_ecos_pol[[k]] <- rel.bio
}

# Example plot of pollock
# ecos_n <- 1:nkept
plot_mat_pol <- matrix(nrow=1200, ncol=nkept)
for (i in 1:nkept) {
  plot_mat_pol[,i] <- relB_EBS_ecos_pol[[i]][,"Walleye pollock"]
}
plot_col <- viridis(nkept)
# par(mfrow=c(1,2))
layout(matrix(c(1,1,1,1,1,1,1,2), nrow = 1, ncol = 8, byrow = TRUE))
plot(601:1200, relB_EBS_ecos_pol[[1]][601:1200,"Walleye pollock"], type='n', xlab="Months",
     ylab="Relative biomass", ylim=c(min(plot_mat_pol[601:1200,]),max(plot_mat_pol[601:1200,])), main="Walleye pollock")
for (i in 1:nkept) {
  lines(601:1200, relB_EBS_ecos_pol[[i]][601:1200,"Walleye pollock"], lwd=2, col=plot_col[i])
}
boxplot(plot_mat_pol[1200,], ylim=c(min(plot_mat_pol[601:1200,]),max(plot_mat_pol[601:1200,])))


#------------------------------------------------------------------------------#
# Boxplot of multiple species
# matrix to store the end biomass for all species (columns) from all the generated
# ecosystems (rows).
pol_perturb_out <- matrix(nrow=nkept, ncol=(EBS_bal$NUM_LIVING+EBS_bal$NUM_DEAD))
for(i in 1:nkept){
  pol_perturb_out[i,] <- relB_EBS_ecos_pol[[i]][1200,]
}
colnames(pol_perturb_out) <- colnames(relB_EBS_ecos_pol[[1]])
par(mfrow=c(1,1), mar=c(9,3,1,1))
boxplot(pol_perturb_out, outline=FALSE, las=2, cex.axis=0.8)
abline(h=1, lty=2)
