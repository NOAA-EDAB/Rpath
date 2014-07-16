################################################################################
# EXAMPLES OF USING ECOPATH AND ECOSIM in R and C 
# written by Kerim Aydin and Sarah Gaichas
# v0.04 3/12/13
################################################################################

# Directory with all the files
  setwd("C:/Users/kaydin/d/src/sims/sim_v04")

# R code files with Ecopath and Ecosim routines
  source("ecopath_r_v0_04.r")
  source("ecosim_r_v0_04.r")
  
# Load ecosim dll 
  ll() 

# Ecopath input parameters are kept in csv files.  See examples for format.
# Note: there is very limited checking to ensure files have the right number
# of columns, likely outcome of poor csv files is an error.

# Example is a reduced version of the Eastern Bering Sea (Aydin et al. 2007)

  # Basic Ecopath parameters (Biomass, Production, Consumption, Fishing)
    BaseFile <- "data/EBS_andre_base.csv"
  # Ecopath Diet Matrix
    DietFile <- "data/EBS_andre_diet.csv"
  # Data pedigree (if no pedigree, use file of shown format, with 1.0 for all values) 
    PedFile  <- "data/EBS_andre_ped.csv"
  # Juvenile adult split groups.  IF NO JUV/ADU IN YOUR MODEL, use file with the 
  # headers in the example, but no rows.  This is required as a placeholder 
  # (minor bug needs fixing)
    JuvFile  <- "data/EBS_andre_juvs.csv"


# EXAMPLES

# EXAMPLE 1:  BALANCE ECOPATH, SENDING RESULTS TO NEW CSV FILE and set of vectors
#             (leave filename blank or FALSE for vectors only)
  
  path <- ecopathR(BaseFile,DietFile,PedFile,"EBS_balanced")  
  summary(path) 
  
# EXAMPLE 2:  Load Ecosim parameters, declaring enough vector length for 100 years.
#             Run once in baseline (equilibrium), run with fishing change.

  # Prepare for run by loading from csv files (includes balancing Ecopath)
    base_sim <- load_ecosim( Years = 100, BaseFile, DietFile, PedFile, JuvFile)
    #long_sim <- load_ecosim( Years = 500, BaseFile, DietFile, PedFile, JuvFile)
  # Base run from year 0 to year 100
    run0 <- ecosim_run(base_sim,0,100)  

  # resulting output biomass
	summary(run0$out_BB)

  # Set up a fishing scenario, set a high F Rate on Pollock (note: this doesn't 
  # remove the fishing on gears) for 90 years.  
	fish_sim <- base_sim 	
	fish_sim$FORCED_FRATE$'W. Pollock_Adu'[10:101] <- 0.6 
	run1 <- ecosim_run(fish_sim,0,100) 
	
  # Now add some lognormal environmental forcing to mortality (monthly from year 10)
    #fish_sim$force_bymort$'W. Pollock_Adu'[120:1201] <- rlnorm(1082,0,1)
    #run2 <- ecosim_run(fish_sim,0,100)
	 
  # Recruitment variabilitiy
    fish_sim <- base_sim   
    fish_sim$FORCED_FRATE$'W. Pollock_Adu'[10:101] <- 0.6   
    fish_sim$force_byrecs$'W. Pollock_Juv'[122:1201] <- (rep(rlnorm(90,0,.1),each=12))
    run2 <- ecosim_run(fish_sim,0,100)
      
  # Run 10 years at base, then remove all fishing from ecosystem by setting Q for all gear to 0
  # and run for another 90 years
    #fish_sim <- base_sim #(resetting to base)
    #run3a <- ecosim_run(fish_sim,0,9)
    #run3a$fishlist$FishQ[] <- 0.0 
    #run3 <-ecosim_run(run3a,9,100)
    fish_sim <- base_sim
    id <- base_sim$spnum; names(id) <-  base_sim$spname

    fish_sim$fishlist$FishQ[which(base_sim$fishind$FishFrom == id['W. Pollock_Adu'])] <- 0.0   
    #fish_sim$FORCED_FRATE$'W. Pollock_Adu'[32:501] <- seq(0,1,length.out=470)
    #fish_sim$FORCED_FRATE$'W. Pollock_Adu'[32:501] <- 0.4
    fish_sim$targlist$TARGET_BIO[3] <- 11.31272
    fish_sim$targlist$TARGET_F[3]   <- 0.4
    fish_sim$targlist$ALPHA[3]      <- 0.2 
    cann <- which(base_sim$ppind$PreyFrom == id['W. Pollock_Juv'] & base_sim$ppind$PreyTo == id['W. Pollock_Adu'])
    fish_sim$pplist[cann,]$VV <- 10 
    fish_sim$force_byrecs$'W. Pollock_Juv'[2:1201] <- (rep(rlnorm(length(2:1201),0,2),each=12))
    run3 <- ecosim_run(fish_sim,0,100)


  # finally, plot output from all the runs for pollock and cod
  par(mfrow=c(2,2))
  
    plot(1:1201,run0$out_BB$'W. Pollock_Adu',ylab="Pollock biomass (t/km2)",xlab="month",type='l', ylim=c(0,30))
   lines(1:1201,run1$out_BB$'W. Pollock_Adu',col="red")  
   lines(1:1201,run2$out_BB$'W. Pollock_Adu',col="blue")  
   lines(1:1201,run3$out_BB$'W. Pollock_Adu',col="green")

    plot(1:1201,run0$out_BB$'P. Cod_Adu',ylab="Cod biomass (t/km2)",xlab="month",type='l', ylim=c(0,10))
   lines(1:1201,run1$out_BB$'P. Cod_Adu',col="red")  
   lines(1:1201,run2$out_BB$'P. Cod_Adu',col="blue")  
   lines(1:1201,run3$out_BB$'P. Cod_Adu',col="green")   
   
    plot(1:1201,run0$out_BB$'Arrowtooth_Adu',ylab="Arrowtooth biomass (t/km2)",xlab="month",type='l', ylim=c(0,2))
   lines(1:1201,run1$out_BB$'Arrowtooth_Adu',col="red")  
   lines(1:1201,run2$out_BB$'Arrowtooth_Adu',col="blue")  
   lines(1:1201,run3$out_BB$'Arrowtooth_Adu',col="green")   
   
   par(mfrow=c(1,2))
   nn <- seq(6,6001,12)
   plot(1:1201,run3$out_BB$'W. Pollock_Adu',ylab="Pollock biomass (t/km2)",xlab="month",type='l', ylim=c(0,35))

   plot(run3$out_SSB$'W. Pollock_Adu'[nn],run3$out_rec$'W. Pollock_Adu'[nn])
   points(run3$out_SSB$'W. Pollock_Adu'[1],run3$out_rec$'W. Pollock_Adu'[1],col="red")
   
   plot(log(run3$out_SSB$'W. Pollock_Adu'[nn]),log(run3$out_rec$'W. Pollock_Adu'[nn]))
   
   plot(run3$out_BB$'W. Pollock_Adu',
        12*run3$out_CC$'W. Pollock_Adu'/run3$out_BB$'W. Pollock_Adu') 
   
      