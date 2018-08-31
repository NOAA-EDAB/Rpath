  
library(Rpath)

# Files containing initial ecosystem data
  Ebase <- "data/EBS_andre_base.csv"  # Base biomass, production, fishing
  Ediet <- "data/EBS_andre_diet.csv"  # Diet matrix
  Eped  <- "data/EBS_andre_ped.csv"   # Data pedigree = quality of input data
  Ejuv  <- "data/EBS_andre_juvs.csv"  # Juvenile/Adult file if needed

# Create an ECOPATH (balanced linear equation) model
  EBS   <- ecopath(Ebase, Ediet, Eped, eco.name = 'E. Bering')

# Create an ECOSIM SCENARIO (TBASE) from an input balanced ecosystem. 
# A Scenario created by ecosim.init is a List of Lists, consisting of:
# params:  
#     List containing the parameters of the system (i.e. the
#     parameter set that defines a particular ecosystem).
# start_state:
#     a List containing start values of state variables.
# forcing:
#     a List of Month x Species matrices (default value 1.0) that
#     are used to apply "environmental forcing" in monthly timesteps.
# fishing:
#     a List of Year x Species matrices for fishing (annual targets
#     for catch, or effort)
  TBASE <- Rsim.scenario(EBS)

# Manually change the fishing rate for species 2, years 1:30
  TBASE$fishing$FRATE[1:30,2]<-0.05

# Run for 50 years using Adams-Basforth, plot result  
  TRUNa <- adams.run(TBASE,50)
  plot(TRUNa$out_BB[,2])

# Run for 50 years Using Rk4, setting stepsize = 2/month, and plot
  TBASE$params$RK4_STEPS <- 2
  TRUNr <- rk4.run(TBASE,50)
  plot(TRUNr$out_BB[,2])

  
# OLDER STUFF (not guaranteed to work)
  require(microbenchmark)
  microbenchmark(adams.run(TBASE,100),times=100L)  
  
  


  








