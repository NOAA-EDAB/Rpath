
setwd("C:/users/kaydin/d/src/sims/sim_r")

# Load functions
  source("ecopathR_v3.r")
  source("ecosimR_v2.r")
  source("ecosim_plots_v1.r")
  require(fields)
# Load ecosim dll; unload is ul()
  ll()

# refresh function
  sr<-function(){source("sim_wrap_tests.r")}

#PF <- "ai/AI_2011_base.csv"
#DF <- "ai/AI_2011_diet.csv"
#PD <- "ai/AI_2011_ped.csv"
#JF <- "ai/AI_2011_stanzas_in.csv"

PF <- "data/EBS_full_base.csv"
DF <- "data/EBS_full_diet.csv"
PD <- "data/EBS_full_ped.csv"
JF <- "data/EBS_full_juvs.csv"


###############################################################################

test_run <- function(){
 
 base   <- load_ecosim( Years = 100, PF, DF, PD, JF, FALSE)
 runs   <- 10
 
 # Thresholds for discards; changing these here only changes what the R code 
 # prints out; these numbers are hardcoded in the C code
   LO_DISCARD <- 1e-4
   HI_DISCARD <- 1e+4
 
# lists to store kept and failed ecosystems
  keptruns <- as.list(NULL)
  failruns <- as.list(NULL)
  nkept <- 0
  nfail <- 0
 
 for (i in 1:runs){
    baseN <- base
    
    # The BURN_YEARS flag tells the C code how long the "burn-in" period is for.
    # If set to -1, the model won't stop if a species dies out.  If set to N,
    # model will stop if any species crashes in the first N years.  After N
    # years, it won't stop running if a species crashes.  NOTE:  the model will
    # also stop running if it encounteras an NaN or inf value, even if it isn't
    # during burn-in. 
      baseN$flags["BURN_YEARS"] <- 50
    
    # perturb Mzero randomly
      baseN$ratelist$MzeroMort <- baseN$ratelist$MzeroMort * (0.5 + 2*runif(length(baseN$ratelist$MzeroMort)))
    
    # run model
      run    <- ecosim_run(baseN)
    
    # This makes a list of which species went below/above thresholds
      ystate    <- run$state$state_BB / run$ratelist$B_BaseRef 
      crashtest <- (run$state$state_BB < run$ratelist$B_BaseRef * LO_DISCARD) | 
                 (run$state$state_BB > run$ratelist$B_BaseRef * HI_DISCARD)
      crashlist <- run$spname[crashtest]
      crashvals <- ystate[crashtest]
    
    # Print status (see example below)
      print(paste("run",i," crash year",run$outflags["CRASH_YEAR"],paste(crashlist,crashvals,collapse=", "))); flush.console()
    
    # Save as failed or kept run.  If run$outflags["CRASH_YEAR"] is -1, the 
    # model didn't crash during burn-in.
      if  (run$outflags["CRASH_YEAR"] < 0){
          nkept <- nkept+1
          keptruns[[nkept]] <- run
       } else {
          nfail <- nfail+1
          failruns[[nfail]] <- run                
       }
  }

} 

## Sample output 
# If the model crashes, this output tells you what species crashed and in what
# year.  The number is Bstate/Bstart.  If the model was "kept" it puts out -1 
# for year.  Note that if it prints out -1 but still prints out some species,
# it means that the species crashed after the burn-in period, but the model was
# still kept.  Use this to adjust burn-in if desired (or just see what's on
# its way out)  

# [1] "run 1  crash year 25 Outside Production 9.77003597119925e-05"
# [1] "run 2  crash year -1 Kamchatka fl._Juv 2.16005893413678e-07, Kamchatka fl. 1.5243797285815e-07, Lg. Sculpins 4.83205444190362e-05"
# [1] "run 3  crash year 32 Sea stars 9.96737312913438e-05"
# [1] "run 4  crash year 15 Outside Production 9.91036327958106e-05"
# [1] "run 5  crash year 15 Sea stars 9.62786767374209e-05"
# [1] "run 6  crash year 15 Outside Production 9.65384305304319e-05"
# [1] "run 7  crash year 27 Atka mackerel_Juv 9.82014784331025e-05"
# [1] "run 8  crash year 36 Atka mackerel_Juv 9.92297160002852e-05"
# [1] "run 9  crash year -1 "
# [1] "run 10  crash year -1 Kamchatka fl._Juv 8.16440623297147e-05, Kamchatka fl. 3.79313804990171e-05, Anemones 1.28865443876152e-07" 

## keptruns now holds a list of the kept ecosystems (in this case runs 2,9, and 10 above)
# > summary(keptruns)
#      Length Class  Mode
# [1,] 29     -none- list
# [2,] 29     -none- list
# [3,] 29     -none- list

## Each one of these is one full run including parameters 
# > summary(keptruns[[3]])
#                Length Class      Mode     
# YEARS             1   -none-     numeric  
# spname          150   -none-     character
# spnum           150   -none-     numeric  
# numlist           8   -none-     numeric  
# flags             1   -none-     numeric  
# outflags          1   -none-     numeric  
# ratelist         13   data.frame list     
# ppind             2   data.frame list     
# pplist            6   data.frame list     
# fishind           3   data.frame list     
# fishlist          1   data.frame list     
# detind            2   data.frame list     
# detlist           1   data.frame list     
# juvind            6   data.frame list     
# juvlist          26   data.frame list     
# statelist         5   data.frame list     
# NageS          6015   -none-     numeric  
# WageS          6015   -none-     numeric  
# WWa            6015   -none-     numeric  
# SplitAlpha     6015   -none-     numeric  
# derivlist        23   data.frame list     
# force_byprey    150   data.frame list     
# force_bymort    150   data.frame list     
# force_byrecs    150   data.frame list     
# force_bysearch  150   data.frame list     
# FORCED_FRATE    150   data.frame list     
# FORCED_CATCH    150   data.frame list     
# out_BB          150   data.frame list     
# out_CC          150   data.frame list     

## Here's the output biomass from one of them
# > summary(keptruns[[3]]$out_BB)
#     Outside  Transient Killers   Sperm and Beaked Whales Resident Killers      Porpoises           Belugas          Gray Whales        Humpbacks           Fin Whales    
#  Min.   :0   Min.   :0.0001351   Min.   :0.0002293       Min.   :0.0002255   Min.   :0.001172   Min.   :0.009806   Min.   :0.01881   Min.   :0.0004649   Min.   :0.4547  
#  1st Qu.:0   1st Qu.:0.0001691   1st Qu.:0.0005298       1st Qu.:0.0003117   1st Qu.:0.001438   1st Qu.:0.010125   1st Qu.:0.01979   1st Qu.:0.0008737   1st Qu.:0.4840  
#  Median :0   Median :0.0001993   Median :0.0013139       Median :0.0004581   Median :0.001821   Median :0.010544   Median :0.02164   Median :0.0017501   Median :0.5016  
#  Mean   :0   Mean   :0.0001969   Mean   :0.0030267       Mean   :0.0005581   Mean   :0.001989   Mean   :0.010700   Mean   :0.02302   Mean   :0.0029165   Mean   :0.4969  
#  3rd Qu.:0   3rd Qu.:0.0002259   3rd Qu.:0.0038431       3rd Qu.:0.0007392   3rd Qu.:0.002426   3rd Qu.:0.011153   3rd Qu.:0.02542   3rd Qu.:0.0039968   3rd Qu.:0.5124  
#  Max.   :0   Max.   :0.0002486   Max.   :0.0176747       Max.   :0.0013502   Max.   :0.003619   Max.   :0.012255   Max.   :0.03267   Max.   :0.0120970   Max.   :0.5187  
# [etc.]