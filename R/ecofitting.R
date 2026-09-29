################################################################################
# KYA REPLACED WITH read.fitting.biomass.cv version below
# read.fitting.biomass <- function(Rsim.scenario, cdat){
#   
#   # TODO KYA 7/24/24: add intelligent warning if missing column (e.g. Type)
#   
#   # Base variables
#   SIM   <- Rsim.scenario
#   years <- as.numeric(row.names(Rsim.scenario$fishing$ForcedFRate))
#   species <- SIM$params$spname
#   
#   
#   #cdat  <-read.csv(filename)
#   missing_sp <- unique(cdat$Group[which(!(cdat$Group %in% species))])
#   if(length(missing_sp)>0){
#     warning("Following species in BIOMASS fit data not in model, dropped: ",missing_sp)
#   }
#   
#   # drop Groups not found in model
#   cmdat <- cdat[cdat$Group %in% species,]
#   # drop lines with NAs in Value or Stdev or Scale
#   c0dat <- cmdat[!is.na(cmdat$Value) & !is.na(cmdat$Stdev) & !is.na(cmdat$Scale),]
#   # drop lines not in model run years
#   c1dat <- c0dat[c0dat$Year %in% years,]
#   # drop lines where Value, StDev or Scale are <=0
#   ccdat <- c1dat[c1dat$Value>0 & c1dat$Stdev>0 & c1dat$Scale>0,]
#   # replace NAs in type with index (maybe with future warning)
#   ccdat$Type[is.na(ccdat$Type)] <- "index"
#   
#   #type <- as.character(rep("absolute",length(ccdat$YEAR)))
#   ccdat$Year  <- as.character(ccdat$Year)
#   ccdat$Group <- as.character(ccdat$Group)
#   ccdat$Source <- as.character(ccdat$Source)
#   ccdat$ID    <- paste(as.character(ccdat$Group),as.character(ccdat$Source),sep=":")
#   
#   obs  <- ifelse(as.numeric(ccdat$Scale)<0, as.numeric(ccdat$Value),
#                  as.numeric(ccdat$Value) * as.numeric(ccdat$Scale))   
#   sd   <- ifelse(as.numeric(ccdat$Scale)<0, as.numeric(ccdat$Stdev),
#                  as.numeric(ccdat$Stdev) * as.numeric(ccdat$Scale))   
#   wt   <- rep(1,length(obs))
#   initial_q   <- rep(1,length(obs))
#   SIM$fitting$Biomass <- cbind(ccdat,obs,sd,initial_q,wt)
#   
#   return(SIM)
#   
# }

################################################################################
#' Read and Process Biomass Fitting Data
#'
#' Reads, cleans, and formats biomass observation data to prepare it for model 
#' fitting within an `Rsim` scenario object.
#'
#' @details
#' Observation data must be provided in the `biomass_data` data frame.  `biomass_data` 
#' must contain the following columns: `Group`, `Year`, `Value`, `Scale`, `Source`, 
#' one of either `Stdev` or `CV` (or both), and optionally `Type`. If both `Stdev`
#' and `CV` are provided, `Stdev` is used for each data point unless it is missing (NA)
#' in which case `CV` is used.
#'
#' @param Rsim.scenario An `Rsim` scenario object containing simulation parameters and setup.
#' @param biomass_data A `data.frame` containing biomass observation data. 
#'
#' @return Modifies and returns a new Rsim scenario object, attaching the cleaned data frame to 
#'   `Rsim.scenario$fitting$Biomass` with additional calculated columns (`obs`, `sd`, `initial_q`, `wt`).
#'
#'@export
read.fitting.biomass <- function(Rsim.scenario, biomass_data){
  #originally read.fitting.biomass.cv  
  # TODO KYA 7/24/24: add intelligent warning if missing column (e.g. Type)
  
  # Base variables
  SIM   <- Rsim.scenario
  years <- as.numeric(row.names(Rsim.scenario$fishing$ForcedFRate))
  species <- SIM$params$spname
  
  # check if CV column exists to prevent script errors
  if(!"CV" %in% colnames(biomass_data)) biomass_data$CV <- NA
  
  #biomass_data  <-read.csv(filename)
  missing_sp <- unique(biomass_data$Group[which(!(biomass_data$Group %in% species))])
  if(length(missing_sp)>0){
    warning("Following species in BIOMASS fit data not in model, dropped: ",missing_sp)
  }
  
  # drop Groups not found in model
  cmdat <- biomass_data[biomass_data$Group %in% species,]
  # drop lines with NAs in Value or Stdev or Scale or CV
  c0dat <- cmdat[!is.na(cmdat$Value) & !is.na(cmdat$Scale) & 
                   (!is.na(cmdat$Stdev) | !is.na(cmdat$CV)), ]
  # drop lines not in model run years
  c1dat <- c0dat[c0dat$Year %in% years,]
  # drop lines where Value, StDev or Scale are <=0
  ccdat <- c1dat[c1dat$Value>0 & c1dat$Scale!=0,]
  # replace NAs in type with index (maybe with future warning)
  ccdat$Type[is.na(ccdat$Type)] <- "index"
  
  #type <- as.character(rep("absolute",length(ccdat$YEAR)))
  ccdat$Year  <- as.character(ccdat$Year)
  ccdat$Group <- as.character(ccdat$Group)
  ccdat$Source <- as.character(ccdat$Source)
  ccdat$ID    <- paste(as.character(ccdat$Group),as.character(ccdat$Source),sep=":")
  
  # Select/Calculate SD
  # give priotity to Stdev; if NA or 0, use Value * CV
  use_cv <- (is.na(ccdat$Stdev) | ccdat$Stdev <= 0) & (!is.na(ccdat$CV) & ccdat$CV > 0)
  if(any(use_cv)){
    cv_sources <- unique(ccdat$ID[use_cv])
    message("Note: Calculating Stdev from CV for: ", paste(cv_sources, collapse=", "))
  }
  final_stdev <- ifelse(use_cv, 
                        ccdat$Value * ccdat$CV, 
                        ccdat$Stdev)
  
  obs  <- ifelse(as.numeric(ccdat$Scale)<0, as.numeric(ccdat$Value),
                 as.numeric(ccdat$Value) * as.numeric(ccdat$Scale))   
  sd   <- ifelse(as.numeric(ccdat$Scale)<0, as.numeric(final_stdev),
                 as.numeric(final_stdev) * as.numeric(ccdat$Scale)) 
  # adding a safety net to find if any sd are NA or <=0 after calculation
  valid_indices <- !is.na(sd) & sd > 0
  ccdat <- ccdat[valid_indices, ]
  obs   <- obs[valid_indices]
  sd    <- sd[valid_indices]
  
  wt   <- rep(1,length(obs))
  initial_q   <- rep(1,length(obs))
  SIM$fitting$Biomass <- cbind(ccdat,obs,sd,initial_q,wt)
  
  return(SIM)
  
}

################################################################################
#' Read and Process Catch Fitting Data
#'
#' Reads, cleans, and formats catch observation data to prepare it for model 
#' fitting within an `Rpath` simulation scene object.
#'
#' @details
#' The function performs several filtering and data-cleaning steps on `cdat`:
#' * Drops species/groups in the data that are not present in the model.
#' * Removes rows containing `NA` in the `Value` column.
#' * Filters out observations that fall outside the model's defined run years.
#' * Calculates `obs` (observed catch) by multiplying `Value` by `Scale`.
#' * Calculates `sd` (standard deviation) by multiplying `Stdev` by `Scale`.
#'
#' @param Rsim.scenario An `Rsim` scenario object containing simulation parameters and setup.
#' @param catch_data A `data.frame` containing catch observation data. Expected to contain 
#'   the columns: `Group`, `Year`, `Value`, `Stdev`, and `Scale`.
#'
#' @return Modifies and returns the Rsim.scenario, attaching the formatted data frame 
#'   to `Rsim.scenario$fitting$Catch` with additional calculated columns (`obs`, `sd`, `wt`).
#'   
#'@export
read.fitting.catch <- function(Rsim.scenario, catch_data){
  
  # TODO KYA 7/24/24: add intelligent warning if missing column 
  
  SIM   <- Rsim.scenario
  years <- as.numeric(row.names(Rsim.scenario$fishing$ForcedFRate))
  # Columns needed
  #  Group	Year	Value	SD	Scale   
  #catch_data  <- read.csv(filename)
  missing_sp <- unique(catch_data$Group[which(!(catch_data$Group %in% SIM$params$spname))])
  if(length(missing_sp)>0){
    warning("Following species in CATCH fit data not in model, dropped: ",missing_sp)
  }
  # drop Groups not found in model
  cmdat <- catch_data[catch_data$Group %in% SIM$params$spname,]
  # drop lines with NAs in Value or with years outside scenario years range
  ccdat <- cmdat[!is.na(cmdat$Value) & cmdat$Year %in% years,] 
  ccdat$Year  <- as.character(ccdat$Year)
  ccdat$Group <- as.character(ccdat$Group) 
  obs  <- as.numeric(ccdat$Value) * as.numeric(ccdat$Scale)   
  sd   <- as.numeric(ccdat$Stdev) * as.numeric(ccdat$Scale)  
  wt   <- rep(1,length(obs))
  SIM$fitting$Catch <- cbind(ccdat,obs,sd,wt)
  #sdat  <- aggregate(as.numeric(ccdat$Value)*as.numeric(ccdat$Scale),list(ccdat$Year,ccdat$Group),"sum")
  #sd    <- 0.1*sdat$x
  #colnames(SIM$fitting$CATCH) <- c("year","species","obs","sd","wt")

  # Apply fit fishing to matrix
  #SIM$fishing$ForcedEffort[] <- 0
  #SIM$fishing$ForcedCatch[matrix(c(SIM$fitting$Catch$Year, SIM$fitting$Catch$Group),
  #                        length(SIM$fitting$Catch$Year),2)] <- SIM$fitting$Catch$obs
  return(SIM)
}

################################################################################
#' Remove Biomass Timeseries from Fitting Data
#'
#' Removes specific biomass observation timeseries from an `Rsim` scenario 
#' object based on their unique identification strings.
#'
#' @details 
#' This is a utility function used to drop specific observation datasets (e.g., 
#' a specific survey for a specific species) prior to running a model fit. 
#' It filters the `scene$fitting$Biomass` data frame, keeping only the rows 
#' where the `ID` does NOT match the values provided in the `timeseries` argument.
#'
#' @param scene An `Rsim` scenario object containing the fitting data setup.
#' @param timeseries A character vector of `ID` strings representing the timeseries 
#'   to be removed (e.g., `"SpeciesName:SurveySource"`).  Existing timeseries
#'   IDs can be listed using the `rsim.fit.list.bio.series` function.
#'
#' @return Modifies and returns the `scene` object, with the specified timeseries 
#'   removed from the `scene$fitting$Biomass` data frame.
#'@export
rsim.fit.remove.bio.timeseries <- function(Rsim.scenario, timeseries){
  Rsim.scenario$fitting$Biomass <- Rsim.scenario$fitting$Biomass[!(Rsim.scenario$fitting$Biomass$ID %in% timeseries),]
  return(Rsim.scenario)
}

################################################################################
#' Convert Fitting Catch to Forced Catch in Rsim Scenario
#'
#' Performs 3 steps:  Sets effort to 0, sets default F, sets catch.
#' Integrates observed catch data into an `Rpath` simulation scene by zeroing 
#' out effort, setting a baseline fishing mortality rate, and forcing specific 
#' catches for observation years.
#'
#' @details
#' The function modifies the fishing parameters of the simulation by:
#' 1. Setting all gear effort (`ForcedEffort`) to zero.
#' 2. Calculating a baseline Fishing Mortality rate (F) for all living and detrital 
#'    groups based on the balanced `Rpath` model (`(Landings + Discards) / Biomass`).
#' 3. Assigning this calculated F to the `ForcedFRate` matrix for all years.
#' 4. Overriding specific years and groups with empirical catch data (`fitting$Catch$obs`).
#'    For these targeted overrides, the `ForcedFRate` is set to 0 to allow the 
#'    forced catch to drive the simulation.
#'
#' @param Rsim.scenario An `Rsim` scenario object containing the observed fishing data.
#' @param Rpath A balanced `Rpath` model object containing base biomass and landings.
#'
#' @return Modifies and returns the `Rsim.scenario` object with updated 
#'   `ForcedEffort`, `ForcedFRate`, and `ForcedCatch` matrices.
#'
#'@export
fitcatch.to.forcecatch <- function(Rsim.scenario, Rpath){

# TODO: figure out detritus
  scene <- Rsim.scenario
  bal   <- Rpath
  
  # zero out all gear effort
  scene$fishing$ForcedEffort[] <- 0 
  
  # Set forced Frate to Ecopath F
  splist <- c(rpath.living(bal),rpath.detrital(bal))
  flist  <- (rowSums(bal$Landings) + rowSums(bal$Discards))/bal$Biomass
  for (sp in splist){
    scene$fishing$ForcedFRate[,sp] <- flist[sp]
  }
  
  #Put supplied catch in ForcedCatch, and zero Frate for those years
if(!is.null(scene$fitting$Catch)){  
  catchmat <- matrix(c(scene$fitting$Catch$Year, scene$fitting$Catch$Group),
                     length(scene$fitting$Catch$Year),2)
  
  scene$fishing$ForcedCatch[catchmat] <- scene$fitting$Catch$obs
  scene$fishing$ForcedFRate[catchmat] <- 0
}  
  return(scene)

}

#################################################################################
#' List Catch Timeseries Groups
#'
#' Extracts the unique species or functional groups that have observation data 
#' within the catch fitting matrix of an `Rsim` simulation scenario.
#'
#' @details 
#' This is a helper function designed to quickly identify which groups have 
#' empirical catch data loaded into the scenario. It reads the `Group` column 
#' from the `Rsim.scenario$fitting$Catch` data frame and returns the unique values.
#'
#' @param Rsim.scenario An `Rsim` scenario list object containing 
#'   fitting data setup.
#'
#' @return A list containing a single element, `groups`, which holds a vector 
#'   of the unique group names present in the catch data.
#'@export
rsim.fit.list.catch.series <- function(Rsim.scenario){
  
  scene <- Rsim.scenario
  out   <- list()
  out$groups  <- unique(scene$fitting$Catch$Group)
  return(out)  
  
}

#################################################################################
#' List Biomass Timeseries Groups by Survey Source
#'
#' Extracts and groups the unique species or functional groups associated with 
#' each survey source in the biomass fitting data of an `Rsim` simulation scene.
#'
#' @details 
#' This helper function reads the `Source` and `Group` columns from the 
#' `scene$fitting$Biomass` data frame. It iterates through every unique survey 
#' source and compiles a list of all groups (species) that have observation 
#' data attributed to that specific source.
#'
#' @param Rsim.scenario An `Rpath` simulation scene list object containing 
#'   the fitting data setup.
#'
#' @return A named list where each element name corresponds to a unique survey 
#'   `Source`, and the contents are character vectors of the `Group` names 
#'   associated with that source.
#'
#'@export
rsim.fit.list.bio.bysurvey <- function(Rsim.scenario){
  
  scene <- Rsim.scenario
  out   <- list()
  sources <- unique(scene$fitting$Biomass$Source)
  for (i in 1:length(sources)){
   out[[i]] <- unique(scene$fitting$Biomass$Group[scene$fitting$Biomass$Source==sources[i]])
   names(out)[i] <- sources[i]
  }
  
  return(out)  
  
}

#################################################################################
#' List Biomass Timeseries Groups and Sources
#'
#' Extracts all unique groups, survey sources, and combined identifiers from 
#' the biomass fitting data within an `Rsim` simulation scene.
#'
#' @details 
#' This helper function reads the `Group` and `Source` columns from the 
#' `scene$fitting$Biomass` data frame. It returns a summary list containing 
#' the unique functional groups, the unique survey sources, and the unique 
#' combinations of the two formatted as `"Group:Source"`.
#'
#' @param Rsim.scenario An `Rsim` simulation scene list object containing 
#'   the fitting data setup.
#'
#' @return A list containing three elements: 
#'   * `groups`: A character vector of unique species/groups.
#'   * `sources`: A character vector of unique survey sources.
#'   * `all`: A character vector of unique combinations (`"Group:Source"`).
#'
#'@export
rsim.fit.list.bio.series <- function(Rsim.scenario){
  
  scene <- Rsim.scenario
  out   <- list()
  out$groups  <- unique(scene$fitting$Biomass$Group)
  out$sources <- unique(scene$fitting$Biomass$Source)
  out$all     <- unique(paste(scene$fitting$Biomass$Group,scene$fitting$Biomass$Source,sep=':'))
  return(out)  
  
}

#################################################################################
#' Set Catch Fitting Data Weights
#'
#' Updates the statistical weighting (`wt`) for specific catch observation timeseries 
#' prior to running a model fit in an `Rsim` simulation scene. Matches are based 
#' on both the functional group and the survey source.
#' 
#' @details 
#' This function locates all observation rows in the `Rsim.scenario$fitting$Catch` data 
#' frame that match BOTH the specified `group` and `source` names. It replaces 
#' their `"wt"` column with the provided `wt` value. If the specific combination 
#' is not found in the catch data, it issues a warning and returns the scene unmodified.
#'
#' @param Rsim.scenario An Rsim scenario object containing the fitting data.
#' @param group A character string or vector of group names to update.
#' @param source A character string or vector of survey source names to update.
#' @param wt A numeric value to assign as the new weight.
#'
#' @return Modifies and returns the `Rsim.scenario`object with the 
#'   updated weights applied to `Rsim.scenario$fitting$Biomass`.
#'
#'@export
rsim.fit.set.catch.wt <- function(Rsim.scenario, group, wt){
  
  scene  <- Rsim.scenario
  series <- scene$fitting$Catch$Group %in% group  
  
  if (sum(series)==0){
    warning("No catch fitting data found for ",group) 
    return(scene)
  }
  
  scene$fitting$Catch[series,"wt"] <- wt
  
  return(scene)
  
}

#################################################################################
#' Set Biomass Fitting Data Weights
#'
#' Updates the statistical weighting (`wt`) for specific biomass observation timeseries 
#' prior to running a model fit in an `Rsim` simulation scene. Matches are based 
#' on both the functional group and the survey source.
#' 
#' @details 
#' This function locates all observation rows in the `Rsim.scenario$fitting$Biomass` data 
#' frame that match BOTH the specified `group` and `source` names. It replaces 
#' their `"wt"` column with the provided `wt` value. If the specific combination 
#' is not found in the biomass data, it issues a warning and returns the scene unmodified.
#'
#' @param Rsim.scenario An Rsim scenario object containing the fitting data.
#' @param group A character string or vector of group names to update.
#' @param source A character string or vector of survey source names to update.
#' @param wt A numeric value to assign as the new weight.
#'
#' @return Modifies and returns the `Rsim.scenario`object with the 
#'   updated weights applied to `Rsim.scenario$fitting$Biomass`.
#'
#'@export
rsim.fit.set.bio.wt <- function(Rsim.scenario, group, source, wt){

  scene <- Rsim.scenario
  series <- scene$fitting$Biomass$Source %in% source & 
            scene$fitting$Biomass$Group %in% group  

  if (sum(series)==0){
    warning("No biomass fitting data found for ",source," ",group) 
    return(scene)
  }
  
  scene$fitting$Biomass[series,"wt"] <- wt
  
  return(scene)
  
}

#################################################################################
#' Set Catchability (q) for Biomass Fitting Data
#'
#' Updates the catchability coefficient (`initial_q`) and timeseries `Type` for 
#' specific biomass observations within an `Rpath` simulation scene. 
#'
#' @details 
#' This function operates in three distinct modes depending on the arguments provided:
#' 1. **Default Reset (`q` and `years` are `NULL`):** Resets `initial_q` to `1.0` 
#'    and sets `Type` to `"index"`.
#' 2. **Direct Assignment (`q` is provided):** Assigns the provided numeric `q` 
#'    to `initial_q` and assigns the specified `type` (which defaults to `"fixed"`).
#' 3. **Empirical Calculation (`years` is provided, `q` is `NULL`):** Calculates `q` 
#'    by taking the mean observed biomass (`Value * Scale`) over the specified 
#'    `years` and dividing it by the model's base reference biomass (`B_BaseRef`) 
#'    for the specified `group`.
#'
#' @param Rsim.scenario An Rsim scenario object containing the fitting data.
#' @param group A character string (or vector) of the functional group(s) to update.
#' @param source A character string (or vector) of the survey source(s) to update.
#' @param q An optional numeric value to directly set the catchability coefficient.
#' @param years An optional numeric or character vector of years used to calculate `q` empirically.
#' @param type An optional character string to set the timeseries type (defaults to `"fixed"`).
#'
#' @return Modifies and returns the `Rsim.scenario` (scene) object with updated 
#'   `initial_q` and `Type` columns in the `scene$fitting$Biomass` data frame.
#'
#'@export
rsim.fit.set.q <- function(Rsim.scenario, group, source, q=NULL, years=NULL, type=NULL){
  
  scene <- Rsim.scenario
  series <- scene$fitting$Biomass$Source %in% source & 
            scene$fitting$Biomass$Group %in% group
  if(is.null(type)){type="fixed"}
  
  if (sum(series)==0){
    warning("No biomass fitting data found for ",source," ",group) 
    return(scene)
  }
  
  # If both q and years are null, set initial_q to 1.0 and Type to "index"
  if (is.null(q) & is.null(years)){
    scene$fitting$Biomass[series,"initial_q"] <- 1.0
    scene$fitting$Biomass[series,"Type"]      <- "index"
    return(scene)
  }
  
  # If a non-NULL q is supplied, use that
  if (!is.null(q)){
    qq <- as.numeric(q)
    if (!is.na(qq) & qq>0){
      scene$fitting$Biomass[series,"initial_q"] <- qq
      scene$fitting$Biomass[series,"Type"]      <- type
      return(scene)
    } else {
      warning("Supplied q for ",source," ",group," is NA or non-positive - q unchanged.")
      return(scene)
    }
  }
      
  # At this point, use years
  lookup <- scene$fitting$Biomass$Source %in% source & 
            scene$fitting$Biomass$Group %in% group &  
            scene$fitting$Biomass$Year %in% years
  
  dat <- scene$fitting$Biomass[lookup,]
  est <- scene$params$B_BaseRef[group]
  obs <- mean(dat$Value * dat$Scale)
  qq  <- obs/est
  
  if (is.na(qq) || is.nan(qq) || qq<=0){
    warning(source," ",group," survey data in ", years, " is NA, NaN or non-positive - q unchanged.")
  } else {
    scene$fitting$Biomass[series,"initial_q"] <- qq
    scene$fitting$Biomass[series,"Type"]      <- type
  }
  
  return(scene)
}

#################################################################################
#' Calculate Objective Function for Rsim Model Fit
#'
#' Computes the negative log-likelihood (goodness of fit) between an `Rsim` 
#' simulation's dynamic output and empirical observation timeseries data for 
#' both Biomass and Catch. 
#'
#' @details 
#' This function evaluates how well a simulation run (`Rsim.output`) matches the loaded 
#' observation data (`Rsim.scenario$fitting`). It uses a log-normal error distribution 
#' and calculates the negative log-likelihood. 
#' 
#' For biomass observation series marked as `"index"`, the function analytically 
#' calculates a variance-weighted catchability coefficient (`q_est`) to scale 
#' the estimates prior to calculating the log-likelihood error.
#'
#' @param Rsim.scenario An Rsim scenario object containing the fitting data.
#' @param Rsim.output An Rsim simulation output to be evaluated.
#' @param verbose Logical; if `TRUE`, returns a detailed list containing the 
#'   original fitting data appended with model estimates, `q` values, and 
#'   individual fit scores. If `FALSE`, returns only the total sum of the 
#'   negative log-likelihoods (useful for passing to optimizer functions).
#'
#' @return If `verbose = TRUE`, returns a list containing detailed data frames 
#'   (`Biomass`, `Catch`) and the total objective score (`tot`). 
#'   If `verbose = FALSE`, returns a single numeric scalar representing the total 
#'   negative log-likelihood.
#'   
#'@export
rsim.fit.obj <- function(Rsim.scenario, Rsim.output, verbose=TRUE){
  FLOGTWOPI <- 0.5*log(2*pi) #0.918938533204672
  epsilon <- 1e-36
  
  OBJ <- list()
  OBJ$tot <- 0
  
  # BIOMASS to NON-RESCALED "Actual" biomass estimate
  est <- Rsim.output$annual_Biomass[matrix(c(as.character(Rsim.scenario$fitting$Biomass$Year),as.character(Rsim.scenario$fitting$Biomass$Group)),
                                   ncol=2)] + epsilon
  obs <- Rsim.scenario$fitting$Biomass$obs + epsilon
  sd  <- Rsim.scenario$fitting$Biomass$sd  + epsilon
  wt  <- Rsim.scenario$fitting$Biomass$wt
  initial_q <- Rsim.scenario$fitting$Biomass$initial_q
  # Series id (sid) is Source and Group columns combined
  sid <- paste(Rsim.scenario$fitting$Biomass$Source, Rsim.scenario$fitting$Biomass$Group, sep=":")
  
  # We need to get variance-weighted survey means by species, for
  # calculating mean values needed for setting best-fit q
  
  # Formula for weighted average q: 
  # q = exp(sum(w * log(obs/est))/sum(w)) where w is wt/sd  
  logdiff       <- log(obs/est)
  sdlog         <- sqrt(log(1.0+sd*sd/(obs*obs))) # sigma^2 of lognormal dist 
  wt_sd_inverse <- wt/sdlog# sd
  wt_logdiffsum <- tapply(logdiff*wt_sd_inverse, sid ,sum)
  wt_sum        <- tapply(wt_sd_inverse,         sid ,sum)
  q_est         <- exp(wt_logdiffsum/wt_sum) # need ifelse here for 0 weights?
  survey_q      <- ifelse(Rsim.scenario$fitting$Biomass$Type=="index", 
                          q_est[sid], initial_q)
  survey_q      <-ifelse(is.na(survey_q) | is.nan(survey_q),initial_q,survey_q)
  
  ## Jan 2023 incorrect code
  #inv_var <- 1.0/(sd*sd)
  #obs_sum <- tapply(obs*inv_var*wt, as.character(Rsim.scenario$fitting$Biomass$Group),sum)
  #inv_sum <- tapply(inv_var*wt,     as.character(Rsim.scenario$fitting$Biomass$Group),sum)
  #obs_mean <- obs_sum/inv_sum
  #est_mean <- tapply(est,as.character(Rsim.scenario$fitting$Biomass$Group),mean)
  #survey_q <- ifelse(Rsim.scenario$fitting$Biomass$Type=="absolute", 1.0,
  #            #(obs_mean/est_mean)[as.character(Rsim.scenario$fitting$Biomass$Group)])
  #            (est_mean/obs_mean)[as.character(Rsim.scenario$fitting$Biomass$Group)])
  #obs_scaled <-obs*survey_q 
  #sdlog  <- sqrt(log(1.0+sd*sd*survey_q*survey_q/(obs_scaled*obs_scaled)))
  sdiff  <- log((obs/survey_q)/est)/sdlog
  fit    <- wt * (FLOGTWOPI + log(sdlog) + 0.5*sdiff*sdiff)
  
  if (verbose){
    obs_scaled  <- obs/survey_q
    OBJ$Biomass <- cbind(Rsim.scenario$fitting$Biomass,est,survey_q,obs_scaled,sdiff,fit)
  } else {
    OBJ$tot <- OBJ$tot + sum(fit)
  }
  
  # Catch compared (assumes all catch is clean, absolute values)
  est <- Rsim.output$annual_Catch[matrix(c(as.character(Rsim.scenario$fitting$Catch$Year),as.character(Rsim.scenario$fitting$Catch$Group)),
                                 ncol=2)] + epsilon
  obs <- Rsim.scenario$fitting$Catch$obs + epsilon
  sd  <- Rsim.scenario$fitting$Catch$sd  + epsilon
  sdlog  <- sqrt(log(1.0+sd*sd/(obs*obs)))
  sdiff  <- log(obs/est)/sdlog
  fit    <- Rsim.scenario$fitting$Catch$wt * (log(sdlog) + FLOGTWOPI + 0.5*sdiff*sdiff)
  if (verbose){
    OBJ$Catch <- cbind(Rsim.scenario$fitting$Catch,est,sdiff,fit)
  } else {
    OBJ$tot <- OBJ$tot + sum(fit)
  }
  
  # # RATION
  # obs <- Rsim.scenario$fitting$ration$obs + epsilon
  # sd  <- Rsim.scenario$fitting$ration$sd  + epsilon
  # inv_var <- (1.0/sd)*(1.0/sd)
  # obs_sum <- tapply(obs*inv_var,as.character(Rsim.scenario$fitting$ration$Group),sum)
  # inv_sum <- tapply(inv_var,as.character(Rsim.scenario$fitting$ration$Group),sum)
  # obs_mean <- obs_sum/inv_sum
  # est <- Rsim.output$annual_QB[matrix(c(as.character(Rsim.scenario$fitting$ration$Year),as.character(Rsim.scenario$fitting$ration$Group)),
  #                             ncol=2)] + epsilon
  # est_mean <- tapply(est,as.character(Rsim.scenario$fitting$ration$Group),mean)
  # survey_q <- (obs_mean/est_mean)[as.character(Rsim.scenario$fitting$ration$Group)]
  # est_scaled <-est*survey_q 
  # sdlog  <- sqrt(log(1.0+sd*sd/(obs*obs)))
  # sdiff  <- (log(obs)-log(est_scaled))/sdlog
  # fit    <- Rsim.scenario$fitting$ration$wt * (log(sdlog) + FLOGTWOPI + 0.5*sdiff*sdiff)
  # OBJ$ration <- cbind(GOA_Rsim.scenario$fitting$ration,est,survey_q,est_scaled,sdiff,fit)
  # 
  # # Diet proportions estimation
  # linklook   <- matrix(c(as.character(Rsim.scenario$fitting$diets$Year),as.character(Rsim.scenario$fitting$diets$simlink)),ncol=2)
  # totlook    <- matrix(c(as.character(Rsim.scenario$fitting$diets$Year),as.character(Rsim.scenario$fitting$diets$pred)),ncol=2) 
  # dietTot    <- tapply(Rsim.output$annual_Qlink[linklook],list(Rsim.scenario$fitting$diets$Year,Rsim.scenario$fitting$diets$pred),sum)
  # dietProp   <- Rsim.output$annual_Qlink[linklook]/dietTot[totlook]
  # logest     <- log(dietProp)
  # #NEGATIVE log likelihood now
  # fit        <- -Rsim.scenario$fitting$diets$wt * (Rsim.scenario$fitting$diets$log_diff + Rsim.scenario$fitting$diets$alphaM1*logest)  
  # OBJ$diet   <- cbind(Rsim.scenario$fitting$diets,dietProp,logest,fit)
  
  # Final summation and return
  if(verbose){
    OBJ$tot <- sum(OBJ$Biomass$fit, OBJ$Catch$fit)# , OBJ$ration$fit, OBJ$diet$fit)
    return(OBJ)
  }
  else{
    return(OBJ$tot)
  }
}

#################################################################################
#' Extract Summary Table of Model Fits
#'
#' Runs the objective function and aggregates the negative log-likelihood fit scores 
#' by functional group for both Biomass and Catch data, returning them in a 
#' standardized data frame.
#'
#' @details 
#' This function acts as a convenient wrapper around `rsim.fit.obj()`. It runs the 
#' objective function in verbose mode to extract individual observation fit scores, 
#' sums those scores by species/group, and then aligns the results 
#' against the complete species list defined in the model. 
#'
#' @param Rsim.scenario An Rsim scenario object containing the fitting data.
#' @param Rsim.output An Rsim simulation output to be evaluated.
#'
#' @return A `data.frame` where the row names correspond to the model's full species 
#'   list. It contains two columns (`Biomass` and `Catch`) 
#'   representing the summed negative log-likelihood fit scores for each group.
#'
#'@export
rsim.fit.table <- function(Rsim.scenario, Rsim.output){
  fitobj  <- rsim.fit.obj(Rsim.scenario,Rsim.output,verbose=T)
  Btmp <- tapply(fitobj$Biomass$fit,fitobj$Biomass$Group,sum)
  Ctmp <- tapply(fitobj$Catch$fit,fitobj$Catch$Group,sum)
  out <- rep(NA,length(Rsim.scenario$params$spname)); names(out)<- Rsim.scenario$params$spname
  Biomass <- out; Biomass[names(Btmp)] <- Btmp
  Catch <- out;   Catch[names(Ctmp)] <- Ctmp
  return(data.frame(Biomass,Catch))
}

#################################################################################
#' Extract Objective Function Fits by Species
#'
#' Runs the objective function and extracts the detailed, row-by-row fit scores 
#' strictly for specific species or functional groups from both the Biomass 
#' and Catch data.
#'
#' @details 
#' This function acts as a filter on top of `rsim.fit.obj()`. It runs the base 
#' objective function in verbose mode to calculate the individual negative 
#' log-likelihoods for all observation data. It then subsets the resulting 
#' `Biomass` and `Catch` data frames, returning only the rows where the `Group` 
#' matches the provided `species` vector.
#'
#' @param Rsim.scenario An Rsim scenario object containing the fitting data.
#' @param Rsim.output An Rsim simulation output to be evaluated.
#' @param species A character vector of species or group names to extract. Defaults to `NULL`.
#'
#' @return A list containing two data frames (`Biomass` and `Catch`) representing 
#'   the detailed observation-level fit data filtered to the requested species.
#'
#'@export
rsim.fit.obj.species <- function(Rsim.scenario, Rsim.output, species=NULL){
  OBJ <- list()
  fitobj <- rsim.fit.obj(Rsim.scenario,Rsim.output,verbose=T)
  OBJ$Biomass <- fitobj$Biomass[fitobj$Biomass$Group%in%species,] 
  OBJ$Catch   <- fitobj$Catch[fitobj$Catch$Group%in%species,]
  return(OBJ)
}

#################################################################################
#' Get Predator-Prey Link Index
#'
#' Retrieves the numeric index of a specific predator-prey interaction from the 
#' model's internal diet network structure.
#'
#' @details 
#' This helper function scans the internal C++ 0-indexed arrays (`PreyFrom` and 
#' `PreyTo`), converts them to 1-based R indices to look up the corresponding 
#' species names in `spname`, and returns the link index (or indices) where 
#' both the predator and prey match the requested strings.
#'
#' @param Rsim.scenario An `Rpath` simulation scene list object containing model parameters 
#'   (`scene$params$spname`, `scene$params$PreyFrom`, and `scene$params$PreyTo`).
#' @param predator A character string specifying the name of the predator group.
#' @param prey A character string specifying the name of the prey group.
#'
#' @return A numeric scalar (or vector) representing the internal link index 
#'   associated with the specified predator-prey interaction.
#'
#'@export
get.pp.link <- function(Rsim.scenario, predator, prey){
  
  return(as.numeric(which(Rsim.scenario$params$spname[Rsim.scenario$params$PreyFrom + 1] == prey &
               Rsim.scenario$params$spname[Rsim.scenario$params$PreyTo   + 1] == predator
               ))) #+1 here is due to R to C++ array conversion  
  
}

#################################################################################
#' Apply Fitted Parameters to Rpath Simulation Scene
#' Helper function that maps a vector of optimized parameter modifications 
#' back into the `Rpath` simulation parameters object (`scene.params`).  Since
#' this returns a scenario params object instead of a full scenario, it is
#' primarily intended for internal use.
#'
#' @details 
#' This function processes a vector of parameter modifications (`values`) based 
#' on their assigned `vartype`. It handles four parameter types:
#' * `"mzero"`: Modifies non-predation mortality directly.
#' * `"predvul"`: Modifies vulnerability from the predator's perspective.
#' * `"preyvul"`: Modifies vulnerability from the prey's perspective.
#' * `"ppvul"`: Modifies vulnerability for a single, specific predator-prey link.
#' For vulnerabilities, the modifications are applied to the base vulnerability (`VV`) 
#' using a logarithmic transformation: `VV_new = 1 + exp(log(VV_base - 1) + modifications)`.
#' 
#' @param values A numeric vector of parameter values or modifiers.  For 
#' a parameter of vartype `"mzero"`, the unit of the value is a mortality
#' rate.  Values for `"predvul"`, `"preyvul"`, and `"ppvul"` are log-scaled
#' vulnerabilities.   
#' @param species A character vector indicating the species/groups to modify 
#'   (or string representations of numeric link indices for `"ppvul"`).
#' @param vartype A character vector indicating the types of parameters to modify. 
#'   Currently supported parmaters are `"mzero"`, `"predvul"`, `"preyvul"`, and 
#'   `"ppvul"`.
#' @param scene.params An `Rpath` simulation parameters list (usually `scene$params`).
#'
#' @return Modifies and returns the `scene.params` object with updated `MzeroMort` 
#'   and `VV` matrices.
#'
#'@export
rsim.fit.apply <- function(values, species, vartype, scene.params){

# Mzero
  mzerodiff <- values[vartype=="mzero"]
  mzero.sp  <- species[vartype=="mzero"]
  scene.params$MzeroMort[mzero.sp] <- mzerodiff

# Predator contribution to predprey vulnerability   
  predvuls <- values[vartype=="predvul"]
  names(predvuls) <- species[vartype=="predvul"]   
  preddiff <- as.numeric(predvuls[scene.params$spname[scene.params$PreyTo+1]])
  preddiff[is.na(preddiff)] <- 0
  
# Prey contribution to predprey vulnerability
  preyvuls <- values[vartype=="preyvul"]
  names(preyvuls) <- species[vartype=="preyvul"]   
  preydiff <- as.numeric(preyvuls[scene.params$spname[scene.params$PreyFrom+1]])
  preydiff[is.na(preydiff)] <- 0
  
# Single link vulnerability, indexed by link number
  ppvuls   <- values[vartype=="ppvul"]
  pp_index <- as.numeric(species[vartype=="ppvul"])
  ppdiff <- rep(0, scene.params$NumPredPreyLinks+1) # +1 for rccp offset #BIA ####
  ppdiff[pp_index] <- ppvuls
  ppdiff[is.na(ppdiff)] <- 0
  
# Apply pred and prey vuls above to actual pred/prey VV in scenario
  scene.params$VV <- (1 + exp(log(scene.params$VV-1) + preddiff + preydiff + ppdiff))
  
  return(scene.params)
}

#################################################################################
#' Run Fitted Rpath Simulation and Evaluate Objective Function
#'
#' Applies parameter modifications to an `Rpath` scene, runs the simulation, and 
#' either returns the full simulation output or evaluates the fit against 
#' observation data using a penalized negative log-likelihood.
#'
#' @details
#' This function translates the defined fitting variables into the appropriate
#' variables of an Rsim.scenario object and runs an rsim simulation with the
#' resulting scenario. Applying an `"mzero"`
#' parameter to a species replaces the Rsim.scenario$params$MzeroMort paramter
#' for that species. `"predvul"`, `"preyvul"`, and `"ppvul"` are applied using
#' the formula for each predator/prey link \code{Rsim.scenario$params$VV(predator, prey) =   
#' (1 + exp(log(Rsim.scenario$params$VV(predator, prey)-1) + preddiff + preydiff + ppdiff))}
#' 
#' @param values A numeric vector of parameter values or modifiers.  For 
#' a parameter of vartype `"mzero"`, the unit of the value is a mortality
#' rate.  Values for `"predvul"`, `"preyvul"`, and `"ppvul"` are log-scaled
#' vulnerabilities.   
#' @param species A character vector indicating the species/groups to modify 
#'   (or string representations of numeric link indices for `"ppvul"`).
#' @param vartype A character vector indicating the types of parameters to modify. 
#'   Currently supported parmaters are `"mzero"`, `"predvul"`, `"preyvul"`, and 
#'   `"ppvul"`.
#' @param Rsim.scenario An Rsim.scenario list object.
#' @param run_method Numerical integration method for rsim.run. Either \emph{'AB'} 
#' for Adams-Bashforth or \emph{'RK4'} for 4th order Runge-Kutta.
#' @param run_years A numeric vector of years over which to run the simulation.
#' @param verbose Logical; if `FALSE`, returns the NLL fit score. If `TRUE`, returns 
#'   the complete simulation run object. Defaults to `F`.
#' @param penalty_weight Numeric value to scale the L2 penalty term for parameter 
#'   regularization. Defaults to `0`.
#' @param nll_details Logical; if `TRUE` and `verbose = FALSE`, returns a list 
#'   breaking down the total NLL into base error and penalty components.
#' @param ... Additional arguments (currently swallowed and unused by the function).
#'
#' @return 
#' * If `verbose = TRUE`: Returns an `Rsim` simulation run output object.
#' * If `verbose = FALSE` and `nll_details = FALSE`: Returns a single numeric value 
#'   representing the total penalized negative log-likelihood.
#' * If `verbose = FALSE` and `nll_details = TRUE`: Returns a named list containing 
#'   `total_nll`, `base_error`, and `penalty_term`.
#'
#'@export
rsim.fit.run <- function(values,
                         species,
                         vartype,
                         Rsim.scenario,
                         run_method,
                         run_years,
                         verbose = F,
                         penalty_weight = 0,
                         nll_details=FALSE,
                         ...) {
  if (!all(is.na(values))){
    Rsim.scenario$params <- rsim.fit.apply(values, species, vartype, Rsim.scenario$params)
  }
  run.out <- rsim.run(Rsim.scenario, method = run_method, years = run_years)
  if (!verbose) { #return(rsim.fit.obj(Rsim.scenario, run.out, FALSE))
    base_error <- rsim.fit.obj(Rsim.scenario, run.out, FALSE)
    penalty_term <- 0
    total_nll <- base_error
    if(!all(is.na(values))){
      penalty_term <- penalty_weight * sum(values^2)
      total_nll <- base_error + penalty_term
      
      if(nll_details){
        return(list(
          total_nll = total_nll,
          base_error = base_error,
          penalty_term = penalty_term
        ))
      }
      
    }
    return(total_nll)
  }
  else{
    return(run.out)
  }
} 

#################################################################################
#' Update Rpath Scene with Fitted Parameters
#'
#' A function that applies fit variables to an Rsim scenario, returning a
#' new scenario object.
#' 
#' @details 
#' This function translates the defined fitting variables into the appropriate
#' variables of an Rsim.scenario object.  In particular, applying an `"mzero"`
#' parameter to a species replaces the Rsim.scenario$params$MzeroMort paramter
#' for that species. `"predvul"`, `"preyvul"`, and `"ppvul"` are applied using
#' the formula for each predator/prey link \code{Rsim.scenario$params$VV(predator, prey) =   
#' (1 + exp(log(Rsim.scenario$params$VV(predator, prey)-1) + preddiff + preydiff + ppdiff))}
#' 
#' @param values A numeric vector of parameter values or modifiers.  For 
#' a parameter of vartype `"mzero"`, the unit of the value is a mortality
#' rate.  Values for `"predvul"`, `"preyvul"`, and `"ppvul"` are log-scaled
#' vulnerabilities.   
#' @param species A character vector indicating the species/groups to modify 
#'   (or string representations of numeric link indices for `"ppvul"`).
#' @param vartype A character vector indicating the types of parameters to modify. 
#'   Currently supported parmaters are `"mzero"`, `"predvul"`, `"preyvul"`, and 
#'   `"ppvul"`.
#'
#' @return Modifies and returns the complete `Rsim.scenario` object with updated 
#'   vulnerability and mortality parameters.
#'         
#'@export
rsim.fit.update <- function(values, species, vartype, Rsim.scenario){
  Rsim.scenario$params <- rsim.fit.apply(values, species, vartype, Rsim.scenario$params) 
  return(Rsim.scenario)
}

#################################################################################
#' Extract Predator-Prey Interaction Table
#'
#' Generates a data frame of predator-prey linkages and their 
#' associated foraging parameters from an `Rpath` simulation scenario. 
#
#' @details 
#' This function outputs a data frame of predator/prey link parameters contained 
#' in an Rsim scenario for a given set of predators and prey.
#' 
#'@param Rsim.scenario Scenario object that contains all of the rsim rates and 
#'    forcing functions generated by \code{\link{rsim.scenario}()}.
#' @param pred An optional character string or vector of predator names. If not
#' provided, the function will include all model predators in its output.
#' @param prey An optional character string or vector of prey names. If not
#' provided, the function will include all model prey in its output.
#' 
#' @return A `data.frame` containing the columns: 
#' \itemize{
#'   \item{\code{from}, prey name}, 
#'   \item{\code{to} predator name},
#'   \item{\code{QQ} base consumption rate of the predator/prey link},
#'   \item{\code{VV} vulnerability parameter of the predator/prey link},
#'   \item{\code{DD} handling time parameter of the predator/prey link}
#' }  
#'@export
rsim.predprey.table <- function(Rsim.scenario, pred=NULL, prey=NULL){
  pp_all <- data.frame(
    from = Rsim.scenario$params$spname[Rsim.scenario$params$PreyFrom+1],
    to   = Rsim.scenario$params$spname[Rsim.scenario$params$PreyTo+1],
    QQ   = Rsim.scenario$params$QQ,
    VV   = Rsim.scenario$params$VV,
    DD   = Rsim.scenario$params$DD)
  if (!is.null(pred)){pp_all <- pp_all[pp_all$to   %in% pred,]}
  if (!is.null(prey)){pp_all <- pp_all[pp_all$from %in% prey,]}
  return(pp_all)
}

