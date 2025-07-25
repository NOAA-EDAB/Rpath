% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecosim.R
\name{rsim.params}
\alias{rsim.params}
\title{Initial set up for Rsim dynamic simulation}
\usage{
rsim.params(
  Rpath,
  mscramble = 2,
  mhandle = 1000,
  preyswitch = 1,
  scrambleselfwt = 0,
  handleselfwt = 0,
  steps_yr = 12,
  steps_m = 1
)
}
\arguments{
\item{Rpath}{R object containing a static \code{Rpath} model.}

\item{mscramble}{Base value for vulnerability in functional response; default = 2 (mixed response).}

\item{mhandle}{Base value for handling time in functional response; default = 1000 (off).}

\item{preyswitch}{Exponent for prey switching in functional response; default = 1 (off).}

\item{scrambleselfwt}{Value of 1 indicates all predators overlap in the foraging arena while 0 treats predators individually.}

\item{handleselfwt}{Value of 1 indicates all prey overlap in the arena and contribute to saturation while 0 treats prey individually.}

\item{steps_yr}{Number of time steps per year; default = 12.}

\item{steps_m}{Number of time steps per month; default = 1.}
}
\value{
Returns an object of class \code{Rsim.params}, a list of 39 objects that is passed to the \code{rsim.run} 
   function via the \code{rsim.scenario} function.
\itemize{
 \item{\code{NUM_GROUPS}, number of total model groups}
 \item{\code{NUM_LIVING}, number of living model groups}
 \item{\code{NUM_DEAD}, number of detritus model groups}
 \item{\code{NUM_GEARS}, number of fishery model groups}
 \item{\code{NUM_BIO}, number of living + detritus model groups}
 \item{\code{spname}, named \code{NUM_GROUPS}+1 length character vector beginning with "Outside" containing names of all model groups}
 \item{\code{spnum}, named \code{NUM_GROUPS}+1 length character vector beginning with "Outside" numbered 0, containing numbers of all model groups}
 \item{\code{B_BaseRef}, named \code{NUM_GROUPS}+1 length character vector beginning with "Outside" containing Rpath base biomass of all model groups}
 \item{\code{MzeroMort}, named \code{NUM_GROUPS}+1 length character vector beginning with "Outside" containing unaccounted mortality, calculated as PB * (1.0 - EE), of all model groups}
 \item{\code{UnassimRespFrac}, named \code{NUM_GROUPS}+1 length character vector beginning with "Outside" containing the proportion of consumption that goes to detritus of all model groups}
 \item{\code{ActiveRespFrac}, named \code{NUM_GROUPS}+1 length character vector beginning with "Outside" containing the proportion of consumption that is "lost to heat" for all model groups}
 \item{\code{FtimeAdj}, named \code{NUM_GROUPS}+1 length character vector beginning with "Outside" containing rate of change of feeding time, currently set to 0 for all model groups} 
 \item{\code{FtimeQBOpt}, named \code{NUM_GROUPS}+1 length character vector beginning with "Outside" containing base QB for all consumer model groups, or base PB for primary producers}
 \item{\code{PBopt}, named \code{NUM_GROUPS}+1 length character vector beginning with "Outside" containing base PB for all model groups}
 \item{\code{NoIntegrate}, named \code{NUM_GROUPS}+1 length character vector beginning with "Outside" containing flag set to 0 for high turnover model groups and set to \code{spnum} for all others}
 \item{\code{Handleself}, named \code{NUM_GROUPS}+1 length character vector beginning with "Outside" containing flag for handling time influence, set for all model groups with function argument \code{handleselfwt}, default 0 for individual prey handling time}
 \item{\code{Scrambleself}, named \code{NUM_GROUPS}+1 length character vector beginning with "Outside" containing flag for predator density influence, set for all model groups with function argument \code{scrambleselfwt}, default 0 for individual predator density dependent predation}
 \item{\code{PreyFrom}, numeric vector length \code{NumPredPreyLinks}+1, spnum of prey for each predator prey interaction pair in the model}
 \item{\code{PreyTo}, numeric vector length \code{NumPredPreyLinks}+1, spnum of predator for each predator prey interaction pair in the model}
 \item{\code{QQ}, numeric vector length \code{NumPredPreyLinks}+1, base consumption rate for each predator prey interaction pair in the model}
 \item{\code{DD}, numeric vector length \code{NumPredPreyLinks}+1, handling time effect on functional response for each predator prey pair, set for all model groups with function argument \code{mhandle}, default = 1000}
 \item{\code{VV}, numeric vector length \code{NumPredPreyLinks}+1, vulnerability effect on functional response for each predator prey pair, set for all model groups with function argument \code{mscramble}, default = 2}
 \item{\code{HandleSwitch}, numeric vector length \code{NumPredPreyLinks}+1, prey density dependence effect on functional response for each predator prey pair, set for all model groups with function argument \code{preyswitch}, default = 1}
 \item{\code{PredPreyWeight}, numeric vector length \code{NumPredPreyLinks}+1, relative weight of individual predator to total predators for each predator prey pair, used if \code{scrambleselfwt}>0} 
 \item{\code{PreyPreyWeight}, numeric vector length \code{NumPredPreyLinks}+1, relative weight of individual prey to total prey for each predator prey pair, used if \code{handleselfwt}>0}
 \item{\code{NumPredPreyLinks}, number of predator to prey linkages over all groups in the model}
 \item{\code{FishFrom}, numeric vector length \code{NumFishingLinks}+1, spnum of landing and discard for each fishery interaction in the model}
 \item{\code{FishThrough}, numeric vector length \code{NumFishingLinks}+1, spnum of gear type for each fishery interaction in the model}
 \item{\code{FishQ}, numeric vector length \code{NumFishingLinks}+1, landings or discards relative to base fished group biomass for each fishery interaction in the model}
 \item{\code{FishTo}, numeric vector length \code{NumFishingLinks}+1, spnum of sink for each fishery interaction in the model ("Outside" for landings or detritus group for discards)}
 \item{\code{NumFishingLinks}, number of model group landings and discards to fishery links over all groups in the model}
 \item{\code{DetFrac}, numeric vector length \code{NumDetLinks}+1, fraction of detritus going to DetTo for each living and detritus group in the model}
 \item{\code{DetFrom}, numeric vector length \code{NumDetLinks}+1, spnum flowing to detritus for each living and detritus group in the model}    
 \item{\code{DetTo}, numeric vector length \code{NumDetLinks}+1, spnum of detritus sink for each living and detritus group in the model}
 \item{\code{NumDetLinks}, number of model group links to detritus over all groups in the model}
 \item{\code{BURN_YEARS}, number of model run burn-in (spin up) years, default value -1}
 \item{\code{COUPLED}, number to control species interactions, value of 0 allows density dependent non-interacting species, default value 1}
 \item{\code{RK4_STEPS}, number of RK4 integration steps per month, default value 4}
 \item{\code{SENSE_LIMIT}, numeric vector of multipliers on biomass determining acceptable range for continuing a model run within \code{BURN_YEARS} }
}
}
\description{
Converts the outputs from Rpath into rates for use in Rsim.
}
\examples{
# Read in Rpath parameter file and generate model object
Rpath <- rpath(AB.params)
# Create default dynamic parameters from Rpath model
Rsim.params <- rsim.params(Rpath) 

}
\seealso{
Other Rsim functions: 
\code{\link{extract.node}()},
\code{\link{rsim.plot}()},
\code{\link{rsim.run}()},
\code{\link{rsim.scenario}()},
\code{\link{rsim.step}()},
\code{\link{write.Rsim}()}
}
\concept{Rsim functions}
