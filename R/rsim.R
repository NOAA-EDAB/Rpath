# R version of ecosim 
# originally developed by Kerim Aydin
# modified by Sean Lucey

#'Rsim module of \code{Rpath}
#'
#'Uses a balanced \code{Rpath} model and creates a scenario consisting of 5 list objects:
#'\code{params}, \code{start_state}, \code{forcing}, \code{fishing}, and 
#'\code{stanzas}.
#'
#'@family Rsim functions
#'
#'@inheritParams rpath
#'
#'@param Rpath R object containing a balanced \code{Rpath} model.
#'@param years A vector of each year of the simulation.
#'
#'@return Returns an \code{Rsim.scenario} object that can be supplied to the 
#'     \code{rsim.run} function.
#'@import data.table
#'@useDynLib Rpath
#'@importFrom Rcpp sourceCpp
#'@export
rsim.scenario <- function(Rpath, Rpath.params, years = 1:100){
  # KYA 11/1/17 modifying so years can take a vector of actual years
  # for row labels (for fitting and general printing out)
  if (length(years)<2){stop("Years needs to be a vector of numeric year labels.")}
  
  params      <- rsim.params(Rpath)
  start_state <- rsim.state(params)
  forcing     <- rsim.forcing(params, years)  
  fishing     <- rsim.fishing(params, years)
  stanzas     <- rsim.stanzas(Rpath.params, start_state, params)
  
  # Copy stanza base state to start_state
  start_state$SpawnBio   <-stanzas$baseSpawnBio
  start_state$StanzaPred <-stanzas$baseStanzaPred
  start_state$EggsStanza <-stanzas$baseEggsStanza
  start_state$NageS      <-stanzas$baseNageS
  start_state$WageS      <-stanzas$baseWageS
  start_state$QageS        <-stanzas$baseQageS
    
  #Set NoIntegrate Flags
  ieco <- as.vector(stanzas$EcopathCode[which(!is.na(stanzas$EcopathCode))])
  params$NoIntegrate[ieco + 1] <- -1 * ieco 
  
  rsim = list(params      = params, 
              start_state = start_state,
              forcing     = forcing,
              fishing     = fishing,
              stanzas     = stanzas)
  
  class(rsim) <- 'Rsim.scenario'
  attr(rsim, 'eco.name') <- attr(Rpath, 'eco.name')
  attr(rsim, 'Start year') <- years[1]
  
  return(rsim)   
}

 
