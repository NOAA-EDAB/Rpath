#' Rpath parameter set for the fictitious Anchovy Bay
#'
#' This is a parameter set for the example fictitious ecosystem of Anchovy Bay.
#' Anchovy Bay is used as an introduction to Ecopath with Ecosim.
#' 
#' @format An Rpath.params object
#' @family rpathParameters
#'
#' \describe{
#'\item{model}{Main model parameters including Biomass, PB, QB, EE, etc.}
#'\item{diet}{Diet composition matrix with predators as columns and prey as rows}
#'\item{stanzas}{List object containing necessary parameters for multistanza groups.
#'      Note that Anchovy Bay does not have any multistanza groups.}
#'\item{pedigree}{Matrix of certainty surrounding parameter values.  Defaulted to 1
#'      for this fictitious system.}
#'
#' }
#'
#'
"AB.params"