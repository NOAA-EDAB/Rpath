library(methods)
 
########################################################################################
# Set of functions for returning functional group names (character vector)
# using the type input column
#
# Internal supporting function for group names to check type 
#'
#' @importFrom methods is
#'
grouptype <- function(Rpath) {
  if (is(Rpath,"Rpath")) {
    gt<-list(type=Rpath$type, grp=Rpath$Group)
  } else {
    if (is(Rpath,"Rpath.params")) {
      gt<-list(type=Rpath$model$Type, grp=Rpath$model$Group)}
    else{
      stop("Input must be an Rpath (balanced) or Rpath.params (unbalanced) object.")
    }
  }
  names(gt$type)<-NULL; names(gt$grp)<-NULL
  return(gt)
}
####################################
#' Rpath functional group names
#' 
#' Get a character vector of functional group names from an Rpath object (balanced model)  
#' or an Rpath.params object (unbalanced model parameters) based on the 'type' input
#' parameter as follows: (0: consumers, 1: producers, 2: detrital, 3: gears, 0<type<1: mixotrophs).
#' Living groups are consumers + producers.  Note that mixotrophs are not returned as
#' either consumers or producers, only separately.  
#' 
#'@name rpath.groups
#' 
#'@family Rpath functions
#'
#'@param Rpath Balanced Rpath model generated by rpath.
#'
#'@return Returns a character vector containing the names of Rpath functional groups
#'by category (group type).
#'
#'@export
rpath.groups <- function(Rpath){
  gt <- grouptype(Rpath)
  return(gt$grp)
}

#'@rdname rpath.groups
#'@export
rpath.living <- function(Rpath){
  gt <- grouptype(Rpath)
  return(gt$grp[gt$type<2])
}

#'@rdname rpath.groups
#'@export 
rpath.detrital <- function(Rpath){
  gt <- grouptype(Rpath)
  return(gt$grp[gt$type==2])
}

#'@rdname rpath.groups
#'@export 
rpath.gears     <- function(Rpath){
  gt <- grouptype(Rpath)
  return(gt$grp[gt$type==3])
}

#'@rdname rpath.groups
#'@export 
rpath.producers <- function(Rpath){
  gt <- grouptype(Rpath)
  return(gt$grp[gt$type==1])
}

#'@rdname rpath.groups
#'@export 
rpath.consumers <- function(Rpath){
  gt <- grouptype(Rpath)
  return(gt$grp[gt$type==0])
}

#'@rdname rpath.groups
#'@export 
rpath.mixotrophs <- function(Rpath){
  gt <- grouptype(Rpath)
  return(gt$grp[gt$type>0 & gt$type<0])
}



