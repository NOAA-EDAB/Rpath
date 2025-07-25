########################################################################################
# Set of functions for returning functional group names (character vector)
# using the type input column

#' List of Rpath type and group elements
#' 
#' Internal supporting function for group names to check type.
#'
#' @importFrom methods is
#' 
#' @keywords internal
#'
#' @param Rpath An Rpath (balanced) or Rpath.params (unbalanced) object generated by \code{rpath()}
#' 
#' @returns Returns a list of 2 Rpath elements: the type element and the Group element
#' \item{type}{Vector of group type integers}
#' \item{group}{Vector of character string group names}
#'
#'#' @examples
#' # Build the balanced Rpath model and parameter file by calling `rpath`
#' Rpath <- rpath(Ecosense.EBS)
#' # Run the command
#' gt <- grouptype(Rpath)
#' # Print out the first few columns of the gt object
#' head(gt)
#'
#' @noRd
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

#' Rpath functional group names
#' 
#' Get a character vector of functional group names from an Rpath object (balanced model)  
#' or an Rpath.params object (unbalanced model parameters) based on the 'type' input
#' parameter as follows: (0: consumers, 1: producers, 2: detrital, 3: gears, 0<type<1: mixotrophs).
#' Living groups are consumers + producers.  Note that mixotrophs are not returned as
#' either consumers or producers, only separately.  
#'
#' @name rpath.groups
#' 
#' @family Rpath functions
#'
#' @param Rpath Balanced Rpath model generated by \code{rpath()}
#'
#' @returns Returns a string vector containing the names of Rpath functional groups
#'by category (group type).
#'
#' @export
#'
#' @examples
#' # Build the balanced Rpath model and parameter file by calling `rpath`
#' Rpath <- rpath(Ecosense.EBS)
#' # Run the groups command
#' groups = rpath.groups(Rpath)
#' # Print out the first few group names
#' head(groups)
#'
rpath.groups <- function(Rpath){
  gt <- grouptype(Rpath)
  return(gt$grp)
}

#' 
#' List of Rpath living groups
#' 
#' List of living groups from Rpath object. Living groups are those species groups with type < 2.
#' 
#' @family Rpath functions
#'
#' @param Rpath Balanced Rpath model generated by \code{rpath()}
#'
#' @returns Returns a string vector containing the names of Rpath functional groups
#'that are living.
#'
#' @export
#'
#' @examples
#' # Build the balanced Rpath model and parameter file by calling `rpath`
#' Rpath <- rpath(Ecosense.EBS)
#' # Run the rpath.living command
#' livingGroups = rpath.living(Rpath)
#' # Print out the first few living group names
#' head(livingGroups)
#'
rpath.living <- function(Rpath){
  gt <- grouptype(Rpath)
  return(gt$grp[gt$type<2])
}

#' List of Rpath detrital groups
#'
#' List of detrital groups from Rpath object with species type = 2
#' 
#' @family Rpath functions
#'
#' @param Rpath Balanced Rpath model generated by \code{rpath()}
#'
#' @returns Returns a string vector containing the names of Rpath functional groups
#'that are detrital.
#'
#' @export
#'
#' @examples
#' # Build the balanced Rpath model and parameter file by calling `rpath`
#' Rpath <- rpath(Ecosense.EBS)
#' # Run the rpath.detrital command
#' detritalGroups = rpath.detrital(Rpath)
#' # Print out the first few detrital group names
#' head(detritalGroups)
#'
rpath.detrital <- function(Rpath){
  gt <- grouptype(Rpath)
  return(gt$grp[gt$type==2])
}

#' List of Rpath gears groups
#'
#' List of gears groups from Rpath object with species type = 3
#' 
#' @family Rpath functions
#'
#' @param Rpath Balanced Rpath model generated by \code{rpath()}
#'
#' @returns Returns a string vector containing the names of Rpath functional groups
#'that are gear types
#'
#' @export
#'
#' @examples
#' # Build the balanced Rpath model and parameter file by calling `rpath`
#' Rpath <- rpath(Ecosense.EBS)
#' # Run the groups command
#' gearGroups = rpath.gears(Rpath)
#' # Print out the first few gear type group names
#' head(gearGroups)
#'
rpath.gears <- function(Rpath){
  gt <- grouptype(Rpath)
  return(gt$grp[gt$type==3])
}

#' List of Rpath producer groups
#'
#' List of producer groups from Rpath object with species type = 1
#' 
#' @family Rpath functions
#'
#' @param Rpath Balanced Rpath model generated by \code{rpath()}
#'
#' @returns Returns a string vector containing the names of Rpath functional groups
#'that are producers
#'
#' @export
#'
#' @examples
#' # Build the balanced Rpath model and parameter file by calling `rpath`
#' Rpath <- rpath(Ecosense.EBS)
#' # Run the groups command
#' producerGroups = rpath.producers(Rpath)
#' # Print out the first few producer group names
#' head(producerGroups)
#'
rpath.producers <- function(Rpath){
  gt <- grouptype(Rpath)
  return(gt$grp[gt$type==1])
}

#' List of Rpath consumer groups
#'
#' List of consumer groups from Rpath object with species type = 0.
#' 
#' @family Rpath functions
#'
#' @param Rpath Balanced Rpath model generated by \code{rpath()}
#'
#' @returns Returns a string vector containing the names of Rpath functional groups
#'that are consumers
#'
#' @export
#'
#' @examples
#' # Build the balanced Rpath model and parameter file by calling `rpath`
#' Rpath <- rpath(Ecosense.EBS)
#' # Run the consumers command
#' consumerGroups = rpath.consumers(Rpath)
#' # Print out the first few consumer group names
#' head(consumerGroups)
#' 
rpath.consumers <- function(Rpath){
  gt <- grouptype(Rpath)
  return(gt$grp[gt$type==0])
}

#' List of Rpath mixotroph groups
#'
#' List of groups from Rpath object with species type between 0 and 1
#' 
#' @family Rpath functions
#'
#' @param Rpath Balanced Rpath model generated by \code{rpath()}
#'
#' @returns Returns a string vector containing the names of Rpath functional groups
#'that are mixotrophs
#'
#' @export
#'
#' @examples
#' # Build the balanced Rpath model and parameter file by calling `rpath`
#' Rpath <- rpath(Ecosense.EBS)
#' # Run the groups command
#' mixotrophGroups = rpath.mixotrophs(Rpath)
#' # Print out the first few mixotroph group names
#' head(mixotrophGroups)
#'
rpath.mixotrophs <- function(Rpath){
  gt <- grouptype(Rpath)
  return(gt$grp[gt$type>0 & gt$type<1])
}



