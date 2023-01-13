#'Read Rpath parameters from .csv files
#'
#'Creates an Rpath.param object from a series of .csv files.
#'
#'@family Rpath functions
#'
#'@param modfile file location of the flat file containing the model parameters.
#'@param dietfile file location of the flat file containing the diet parameters.
#'@param stanzagroupfile file location of the flat file containing the group parameters
#'  for multistanza groups.  If not specified a blank stanza list will be created.  
#'@param stanzafile file location of the flat file containing the individual stanza 
#'  parameters for multistanza groups.  If not specified a blank stanza list will 
#'  be created.
#'@param pedfile file location of the flat file containg the pedgigree parameters.
#'@return Outputs an Rpath.param object that can be used for Rpath and subsequently
#'  Rsim.  (NOTE: This does function does not ensure data is correct or in the 
#'  correct locations...run check.rpath.param to ensure the appropriate columns are
#'  present).
#'@export
read.rpath.params <- function(modfile, dietfile, pedfile = NA,
                              stanzagroupfile = NA, stanzafile = NA){
  #Need to define variables to eliminate check() note about no visible binding
  Type <- Group <- V1 <- NULL
  
  Rpath.params <- list()
  Rpath.params$model <- as.data.table(read.csv(modfile,  header = T))
  Rpath.params$diet  <- as.data.table(read.csv(dietfile, header = T))
  
  if(!is.na(stanzagroupfile)){
    stanzagroup <- as.data.table(read.csv(stanzagroupfile, header = T))
    Rpath.params$stanzas$NStanzaGroups <- nrow(stanzagroup)
    Rpath.params$stanzas$stgroups      <- stanzagroup
    Rpath.params$stanzas$stindiv       <- as.data.table(read.csv(stanzafile, 
                                                                 header = T))
  } else {
    Rpath.params$stanzas$NStanzaGroups <- 0
    Rpath.params$stanzas$stgroups      <- data.table(StGroupNum  = NA,
                                                     StanzaGroup = NA,
                                                     nstanzas    = NA,
                                                     VBGF_Ksp    = NA,
                                                     VBGF_d      = NA,
                                                     Wmat        = NA,
                                                     RecPower    = NA,
                                                     Wmat001     = NA,
                                                     Wmat50      = NA,
                                                     Amat001     = NA,
                                                     Amat50      = NA)
    Rpath.params$stanzas$stindiv       <- data.table(StGroupNum   = NA,
                                                     StanzaNum   = NA,
                                                     GroupNum    = NA,
                                                     Group       = NA,
                                                     First       = NA,
                                                     Last        = NA,
                                                     Z           = NA,
                                                     Leading     = NA)
  }
  if(!is.na(pedfile)){
    Rpath.params$pedigree <- as.data.table(read.csv(pedfile, header = T))
  } else {
    Rpath.params$pedigree <- data.table(Group = Rpath.params$model$Group,
                                        B     = 1,
                                        PB    = 1,
                                        QB    = 1,
                                        Diet  = 1)
    fleets <- as.character(Rpath.params$model[Type == 3, Group])
    for(i in 1:length(fleets)){
      Rpath.params$pedigree[, V1 := 1]
      setnames(Rpath.params$pedigree, 'V1', fleets[i])
    }
  }
  class(Rpath.params) <- 'Rpath.params'
  return(Rpath.params)
}