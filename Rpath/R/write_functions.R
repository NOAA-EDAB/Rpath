#'Write function for Ecopath object
#'
#'Outputs basic parameters or mortalities to a .csv file.
#'
#'@family Rpath functions
#'
#'@param x Rpath model created by the ecopath() function.
#'@param file file name for resultant .csv file.  Be sure to include ".csv".
#'@param morts Logical value whether to output basic parameters or mortalities.  
#'
#'@return Writes a .csv file with the basic parameters or mortalities from an Rpath object.
#'@export
#Write -- note, not a generic function
write.Rpath <- function(x, file, morts = F, ...){
  if(morts == F){
    removals <- rowSums(x$Catch) + rowSums(x$Discard)
    out <- data.frame(Group    = x$Group,
                      type     = x$type,
                      TL       = x$TL,
                      Biomass  = x$BB,
                      PB       = x$PB,
                      QB       = x$QB,
                      EE       = x$EE,
                      GE       = x$GE,
                      Removals = removals)
  }
  if(morts == T){
    ngroup <- x$NUM_LIVING + x$NUM_DEAD
    out <- data.frame(Group    = x$Group[1:ngroup],
                      type     = x$type [1:ngroup],
                      PB       = x$PB   [1:ngroup])
    #Calculate M0
    M0  <- c(x$PB[1:x$NUM_LIVING] * (1 - x$EE[1:x$NUM_LIVING]), 
             x$EE[(x$NUM_LIVING + 1):ngroup])
    out <- cbind(out, M0)
    #Calculate F mortality
    totcatch <- x$Catch + x$Discards
    Fmort    <- as.data.frame(totcatch / x$BB[row(as.matrix(totcatch))])
    setnames(Fmort, paste('V',  1:x$NUM_GEARS,                     sep = ''), 
             paste('F.', x$Group[(ngroup +1):x$NUM_GROUPS], sep = ''))
    out  <- cbind(out, Fmort[1:ngroup, ])
    #Calculate M2
    bio  <- x$BB[1:x$NUM_LIVING]
    BQB  <- bio * x$QB[1:x$NUM_LIVING]
    diet <- as.data.frame(x$DC)
    nodetrdiet <- diet[1:x$NUM_LIVING, ]
    detrdiet   <- diet[(x$NUM_LIVING +1):ngroup, ]
    newcons    <- nodetrdiet * BQB[col(as.matrix(nodetrdiet))]
    predM      <- newcons / bio[row(as.matrix(newcons))]
    detcons    <- detrdiet * BQB[col(as.matrix(detrdiet))]
    predM      <- rbind(predM, detcons)
    setnames(predM, paste('V',  1:x$NUM_LIVING,    sep = ''), 
             paste('M2.', x$Group[1:x$NUM_LIVING], sep = ''))
    out <- cbind(out, predM)
  }
  write.csv(out, file = file)
}


#'Write function for Ecosim object
#'
#'Outputs start/end biomass and catch to a .csv file.
#'
#'@family Rpath functions
#'
#'@param Rsim.output object created by the rsim.run() function.
#'@param file file name for resultant .csv file.  Be sure to include ".csv".  
#'
#'@return Writes a .csv file with the start and end biomass and catch per group
#' from an Rpath.sim object.
#'@export
#Write -- note, not a generic function
write.Rsim <- function(Rsim.output, file, ...){
  out <- data.frame(Group      = Rsim.output$params$spname,
                    StartBio   = Rsim.output$start_state$BB,
                    EndBio     = Rsim.output$end_state$BB,
                    StartCatch = Rsim.output$start_state$CC,
                    EndCatch   = Rsim.output$end_state$CC)
  write.csv(out, file = file)
}
