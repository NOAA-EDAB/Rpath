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
#'@import utils
#'@export
#Write -- note, not a generic function
write.Rsim <- function(Rsim.output, file = NA){
  #Check extention for output type - if not supplied a .csv or .RData the output
  #will be saved as a list object
  ext <- 'list'
  if(grepl('.csv', file) == T) ext <- 'csv'
  if(grepl('.RData', file, ignore.case = T) == T) ext <- 'RData'
  
  gear.zero <- rep(0, Rsim.output$params$NUM_GEARS)
  start_Catch <- c(Rsim.output$out_Catch[2, ], gear.zero)
  end_Catch   <- c(Rsim.output$out_Catch[nrow(Rsim.output$out_Catch) - 1, ], gear.zero)
  out <- data.frame(Group      = Rsim.output$params$spname,
                    StartBio   = Rsim.output$start_state$Biomass,
                    EndBio     = Rsim.output$end_state$Biomass,
                    BioES      = Rsim.output$end_state$Biomass / 
                      Rsim.output$start_state$Biomass,
                    StartCatch = start_Catch * 12,
                    EndCatch   = end_Catch * 12,
                    CatchES    = (end_Catch * 12) / (start_Catch * 12))
  if(ext == 'csv')   write.csv(out, file = file)
  if(ext == 'RData') save(out, file = file)
  if(ext == 'list')  return(out)
}