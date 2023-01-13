#'Write Rpath parameters to .csv files
#'
#'Creates a series of .csv files from an Rpath.params object.
#'
#'@family Rpath functions
#'
#'@param Rpath.params R object containing the Rpath parameters.  Most likely this
#'  was created using create.rpath.params or read.rpath.params.
#'@param eco.name ecosystem name that will be included in all the file names.
#'@param path location for the output files.  
#'@return Outputs a series of .csv files named by the provided eco.name and the 
#'  parameters they represent.  For example the model parameters will be named 
#'  "eco.name_model.csv".
#'@export
write.rpath.params <- function(Rpath.params, eco.name, path = ''){
  
  write.csv(Rpath.params$model,
            file = file.path(path, paste(eco.name, '_model.csv', sep = '')), 
            row.names = F)
  
  write.csv(Rpath.params$diet,
            file = file.path(path, paste(eco.name, '_diet.csv', sep = '')), 
            row.names = F)
  
  write.csv(Rpath.params$pedigree, 
            file = file.path(path, paste(eco.name, '_pedigree.csv', sep = '')), 
            row.names = F)
  
  #Multistanza parameters are in several different files
  write.csv(Rpath.params$stanzas$stgroups, 
            file = file.path(path, paste(eco.name, '_stanza_groups.csv', sep = '')),
            row.names = F)
  
  write.csv(Rpath.params$stanzas$stindiv,  
            file = file.path(path, paste(eco.name, '_stanzas.csv', sep = '')),
            row.names = F)
  
  if(Rpath.params$stanzas$NStanzaGroups > 0){
    for(isp in 1:Rpath.params$stanzas$NStanzaGroups){
      write.csv(Rpath.params$stanzas$StGroup[[isp]], 
                file = file.path(path, paste(eco.name, '_', 
                                             Rpath.params$stanzas$stgroups$StanzaGroup[isp], 
                                             '.csv', sep = '')), row.names = F)
    }
  }
}