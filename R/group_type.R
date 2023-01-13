########################################################################################
# Set of functions for returning functional group names (character vector)
# using the type input column
#
# Internal supporting function for group names to check type 
grouptype <- function(Rpath){
  if(class(Rpath)=="Rpath"){gt<-list(type=Rpath$type, grp=Rpath$Group)}
  else{
    if(class(Rpath)=="Rpath.params"){
      gt<-list(type=Rpath$model$Type, grp=Rpath$model$Group)}
    else{
      stop("Input must be an Rpath (balanced) or Rpath.params (unbalanced) object.")
    }
  }
  names(gt$type)<-NULL; names(gt$grp)<-NULL
  return(gt)
}