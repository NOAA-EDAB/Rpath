#Functions for Rpath objects
summary.Rpath <- function(object, ...){
  z <- object
  EE <- z$EE
  if(max(EE, na.rm = T) > 1){
    balance <- "Unbalanced"
    unbalanced.groups <- z$Group[which(EE > 1)]
  }else{
    balance <- "Balanced"
  }
  totbiomass <- sum(z$BB,    na.rm = T)
  totcatch   <- sum(z$Catch, na.rm = T)
  table <- data.table(Num.Groups   = z$NUM_GROUPS,
                      Num.Living   = z$NUM_LIVING,
                      Num.Detritus = z$NUM_DEAD,
                      Num.Fleets   = z$NUM_GEARS,
                      TotBiomass   = totbiomass,
                      TotCatch     = totcatch)
  ans <- list(model = balance,
              stats = table)
  return(ans)
}
