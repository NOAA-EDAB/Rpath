#New function to display foobweb
#Takes an ecopath object as input
foodweb <- function(x){
  require(diagram)
  diet      <- x$DC
  ngroups   <- nrow(diet)
  detgroups <- ngroups - ncol(diet)
  diet      <- cbind(diet, matrix(rep(0, ngroups * detgroups), ngroups, detgroups))
  TL        <- data.table(xpos = as.numeric(NA), TL = x$TL[1:ngroups])
  TL[, TLbin := round(TL)]
  TLmax     <- max(TL[, TLbin])
  counts    <- table(TL[, TLbin])
  TL[, xgap := 1/(counts[names(counts)==TLbin] + 1), by = TLbin]
  posmat    <- TL[TLbin == TLmax, ]
  
  for(i in 1:nrow(posmat)) posmat[i, n := i]
  
  for(i in (TLmax - 1):1){
    TL.bin <- TL[TLbin == i, ]
    for(j in 1:nrow(TL.bin)) TL.bin[j, n := j]
    posmat <- rbindlist(list(posmat, TL.bin))
    }
  
  posmat[, xpos := n * xgap]
  posmat[, c('TLbin', 'xgap', 'n') := NULL, with = F]
  posmat[, TL2 := (TL / TLmax) + (TL - TLmax / 2) * ((1 / TLmax) / 3)]
  plotmat(diet, pos = as.matrix(posmat), name = x$spname, relsize = 0.8, arr.lwd = diet * 5,
          arr.pos = 0.75, box.cex = 1.2, cex.txt = 0.0, box.size = x$BB[1:ngroups]/1000)
  }
  
  plotmat(diet, pos        = as.matrix(posmat), 
                box.size   = 0.001, 
                name       = 1:ngroups,
                lcol       = 'grey',
                cex.txt    = 0.0,
                arr.lwd    = diet * 4,
                arr.length = 0.0,
                box.cex    = .8)
  axis(2, at = axTicks(2), labels = axTicks(2) * TLmax)