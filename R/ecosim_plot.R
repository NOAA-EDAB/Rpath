#'Plot routine for Ecosim runs
#'
#'Plots the relative biomass of each group from a run of ecosim.
#'
#'@family Rpath functions
#'
#'@param Rsim.output Rpath ecosim run created by the rsim.run() function.
#'
#'@return Creates a figure of relative biomass.
#'@import data.table
#'@importFrom grDevices rainbow
#'@export
rsim.plot <- function(Rsim.output, spname="all", indplot = F, ...){
  opar <- par(no.readonly = T)
  if(indplot == F){
    # KYA April 2020 this seems incorrect?
    #biomass <- Rsim.output$out_Biomass[, 2:ncol(Rsim.output$out_Biomass)]
    if("all" %in% spname){spname <- colnames(Rsim.output$out_Biomass[,2:ncol(Rsim.output$out_Biomass)])}
    biomass <- Rsim.output$out_Biomass[, spname]
    n <- ncol(biomass)
    start.bio <- biomass[1, ]
    start.bio[which(start.bio == 0)] <- 1
    rel.bio <- matrix(NA, dim(biomass)[1], dim(biomass)[2])
    for(isp in 1:n) rel.bio[, isp] <- biomass[, isp] / start.bio[isp]
  }
  if(indplot == T){
    spnum <- which(Rsim.output$params$spname == spname)
    biomass <- Rsim.output$out_Biomass[, spnum]
    n <- 1
    rel.bio <- biomass / biomass[1]
  }
  
  ymax <- max(rel.bio) + 0.1 * max(rel.bio)
  ymin <- min(rel.bio) - 0.1 * min(rel.bio)
  ifelse(indplot, xmax <- length(biomass), xmax <- nrow(biomass))
  
  #Plot relative biomass
  opar <- par(mar = c(4, 6, 2, 0))
  
  #Create space for legend
  plot.new()
  l <- legend(0, 0, bty='n', spname, 
              plot=FALSE, fill = line.col, cex = 0.6)
  # calculate right margin width in ndc
  w <- grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
  
  par(omd=c(0, 1-w, 0, 1))
  plot(0, 0, ylim = c(ymin, ymax), xlim = c(0, xmax), 
       axes = F, xlab = '', ylab = '', type = 'n', ...)
  axis(1)
  axis(2, las = T)
  box(lwd = 2)
  mtext(1, text = 'Months', line = 2.5, cex = 1.8)
  mtext(2, text = 'Relative Biomass', line = 3, cex = 1.8)
  
  line.col <- rainbow(n)
  for(i in 1:n){
    if(indplot == T) lines(rel.bio,      col = line.col[i], lwd = 3)
    if(indplot == F) lines(rel.bio[, i], col = line.col[i], lwd = 3)
  }
  
  legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
         spname, fill = line.col, cex = 0.6)
  
  par(opar)
}
