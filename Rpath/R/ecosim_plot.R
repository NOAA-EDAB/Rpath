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
#'@export
ecosim.plot <- function(Rsim.output, spname){
  biomass <- Rsim.output$out_BB[, 2:ncol(Rsim.output$out_BB)]
  n <- ncol(biomass)
  start.bio <- biomass[1, ]
  rel.bio <- matrix(NA, dim(biomass)[1], dim(biomass)[2])
  for(isp in 1:n) rel.bio[, isp] <- biomass[, isp] / start.bio[isp]

  ymax <- max(rel.bio) + 0.1 * max(rel.bio)
  ymin <- min(rel.bio) - 0.1 * min(rel.bio)
  xmax <- nrow(biomass)
  
  #Plot relative biomass
  layout(matrix(c(1, 2), 1, 2), widths = c(4, 1))
  opar <- par(mar = c(4, 6, 2, 0))
  plot(0, 0, ylim = c(ymin, ymax), xlim = c(0, xmax), 
       axes = F, xlab = '', ylab = '', type = 'n')
  axis(1)
  axis(2, las = T)
  box(lwd = 2)
  mtext(1, text = 'Months', line = 2.5, cex = 1.8)
  mtext(2, text = 'Relative Biomass', line = 3, cex = 1.8)
  
  line.col <- rainbow(n)
  for(i in 1:n) lines(rel.bio[, i], col = line.col[i], lwd = 3)
  
  opar <- par(mar = c(0, 0, 0, 0))
  plot(0, 0, xlab = '', ylab = '', axes = F)
  legend('center', legend = spname, fill = line.col, cex = .6)
}
