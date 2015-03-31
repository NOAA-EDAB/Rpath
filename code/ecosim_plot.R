ecosim.plot <- function(x){
  biomass <- as.data.table(x$out_BB[, 2:23])
  rel.bio <- data.table(Month = 1:nrow(biomass))
  for(i in 1:ncol(biomass)){
    sp.bio.start <- biomass[1, i, with = F]
    sp.rel.bio   <- biomass[ , i, with = F] / as.numeric(sp.bio.start)
    rel.bio <- cbind(rel.bio, sp.rel.bio)
  }
  
  ymax <- max(rel.bio[, 2:ncol(rel.bio), with = F])
  ymax <- ymax + 0.1 * ymax
  ymin <- min(rel.bio[, 2:ncol(rel.bio), with = F])
  ymin <- ymin - 0.1 * ymin
  xmax <- max(rel.bio[, Month])
  
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
  
  line.col <- rainbow(ncol(rel.bio) - 1)
  for(i in 1:(ncol(rel.bio) - 1)){
    lines(rel.bio[, c(1, i + 1), with = F], col = line.col[i], lwd = 3)
  }
  
  opar <- par(mar = c(0, 0, 0, 0))
  plot(0, 0, xlab = '', ylab = '', axes = F)
  legend('center', legend = x$spname[2:23], fill = line.col, cex = .6)
}
