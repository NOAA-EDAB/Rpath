#'Plot routine for Ecopath multistanzas
#'
#'Plots the biomass composition of multistanza groups from an Rpath.stanzas object.
#'
#'@family Rpath functions
#'
#'@inheritParams rpath
#'@param StanzaGroup The Stanza group's name to be plotted.
#'@param line.cols A vector of four colors used to represent the population biomass,
#'relative number, indvidual weights, and stanza separation lines.
#'
#'@return Creates a figure showing the break down of biomass and number per stanza.
#'@import data.table
#'@import graphics
#'@export
stanzaplot <- function(Rpath.params, StanzaGroup, line.cols = c('black', 'green', 
                                                                'blue', 'red')){
  #Need to define variables to eliminate check() note about no visible binding
  B <- NageS <- WageS <- age <- B.scale <- NageS.scale <- WageS.scale <- StGroupNum <- Last <- NULL
  
  opar <- par(no.readonly = T)
  
  #Convert StanzaGroup to number
  if(is.character(StanzaGroup)){
    SGNum <- which(Rpath.params$stanzas$stgroups$StanzaGroup == StanzaGroup)
  } else {SGNum <- StanzaGroup}
  
  StGroup <- Rpath.params$stanzas$StGroup[[SGNum]]
  stanza.data <- data.table(B     = StGroup[, B],
                            NageS = StGroup[,NageS],
                            WageS = StGroup[,WageS])
  stanza.data[, age := 0:(nrow(stanza.data) - 1)]
  
  #Scale data between 0 and 1
  stanza.data[, B.scale     := (B     - min(B))     / (max(B)     - min(B))] 
  stanza.data[, NageS.scale := (NageS - min(NageS)) / (max(NageS) - min(NageS))]
  stanza.data[, WageS.scale := (WageS - min(WageS)) / (max(WageS) - min(WageS))]
  
  #Plot the total biomass
  plot(stanza.data[, age], stanza.data[, B.scale], xlab = '', ylab = '', 
       type = 'l', lwd = 3, axes = F, col = line.cols[1])
  
  #Add total number line and weight at age line
  lines(stanza.data[, age], stanza.data[, NageS.scale], lwd = 3, col = line.cols[2])
  lines(stanza.data[, age], stanza.data[, WageS.scale], lwd = 3, col = line.cols[3])
  
  #Add Stanza breaks
  breaks <- Rpath.params$stanzas$stindiv[StGroupNum == SGNum, Last]
  breaks <- breaks[1:(length(breaks) - 1)]
  abline(v = breaks + 1, lwd = 3, col = line.cols[4])
  
  #Add axes, labels, and legend
  axis(1)
  axis(2, las = T)
  box(lwd = 2)
  mtext(1, text = 'Age in Months', line = 2.5)
  mtext(2, text = 'Normalized value', line = 2.5)
  mtext(3, text = Rpath.params$stanzas$stgroups[StGroupNum == SGNum, StanzaGroup], 
        line = 2.3, cex = 2)
  legend('top', legend = c('Population Biomass', 'Number', 'Individual Weight', 
                           'Stanza Separation'), 
         lwd = 2, bty = 'n', col = line.cols, xpd = T, inset = -.15, ncol = 4)
  
  par(opar)
}