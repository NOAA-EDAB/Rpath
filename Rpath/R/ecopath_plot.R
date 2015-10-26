#'Plot routine for Ecopath food web
#'
#'Plots the food web associated with an Rpath object.
#'
#'@family Rpath functions
#'
#'@param Rpath.obj Rpath model created by the ecopath() function.
#'@param highlight Box number to highlight connections.
#'@param eco.name Optional name of the ecosystem.  Default is the eco.name attribute from the
#'    rpath object.
#'@param highlight Set to the group number or name to highlight the connections of that group.
#'@param highlight.col Color of the connections to the highlighted group.
#'@param labels Logical whether or not to display group names.  If True and label.pos is Null, no 
#'  points will be ploted, just label names.
#'@param label.pos A position specifier for the labels.  Values of 1, 2, 3, 4, respectively 
#'  indicate positions below, to the left of, above, and to the right of the points. A null 
#'  value will cause the labels to be ploted without the points (Assuming that labels = TRUE).
#'@param label.num Logical value indication whether group numbers should be used for labels 
#'  instead of names.
#'@param line.col The color of the lines between nodes of the food web.
#'@param fleets Logical value indicating whether or not to include fishing fleets in the food web.
#'@param type.col The color of the points cooresponding to the types of the group.  Can either be 
#'  of length 1 or 4.  Color order will be living, primary producers, detrital, and fleet groups.  
#'@param box.order Vector of box numbers to change the default plot order.  Must include all box numbers
#'@param label.cex The relative size of the labels within the plot.
#'
#'@return Creates a figure of the food web.
#'@import data.table
#'@export
webplot <- function(Rpath.obj, eco.name = attr(Rpath.obj, 'eco.name'), line.col = 'grey',
                    highlight = NULL, highlight.col = c('black', 'red', 'orange'), 
                    labels = FALSE, label.pos = NULL, label.num = FALSE, label.cex = 1,
                    fleets = FALSE, type.col = 'black', box.order = NULL){
  pointmap <- data.table(GroupNum = 1:length(Rpath.obj$TL), 
                         Group    = Rpath.obj$Group, 
                         type     = Rpath.obj$type, 
                         TL       = Rpath.obj$TL, 
                         Biomass  = Rpath.obj$BB)
  pointmap[TL < 2,               TLlevel := 1]
  pointmap[TL >= 2.0 & TL < 3.0, TLlevel := 2]
  pointmap[TL >= 3.0 & TL < 3.5, TLlevel := 3]
  pointmap[TL >= 3.5 & TL < 4.0, TLlevel := 4]
  pointmap[TL >= 4.0 & TL < 4.5, TLlevel := 5]
  pointmap[TL >= 4.5 & TL < 5.0, TLlevel := 6]
  pointmap[TL >= 5.0,            TLlevel := 7]
  
  if(!is.null(box.order)) pointmap <- pointmap[box.order, ]
  
  if(fleets == F) pointmap <- pointmap[type < 3, ]
  nTL <- table(pointmap[, TLlevel])
  pointmap[, n := nTL[which(names(nTL) == TLlevel)], by = TLlevel]
  pointmap[, x.space  := 1 / n]
  pointmap[, x.offset := x.space / 2]
  x.count.all <- c()
  for(i in 1:max(pointmap[, TLlevel])){
    x.count <- pointmap[TLlevel == i, list(Group)]
    for(j in 1:nrow(x.count)){
      x.count[j, x.count := j]  
    }
    x.count.all <- rbind(x.count.all, x.count)
  }
  pointmap <- merge(pointmap, x.count.all, by = 'Group', all.x = T)
  pointmap[x.count == 1, x.pos := x.offset + rnorm(1, 0, 0.01)]
  pointmap[x.count != 1, x.pos := x.offset + x.space * (x.count - 1) + rnorm(1, 0, 0.01)]
  pointmap[, c('TLlevel', 'n', 'x.offset', 'x.space', 'x.count') := NULL]
  
  ymin <- min(pointmap[, TL]) - 0.1 * min(pointmap[, TL])
  ymax <- max(pointmap[, TL]) + 0.1 * max(pointmap[, TL])
  plot(0, 0, ylim = c(ymin, ymax), xlim = c(0, 1), typ = 'n', xlab = '', 
       ylab = '', axes = F)
  if(!is.null(eco.name)) mtext(3, text = eco.name, cex = 1.5)
  axis(2, las = T)
  box()
  mtext(2, text = 'Trophic Level', line = 2)
  
  #Web connections
  tot.catch <- Rpath.obj$Catch + Rpath.obj$Discards
  pred      <- pointmap[type %in% c(0, 3), GroupNum]
  
  for(i in pred){
    pred.x <- pointmap[GroupNum == i, x.pos] 
    pred.y <- pointmap[GroupNum == i, TL]
    if(pointmap[GroupNum == i, type] == 0){
      prey <- which(Rpath.obj$DC[, i] > 0)
    }
    if(pointmap[GroupNum == i, type] == 3){
      gear.num <- i - (Rpath.obj$NUM_GROUPS - Rpath.obj$NUM_GEARS)
      prey <- which(tot.catch[, gear.num] > 0)
    }
    prey.x <- pointmap[GroupNum %in% prey, x.pos]
    prey.y <- pointmap[GroupNum %in% prey, TL]
    for(j in 1:length(prey)){
      lines(c(pred.x, prey.x[j]), c(pred.y, prey.y[j]), col = line.col)
    }
  }
  
  if(!is.null(highlight)){
    if(is.character(highlight)) highlight <- which(Rpath.obj$Group == highlight)
    pred.x <- pointmap[GroupNum == highlight, x.pos]
    pred.y <- pointmap[GroupNum == highlight, TL]
    if(pointmap[GroupNum == highlight, type] == 0){
      prey       <- which(Rpath.obj$DC[, highlight] > 0)
      group.pred <- which(Rpath.obj$DC[highlight, ] > 0)
      fleet.pred <- which(tot.catch[highlight, ] > 0)
    }
    if(pointmap[GroupNum == highlight, type] %in% c(1:2)){
      prey       <- NULL
      group.pred <- which(Rpath.obj$DC[highlight, ] > 0)
      fleet.pred <- which(tot.catch[highlight, ] > 0)
    } 
    if(pointmap[GroupNum == highlight, type] == 3){
      gear.num   <- highlight - (Rpath.obj$NUM_GROUPS - Rpath.obj$NUM_GEARS)
      prey       <- which(tot.catch[, gear.num] > 0)
      group.pred <- NULL
      fleet.pred <- NULL
    }
    if(!is.null(prey)){
      prey.x <- pointmap[GroupNum %in% prey, x.pos]
      prey.y <- pointmap[GroupNum %in% prey, TL]
      for(j in 1:length(prey)){
        lines(c(pred.x, prey.x[j]), c(pred.y, prey.y[j]), col = highlight.col[1], lwd = 2)
      }
    }
    if(!is.null(group.pred)){
      group.pred.x <- pointmap[GroupNum %in% group.pred, x.pos]
      group.pred.y <- pointmap[GroupNum %in% group.pred, TL]
      for(j in 1:length(group.pred)){
        lines(c(pred.x, group.pred.x[j]), c(pred.y, group.pred.y[j]), 
              col = highlight.col[2], lwd = 2)
      }
    }
    if(length(fleet.pred) > 0){
      gear.num <- fleet.pred + (Rpath.obj$NUM_GROUPS - Rpath.obj$NUM_GEARS)
      fleet.pred.x <- pointmap[GroupNum %in% gear.num, x.pos]
      fleet.pred.y <- pointmap[GroupNum %in% gear.num, TL]
      for(j in 1:length(fleet.pred)){
        lines(c(pred.x, fleet.pred.x[j]), c(pred.y, fleet.pred.y[j]), 
              col = highlight.col[3], lwd = 2)
      }
    }
    legend('bottomleft', legend = c('prey', 'predator', 'fleet'), lty = 1, col = highlight.col,
           lwd = 2, ncol = 3, xpd = T, inset = c(0, -.1))
    legend('topright', legend = pointmap[GroupNum == highlight, Group], bty = 'n')
  }
  
  #Group points
  if(!is.null(label.pos) | labels == F){
    if(length(type.col) ==4){
      legend('bottomright', legend = c('living', 'primary', 'detrital', 'fleet'), 
             pch = 16, col = type.col, ncol = 4, xpd = T, inset = c(0, -.1))
    }
    if(length(type.col) < 4) type.col <- rep(type.col[1], 4)
    points(pointmap[type == 0, x.pos], pointmap[type == 0, TL], pch = 16, col = type.col[1])
    points(pointmap[type == 1, x.pos], pointmap[type == 1, TL], pch = 16, col = type.col[2])
    points(pointmap[type == 2, x.pos], pointmap[type == 2, TL], pch = 16, col = type.col[3])
    points(pointmap[type == 3, x.pos], pointmap[type == 3, TL], pch = 16, col = type.col[4])
  }
  
  if(labels == T){
    if(label.num == F){
      text(pointmap[, x.pos], pointmap[, TL], pointmap[, Group], 
           pos = label.pos, cex = label.cex)
    }
    if(label.num == T){
      text(pointmap[, x.pos], pointmap[, TL], pointmap[, GroupNum], 
           pos = label.pos, cex = label.cex)
    }
  }
  
}

#'Plot routine for Ecopath multistanzas
#'
#'Plots the biomass composition of multistanza groups from an Rpath.stanzas object.
#'
#'@family Rpath functions
#'
#'@param Rpath.stanza Rpath.stanza object created by the rpath.stanza() function.
#'@param stanzafile Object containing the stanza specific characteristics use in the
#'rpath.stanza.stanza() function.
#'@param StanzaNum The Stanza group's number in the stanzafile.
#'@param line.cols A vector of four colors used to represent the population biomass,
#'relative number, indvidual weights, and stanza separation lines.
#'
#'@return Creates a figure of the food web.
#'@import data.table
#'@export
stanzaplot <- function(Rpath.stanza, stanzafile, StanzaNum, 
                       line.cols = c('black', 'green', 'blue', 'red')){
  stanza.data <- copy(Rpath.stanza[[StanzaNum]])
  #Scale data between 0 and 1
  stanza.data[, lawa.scale := (lawa - min(lawa)) / (max(lawa) - min(lawa))] 
  stanza.data[, la.scale := (la - min(la)) / (max(la) - min(la))]
  stanza.data[, wa.scale := (wa - min(wa)) / (max(wa) - min(wa))]
  
  #Plot the total biomass
  plot(stanza.data[, age], stanza.data[, lawa.scale], xlab = '', ylab = '', 
       type = 'l', lwd = 3, axes = F, col = line.col[1])
  
  #Add total number line and weight at age line
  lines(stanza.data[, age], stanza.data[, la], lwd = 3, col = line.col[2])
  lines(stanza.data[, age], stanza.data[, wa], lwd = 3, col = line.col[3])
  
  #Add Stanza breaks
  breaks <- stanzafile[StGroupNum == StanzaNum, Last]
  breaks <- breaks[1:(length(breaks) - 1)]
  abline(v = breaks, lwd = 3, col = line.col[4])
  
  #Add axes, labels, and legend
  axis(1)
  axis(2, las = T)
  box(lwd = 2)
  mtext(1, text = 'Age in Months', line = 2.5)
  mtext(2, text = 'Normalized value', line = 2.5)
  legend('top', legend = c('Population Biomass', 'Number', 'Individual Weight', 
                                'Stanza Separation'), 
         lwd = 2, bty = 'n', col = line.col, xpd = T, inset = -.2, ncol = 4)
}
