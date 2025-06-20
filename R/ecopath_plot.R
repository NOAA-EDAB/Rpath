#'Plot routine for Ecopath food web
#'
#'Plots the food web associated with an Rpath object.
#'
#'@family Rpath functions
#'
#'@param Rpath.obj Rpath model created by \code{rpath()}.
#'@param eco.name Optional name of the ecosystem.  Default is the `eco.name` attribute from the
#'  rpath object created from running \code{rpath()}.
#'@param highlight Group number or name to highlight the connections of that group. Valid values are found in the 
#' `Group` field of the object created from running \code{rpath()}.
#'@param highlight.col Color of the connections to the highlighted group, vector of length 3.
#'  Defaults to black for prey, red for predator, orange for fleet.
#'@param labels Logical whether or not to display group names.  If \code{TRUE} and \code{label.pos} = NULL, no 
#'  points will be plotted, just label names.
#'@param label.pos A position specifier for the labels.  Values of 1, 2, 3, 4
#'  indicate positions below, to the left of, above, and to the right of the points, respectively. A null 
#'  value will cause the labels to be plotted without the points, assuming that \code{labels} = TRUE.
#'@param label.num Logical value indication whether group numbers should be used for labels 
#'  instead of names.
#'@param line.col The color of the lines between nodes of the food web.
#'@param fleets Logical value indicating whether or not to include fishing fleets in the food web plot.
#'@param type.col The color of the points corresponding to the types of the group.  Must be 
#'  of length 1 or 4.  Color order will be consumers, primary producers, detrital, and fleet groups.  
#'@param box.order Vector of box numbers to change the default plot order.  Must include all box numbers
#'@param label.cex Numeric value of the relative size of the labels within the plot.
#'
#'@return Returns a plot visualization of the food web.
#'
#'@examples
#' # Read in Rpath parameter file, generate and name model object
#' Rpath.obj <- rpath(AB.params, eco.name = "Anchovy Bay")
#' # Plot food web diagram with all groups labeled, including fleets
#' webplot(Rpath.obj, labels = TRUE, fleets = TRUE)
#' # Plot food web diagram without labels, highlighting connections of cod group
#' webplot(Rpath.obj, highlight = "cod",fleets = TRUE)
#'
#'@import data.table
#'@import graphics
#'@import utils
#'@export
webplot <- function(Rpath.obj, eco.name = attr(Rpath.obj, 'eco.name'), line.col = 'grey',
                    highlight = NULL, highlight.col = c('black', 'red', 'orange'), 
                    labels = FALSE, label.pos = NULL, label.num = FALSE, label.cex = 1,
                    fleets = FALSE, type.col = 'black', box.order = NULL) {
  web.info <- summarize.for.webplot(Rpath.obj, fleets, box.order)
  
  list2env(web.info, sys.frame(sys.nframe()))
  
  opar <- par(no.readonly = T)
  ymin <- min(pointmap[, TL]) - 0.1 * min(pointmap[, TL])
  ymax <- max(pointmap[, TL]) + 0.1 * max(pointmap[, TL])
  plot(0, 0, ylim = c(ymin, ymax), xlim = c(0, 1), typ = 'n', xlab = '', 
       ylab = '', axes = F)
  if(!is.null(eco.name)) mtext(3, text = eco.name, cex = 1.5)
  axis(2, las = T)
  box()
  mtext(2, text = 'Trophic Level', line = 2)
  
  with(connections,
       segments(x0 = pred.x, y0 = pred.y, x1 = prey.x, y1 = prey.y, col = line.col))
  
  if(!is.null(highlight)){
    if(is.character(highlight)) highlight <- which(Rpath.obj$Group == highlight)
    pred.x <- pointmap[GroupNum == highlight, x.pos]
    pred.y <- pointmap[GroupNum == highlight, TL]
    tot.catch <- Rpath.obj$Landings + Rpath.obj$Discards
    
    if(pointmap[GroupNum == highlight, type] < 1){
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
    points(pointmap[type <  1, x.pos], pointmap[type <  1, TL], pch = 16, col = type.col[1])
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
  
  par(opar)  
}

#' ggplot routine for Rpath food web
#' 
#' Plots the food web associated with an Rpath object using ggplot functions
#'
#' @param Rpath.obj Rpath model created by the \code{rpath()} function.
#' @param eco.name Optional name of the ecosystem. Default is the `eco.name` 
#'                 attribute from the rpath object.
#' @param line.col The color of the lines between nodes of the food web.
#' @param highlight Set to the group number or name to highlight the connections of that group.
#'  Valid values are found in the `Group` field of the object created from running \code{rpath()}.
#' @param highlight.col Color of the connections to the highlighted group, vector of length 3.
#'  Defaults to black = predator, red = prey, orange = fleet.
#' @param labels Logical whether or not to display group names.  
#' @param label.num Logical whether or not to display group numbers instead of points at nodes. 
#' If \code{TRUE}, \code{type.col} must be length 1, not 4.
#' @param label.cex Numeric value of the relative size of the labels within the plot.
#' @param fleets Logical value indicating whether or not to include fishing fleets in the food web.
#' @param type.col The color of the points corresponding to the types of the group.  Can either be 
#'  of length 1 or 4. Color order will be living, primary producers, detrital, and fleet groups.  
#' @param box.order Vector of box numbers to change the default plot order. Must include all box numbers. 
#' Passed to \code{summarize.for.webplot()}
#' @param line.alpha Transparency of lines between nodes of the food web.
#' @param point.size Size of points at nodes. 
#' @param text.size Size of text
#' @param max.overlaps Maximum number of overlaps allowed for group labels by \code{ggrepel}
#'
#' @return Returns a ggplot object visualizing the food web
#' 
#' @examples
#' # Read in Rpath parameter file, generate and name model object
#' Rpath.obj <- rpath(AB.params, eco.name = "Anchovy Bay")
#' # Plot food web diagram with all groups labeled, including fleets, using ggplot
#' ggwebplot(Rpath.obj, labels = TRUE, fleets = TRUE)
#' # Plot food web diagram without labels, highlighting connections of cod group
#' ggwebplot(Rpath.obj, highlight = "cod",fleets = TRUE)
#' 
#' @import data.table
#' @import ggplot2
#' @import utils
#' @export
#'
ggwebplot <- function(Rpath.obj, eco.name = attr(Rpath.obj, 'eco.name'), line.col = 'grey',
                      highlight = NULL, highlight.col = c('black', 'red', 'orange'), 
                      labels = FALSE, label.num = FALSE, label.cex = 1,
                      fleets = FALSE, type.col = 'grey50', box.order = NULL, 
                      line.alpha = 0.5, point.size = 1, text.size = 5,
                      max.overlaps = 10) {
  # build pointmap and set of connections
  web.info <- summarize.for.webplot(Rpath.obj, fleets, box.order)
  list2env(web.info, sys.frame(sys.nframe()))
  
  # plot connections
  p <- ggplot() +
    geom_segment(aes(x = pred.x, y = pred.y, xend = prey.x, yend = prey.y),
                 col = line.col, alpha = line.alpha, data = connections) +
    labs(x = '', y = 'Trophic position') +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) 
  
  if(!is.na(eco.name)) p <- p + ggtitle(eco.name)
  
  # plot nodes
  if(length(type.col) == 1) {
    if(label.num) {
      p <- p + geom_text(aes(x = x.pos, y = TL, label = GroupNum), size = text.size,
                         col = type.col, data = pointmap)
    } else {
      p <- p + geom_point(aes(x = x.pos, y = TL), col = type.col, size = point.size,
                          data = pointmap)
    }
  } else {
    p <- p + geom_point(aes(x = x.pos, y = TL, fill = factor(type)),
                        pch = 21, size = point.size, data = pointmap) +
      scale_fill_manual(values = type.col, name = 'Type', 
                        labels = c('Living', 'Primary', 'Detrital', 'Fleet'))
  }
  
    # label nodes
    if(labels) {
      p <- p + ggrepel::geom_text_repel(aes(x = x.pos, y = TL, label = Group), 
                                        max.overlaps = max.overlaps, 
                                        size = text.size, data = pointmap)
    }
    
    # highlight group
    if(!is.null(highlight)) {
      if(is.character(highlight)) highlight <- which(Rpath.obj$Group == highlight)
      
      # make the highlight data frame
      highlight.df <- connections[prey.id == highlight | pred.id == highlight][,
                                  type := factor(ifelse(pred.id == highlight, 'predator',
                                                        ifelse(pred.id %in% pointmap[type==3, GroupNum], 
                                                               'fleet', 'prey')),
                                                 levels = c('predator', 'prey', 'fleet'))]
      p <- p + geom_segment(aes(x = pred.x, y = pred.y, xend = prey.x, yend = prey.y,
                                col = type),
                            data = highlight.df) +
        scale_color_manual(values = highlight.col, name = 'Connection')
    }
    return(p)
}
  
#' Prepare food web model to plot
#' 
#' Prepare and summarize food web model for plotting by generating x-values for node 
#' locations and creating a data frame of all trophic connections. This function could also
#' be used to build custom plotting functions. 
#' 
#' @keywords internal
#'
#' @param Rpath.obj Rpath model created by the \code{rpath()} function.
#' @param fleets Logical value indicating whether or not to include fishing fleets in the food web.
#' @param box.order Vector of box numbers to change the default plot order. Must include all box numbers
#'
#' @return Returns a list of length 2. The first element, \code{pointmap} is a data table of
#' the functional group names, group numbers, type (0- living, 1- producer, 2- detritus, 3- fleet),
#' trophic level, biomass, and \code{x.pos} which determines the x coordinate of the node 
#' when plotting. The second element is a data table with one row for each trophic connection, 
#' and columns for predator and prey trophic levels, predator and prey \code{x.pos} from \code{pointmap},
#' and predator and prey group numbers (\code{id}). 
#' @import data.table
#' @import utils
#' @export
#'
summarize.for.webplot <- function(Rpath.obj, fleets = FALSE, box.order = NULL){
  #Need to define variables to eliminate check() note about no visible binding
  TL <- TLlevel <- type <- n <- x.space <- x.offset <- Group <- x.pos <- GroupNum <- NULL
  
  pointmap <- data.table(GroupNum = 1:length(Rpath.obj$TL), 
                         Group    = Rpath.obj$Group, 
                         type     = Rpath.obj$type, 
                         TL       = Rpath.obj$TL, 
                         Biomass  = Rpath.obj$Biomass)
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
    if(length(x.count[, Group] > 0)){
      for(j in 1:nrow(x.count)){
        x.count[j, x.count := j] 
      }
      x.count.all <- rbind(x.count.all, x.count)
    }
  }
  pointmap <- merge(pointmap, x.count.all, by = 'Group', all.x = T)
  pointmap[x.count == 1, x.pos := x.offset + rnorm(1, 0, 0.01)]
  pointmap[x.count != 1, x.pos := x.offset + x.space * (x.count - 1) + rnorm(1, 0, 0.01)]
  pointmap[, c('TLlevel', 'n', 'x.offset', 'x.space', 'x.count') := NULL]
  
  #Web connections
  tot.catch <- Rpath.obj$Landings + Rpath.obj$Discards
  pred      <- pointmap[!type %in% 1:2, GroupNum]
  
  connection.mat <- matrix(NA, nrow = sum(Rpath.obj$DC != 0) + sum(tot.catch != 0) * fleets, 
                           ncol = 6,
                           dimnames = list(pair = NULL, 
                                           position = c('pred.x', 'pred.y',
                                                        'prey.x', 'prey.y',
                                                        'pred.id', 'prey.id')))
  row.num <- 0
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
    nprey <- length(prey)
    prey.x <- pointmap[GroupNum %in% prey, x.pos]
    prey.y <- pointmap[GroupNum %in% prey, TL]
    prey.id <- pointmap[GroupNum %in% prey, GroupNum]
    
    connection.mat[row.num + 1:nprey, 'pred.x'] <- pred.x
    connection.mat[row.num + 1:nprey, 'pred.y'] <- pred.y
    connection.mat[row.num + 1:nprey, 'prey.x'] <- prey.x
    connection.mat[row.num + 1:nprey, 'prey.y'] <- prey.y
    connection.mat[row.num + 1:nprey, 'pred.id'] <- i
    connection.mat[row.num + 1:nprey, 'prey.id'] <- prey.id
    
    row.num <- row.num + nprey
  }
  
  return(list(pointmap = pointmap, 
              connections = data.table(connection.mat)))  
}

#'Plot routine for Ecopath multistanzas
#'
#'Plots the biomass composition of multistanza groups from an Rpath.stanzas object.
#'
#'@family Rpath functions
#'
#'@inheritParams rpath
#'@param StanzaGroup The Stanza group's name to be plotted.Valid values are found in the 
#' `stanzas` field of the Rpath parameter file.
#'@param line.cols A vector of four colors used to represent the population biomass,
#'relative number, individual weights, and stanza separation lines.
#'
#'@return Creates a figure showing the break down of biomass and number per stanza.
#'
#'@examples 
#' # Choose group with multiple stanzas
#' params <- REco.params
#' stanzaplot(params, StanzaGroup = "Roundfish1")
#' 
#'
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
