source("test-constants.R")

#' Plot test results superimposed on the same plot
#' 
#' Plot specific columns from that passed in matrices and plot them as groups of plots per page.
#' The plots will be of the baseline run and the current run plot superimposed on the same plot.
#'
#' @param BaseData : Baseline data to compare current data to
#' @param CurrData : Current data to compare against pre-run baseline data
#' @param baseAlg : The baseline algorithm used (currently AB or RK4)
#' @param currAlg : The current algorithm used (currently AB or RK4)
#' @param tableName : Table name used for the plot's title
#' @param forcedData : The data that's being forced (i.e., Biomass, Catch)
#' @param forcedType : The type of forcing being done (.e., Random (Jitter) or Stair-Stepped)
#' @param species : A vector of fish species
#'
#' @return Returns the final, combined plot
#' 
plotResultsSuperimposed <- function(BaseData,CurrData,baseAlg,currAlg,tableName,forcedData,forcedType,species) {
  plots  <- list()
  group  <- species
  yLabel <- "Biomass (mt/km²)"
  currDf <- data.frame()
  baseDf <- data.frame()
  
  for (member in group) {
    xvalues     <- c(1:length(CurrData[,member]))
    numMonths   <- length(CurrData[,member])
    currYvalues <- CurrData[,member]
    baseYvalues <- BaseData[,member]
    currDf    <- data.frame(Legend=replicate(numMonths,paste0('Current ', currAlg)), xvalues,currYvalues)
    baseDf    <- data.frame(Legend=replicate(numMonths,paste0('Baseline ',baseAlg)), xvalues,baseYvalues)
    
    aPlot <- ggplot() + 
      geom_line(data=baseDf, aes(x=xvalues, y=baseYvalues, color=Legend)) +
      geom_line(data=currDf, aes(x=xvalues, y=currYvalues, color=Legend)) +
      labs(x="Months",y=yLabel,
           title=paste0('Sim Run (',forcedData,' w/ ',forcedType,' Noise) - ',member),
           subtitle=paste0("Dataset: ",tableName)) +
      scale_color_manual(name="Legend:",values=c('red','darkblue')) +
      theme( plot.title = element_text(hjust=0.5,size=7,face="bold"),
             plot.subtitle = element_text(hjust=0.5,size=7),
             axis.text = element_text(size=8),
             axis.title = element_text(size=8),
             legend.text = element_text(size=8),
             legend.title = element_text(size=10),
             legend.position='bottom',
             legend.spacing.y = unit(0.0,'cm'),
             legend.background = element_rect(fill='#f7f7f7'),
             # legend.box.background = element_rect(color = 'black'),
             plot.background = element_rect(color='black',fill=NA,linewidth=1)
      )
    plots <- list.append(plots,aPlot)
    
  }
  combinedPlot <- ggarrange(plotlist=plots,nrow=3,ncol=2)
  # annotate_figure(combinedPlot, top = text_grob("Sample main title here", color = "red", face = "bold", size = 14))
  # saveWidget(ggplotly(combinedPlot), file = "Rplots.html");
  # print(ggplotly(combinedPlot))
  print(combinedPlot)
}

#' Plot the difference of the two runs
#' 
#' Plot specific columns from that passed in matrices and plot them as groups of plots per page.
#' The plots will be of the current run - baseline run. So if there's a perfect match, the plot should
#' be all zeros.
#'
#' @param BaseData : Baseline data to compare current data to
#' @param CurrData : Current data to compare against pre-run baseline data
#' @param baseAlg : The baseline algorithm used (currently AB or RK4)
#' @param currAlg : The current algorithm used (currently AB or RK4)
#' @param tableName : Table name used for the plot's title
#' @param forcedData : The data that's being forced (i.e., Biomass, Catch)
#' @param forcedType : The type of forcing being done (.e., Random (Jitter) or Stair-Stepped)
#' @param species : A vector of fish species
#'
#' @return Returns the final, combined plot
#' 
plotResultsDifference <- function(BaseData,CurrData,baseAlg,currAlg,tableName,forcedData,forcedType,species) {
  plots  <- list()
  group  <- species
  yLabel <- "Biomass (mt/km²)"
  diffDf <- data.frame()
  
  for (member in group) {
    xvalues     <- c(1:length(CurrData[,member]))
    numMonths   <- length(CurrData[,member])
    currYvalues <- CurrData[,member]
    baseYvalues <- BaseData[,member]
    diffDf <- data.frame(Legend=replicate(numMonths,paste0('Current(',currAlg,')-Baseline(',baseAlg,') ')), xvalues,currYvalues-baseYvalues)
    aPlot  <- ggplot() + 
      geom_line(data=diffDf, aes(x=xvalues, y=currYvalues-baseYvalues, color=Legend)) +
      labs(x="Months",y=yLabel,
           title=paste0('Sim Run using ',forcedData,' with ',forcedType,' Noise - ',member),
           subtitle=paste0("Dataset: ",tableName)) +
      scale_color_manual(name="Legend:",values=c('darkblue')) +
      coord_cartesian(ylim = c(-YLIMIT_DIFFERENCE_PLOTS, YLIMIT_DIFFERENCE_PLOTS)) + # RSK
      theme( plot.title = element_text(hjust=0.5,size=7,face="bold"),
             plot.subtitle = element_text(hjust=0.5,size=7),
             axis.text = element_text(size=8),
             axis.title = element_text(size=8),
             legend.text = element_text(size=8),
             legend.title = element_text(size=10),
             legend.position='bottom',
             legend.spacing.y = unit(0.0,'cm'),
             legend.background = element_rect(fill='#f7f7f7'),
             # legend.box.background = element_rect(color = 'black'),
             plot.background = element_rect(color='black',fill=NA,linewidth=1)
      )
    plots <- list.append(plots,aPlot)
  }
  combinedPlot <- ggarrange(plotlist=plots,ncol=1)
  # annotate_figure(combinedPlot, top = text_grob("Sample main title here", color = "red", face = "bold", size = 14))
  print(combinedPlot)
} 