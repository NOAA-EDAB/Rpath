#R Ecosystem
#Comparison to EwE

#User parameters - file locations
if(Sys.info()['sysname']=="Windows"){
  #data.dir <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\data\\"
  #out.dir  <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\outputs\\"
}

if(Sys.info()['sysname']=="Linux"){
  rpath.dir <- "/home/slucey/slucey/Rpath/outputs"
  ewe.dir   <- "/home/slucey/slucey/Manuscripts/Rpath/EwE_files"
  out.dir   <- "/home/slucey/slucey/Manuscripts/Rpath/Comparisons"
}

library(Rpath); library(data.table)

#-------------------------------------------------------------------------------
#User functions
EwE.summary <- function(directory){
  EwE.bio   <- as.data.table(read.csv(file.path(directory, 'Biomass_annual.csv'), 
                                      skip = 9))
  EwE.catch <- as.data.table(read.csv(file.path(directory, 'Yield_annual.csv'), 
                                      skip = 9))
  #Create vectors of values
  Groups     <- as.data.table(names(EwE.bio))
  StartBio   <- transpose(EwE.bio[1, ])
  EndBio     <- transpose(EwE.bio[nrow(EwE.bio) - 1, ])
  BioES      <- EndBio / StartBio
  StartCatch <- transpose(EwE.catch[1, ])
  EndCatch   <- transpose(EwE.catch[nrow(EwE.catch) - 1, ])
  CatchES    <- EndCatch / StartCatch
  
  setnames(Groups,     1, 'Group')
  setnames(StartBio,   1, 'StartBio')
  setnames(EndBio,     1, 'EndBio')
  setnames(BioES,      1, 'BioES')
  setnames(StartCatch, 1, 'StartCatch')
  setnames(EndCatch,   1, 'EndCatch')
  setnames(CatchES,    1, 'CatchES')
  
  EwE.sim <- as.data.table(cbind(Groups, StartBio, EndBio, BioES,
                                 StartCatch, EndCatch, CatchES))
  return(EwE.sim)
}
#-------------------------------------------------------------------------------
#Compare Ecopath
Rpath <- as.data.table(read.csv(file.path(rpath.dir, 
                                          'R_Ecosystem_Parameters.csv')))
EwE   <- as.data.table(read.csv(file.path(ewe.dir, 'REco-Basic estimates.csv')))

#Fix EwE columns
setnames(EwE, c('Group.name', 'Trophic.level', 'Biomass..t.km..', 
                'Production...biomass...year.', 'Consumption...biomass...year.',
                'Ecotrophic.efficiency', 'Production...consumption'),
                c('Group', 'TL', 'Biomass', 'PB', 'QB', 'EE', 'GE'))
EwE[is.na(PB), PB := Z...year.]

#Remove extra rows/columns
EwE <- EwE[!is.na(X), ]
EwE[, c('Habitat.area..fraction.', 'Biomass.in.habitat.area..t.km..',
        'Z...year.', 'X.1') := NULL]
Rpath <- Rpath[type != 3, ]
Rpath[, c('type', 'Removals') := NULL]

TL  <- ((EwE[, TL] - Rpath[, TL]) / ((EwE[, TL] + Rpath[, TL]) / 2)) * 100
Bio <- ((EwE[, Biomass] - Rpath[, Biomass]) / 
          ((EwE[, Biomass] + Rpath[, Biomass]) / 2)) * 100
PB <- ((EwE[, PB] - Rpath[, PB]) / ((EwE[, PB] + Rpath[, PB]) / 2)) * 100
QB <- ((EwE[, QB] - Rpath[, QB]) / ((EwE[, QB] + Rpath[, QB]) / 2)) * 100
EE <- ((EwE[, EE] - Rpath[, EE]) / ((EwE[, EE] + Rpath[, EE]) / 2)) * 100
GE <- ((EwE[, GE] - Rpath[, GE]) / ((EwE[, GE] + Rpath[, GE]) / 2)) * 100

ecopath.outputs <- cbind(TL, Bio, PB, QB, EE, GE)

png(file = file.path(out.dir, 'Ecopath_attributes.png'), height = 1700,
    width = 2000, res = 200)
opar <- par(mar = c(5, 7, 1, 1))
boxplot(ecopath.outputs, axes = F, cex = 3)
box(lwd = 3)
axis(1, at = axTicks(1), labels = c('TL','Bio', 'PB', 'QB', 'EE', 'GE'),
     cex.axis = 2)
axis(2, las = T, cex.axis = 2)
mtext(1, text = 'Ecopath attributes', line = 3, cex = 2)
mtext(2, text = '% difference', line = 5 , cex = 2)
dev.off()

#-------------------------------------------------------------------------------
#Ecosim
#Base Effort
#Rpath
Rpath.sim <- as.data.table(read.csv(file.path(rpath.dir, 'REco_baserun_AB.csv')))
#Remove extra rows
Rpath.sim <- Rpath.sim[!Group %in% c('Outside', 'Detritus', 'Discards', 
                                     'Trawlers', 'Midwater', 'Dredgers')]

#EwE
base.dir <- file.path(ewe.dir, 'Ecosim_REco Base')
EwE.sim  <- EwE.summary(base.dir)
#Remove extra rows
EwE.sim <- EwE.sim[!Group %in% c('Detritus', 'Discards'), ]

Start.bio   <- ((EwE.sim[, StartBio] - Rpath.sim[, StartBio]) / 
                  ((EwE.sim[, StartBio] + Rpath.sim[, StartBio]) / 2)) * 100
End.bio     <- ((EwE.sim[, EndBio] - Rpath.sim[, EndBio]) / 
                  ((EwE.sim[, EndBio] + Rpath.sim[, EndBio]) / 2)) * 100
Start.catch <- ((EwE.sim[, StartCatch] - Rpath.sim[, StartCatch]) / 
                  ((EwE.sim[, StartCatch] + Rpath.sim[, StartCatch]) / 2)) * 100
End.catch   <- ((EwE.sim[, EndCatch] - Rpath.sim[, EndCatch]) / 
                  ((EwE.sim[, EndCatch] + Rpath.sim[, EndCatch]) / 2)) * 100

ecosim.outputs <- cbind(Start.bio, End.bio, Start.catch, End.catch)

png(file = file.path(out.dir, 'Ecosim_attributes_baserun.png'), height = 1700,
    width = 2000, res = 200)
opar <- par(mar = c(5, 7, 1, 1))
boxplot(ecosim.outputs, axes = F, cex = 3)
box(lwd = 3)
axis(1, at = axTicks(1), labels = c('Start Bio','End Bio', 'Start Catch', 
                                    'End Catch'), cex.axis = 2)
axis(2, las = T, cex.axis = 2)
mtext(1, text = 'Ecosim statistics', line = 3, cex = 2)
mtext(2, text = '% difference', line = 5 , cex = 2)
dev.off()

#Double Trawling Effort
#Rpath
Rpath.sim <- as.data.table(read.csv(file.path(rpath.dir, 'REco_doubletrawl_AB.csv')))
#Remove extra rows
Rpath.sim <- Rpath.sim[!Group %in% c('Outside', 'Detritus', 'Discards', 
                                     'Trawlers', 'Midwater', 'Dredgers')]

#EwE
double.dir <- file.path(ewe.dir, 'Ecosim_Doubletrawl')
EwE.sim  <- EwE.summary(double.dir)
#Remove extra rows
EwE.sim <- EwE.sim[!Group %in% c('Detritus', 'Discards'), ]

Start.bio   <- ((EwE.sim[, StartBio] - Rpath.sim[, StartBio]) / 
                  ((EwE.sim[, StartBio] + Rpath.sim[, StartBio]) / 2)) * 100
End.bio     <- ((EwE.sim[, EndBio] - Rpath.sim[, EndBio]) / 
                  ((EwE.sim[, EndBio] + Rpath.sim[, EndBio]) / 2)) * 100
Start.catch <- ((EwE.sim[, StartCatch] - Rpath.sim[, StartCatch]) / 
                  ((EwE.sim[, StartCatch] + Rpath.sim[, StartCatch]) / 2)) * 100
End.catch   <- ((EwE.sim[, EndCatch] - Rpath.sim[, EndCatch]) / 
                  ((EwE.sim[, EndCatch] + Rpath.sim[, EndCatch]) / 2)) * 100

ecosim.outputs <- cbind(Start.bio, End.bio, Start.catch, End.catch)

png(file = file.path(out.dir, 'Ecosim_attributes_doubletrawl.png'), 
    height = 1700, width = 2000, res = 200)
opar <- par(mar = c(5, 7, 1, 1))
boxplot(ecosim.outputs, axes = F, cex = 3)
box(lwd = 3)
axis(1, at = axTicks(1), labels = c('Start Bio','End Bio', 'Start Catch', 
                                    'End Catch'), cex.axis = 2)
axis(2, las = T, cex.axis = 2)
mtext(1, text = 'Ecosim statistics', line = 3, cex = 2)
mtext(2, text = '% difference', line = 5 , cex = 2)
dev.off()
