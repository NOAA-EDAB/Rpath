#R Ecosystem
#Comparison to EwE

#User parameters - file locations
if(Sys.info()['sysname']=="Windows"){
  #data.dir <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\data\\"
  #out.dir  <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\outputs\\"
}

if(Sys.info()['sysname']=="Linux"){
  rpath.dir <- "/home/slucey/slucey/Rpath/outputs/"
  ewe.dir   <- "/home/slucey/slucey/Manuscripts/Rpath/EwE_files/"
  out.dir   <- "/home/slucey/slucey/Manuscripts/Rpath/Comparisons/"
}

library(Rpath); library(data.table)

#-------------------------------------------------------------------------------
#Compare Ecopath
Rpath <- as.data.table(read.csv(paste(rpath.dir, 'R_Ecosystem_Parameters.csv', 
                                      sep = '')))
EwE   <- as.data.table(read.csv(paste(ewe.dir,   'REco-Basic estimates.csv',   
                                      sep = '')))

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

png(file = paste(out.dir, 'Ecopath_attributes.png', sep = ''), height = 1700,
    width = 2000, res = 200)
boxplot(ecopath.outputs, axes = F)
box()
axis(1, at = axTicks(1), labels = c('TL','Bio', 'PB', 'QB', 'EE', 'GE'))
axis(2, las = T)
mtext(1, text = 'Ecopath attribute', line = 2.5)
mtext(2, text = '% difference', line = 3)
dev.off()

#-------------------------------------------------------------------------------
#Ecosim
#Base scenario

