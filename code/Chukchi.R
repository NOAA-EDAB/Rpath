#Chukchi Sea test

if(Sys.info()['sysname']=="Windows"){
  data.dir <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\data\\"
  out.dir  <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\outputs\\"
}

if(Sys.info()['sysname']=="Linux"){
  data.dir <- "/home/slucey/slucey/Rpath/data/"
  out.dir  <- "/home/slucey/slucey/Rpath/outputs/"
}

#To download Rpath
#This only needs to be done the first time you run the script
# library(devtools)
# devtools::install_github('slucey/Rpath/Rpath', 
#                          auth_token = 'd95526d2fb3f6e34f9c8481b1740e0033ac1d623')

library(Rpath); library(data.table)

#files
modfile  <- paste(data.dir, 'ECS_eis_base_Jun2015.csv',  sep = '')
dietfile <- paste(data.dir, 'ECS_eis_diet_Jun2015.csv', sep = '')
pedfile  <- paste(data.dir, 'ECS_eis_ped_Jun2015.csv',  sep = '')
juvfile  <- paste(data.dir, 'ECS_eis_juv_Jun2015.csv',  sep = '')

#Run ecosim on the R Ecosystem parameter files
chu <- ecopath(modfile, dietfile, pedfile, 'Chukchi Sea')

webplot(ocean, labels = T, fleets = T, label.cex = 0.65)

#Ecosim
ocean.init <- ecosim.init(ocean, YEARS = 100, juvfile)

#Double trawling effort
ocean.i1  <- copy(ocean.init)
ocean.s1a <- ecosim.run(ocean.i1, 0, 25)
ocean.s1a$fish_Effort[14] <- 2

ocean.s1 <- ecosim.run(ocean.s1a, 25, 100)

write.Rpath.sim(ocean.s1, file = paste(out.dir, 'Ocean_Ecosim_s1.csv', sep = ''))

ecosim.plot(ocean.s1)

ocean.summary <- as.data.table(read.csv(paste(out.dir, 'Ocean_Ecosim_s1.csv', sep = '')))

#Compare with EwE
ewe.2.rpath <- function(x, numliv, numdead){
  ewe <- as.data.table(read.csv(x, skip = 9))
  ewe[, Outside := 0]
  setcolorder(ewe, c('Outside', names(ewe)[which(names(ewe) != 'Outside')]))
  Bmat <- as.matrix(ewe)
  out <- list(spname = names(ewe), out_BB = Bmat, NUM_LIVING = numliv, NUM_DEAD = numdead)
  return(out)
}

ewe       <- ewe.2.rpath(paste(out.dir, 'Ocean_Test_Biomass.csv', sep = ''), numliv = 10, numdead = 1)
ewe.catch <- ewe.2.rpath(paste(out.dir, 'Ocean_Test_Yield.csv', sep = ''), numliv = 10, numdead = 1)

ewe.summary <- data.table(Group = ocean.summary[, Group])
for(i in 2:nrow(ewe.summary)){
  group <- as.character(ocean.summary[i, Group])
  group <- sub(' ', '.', group)
  ewe.summary[i, StartBio := ewe[[2]][1, i]]
  ewe.summary[i, EndBio   := ewe[[2]][nrow(ewe[[2]]) - 1, i]]
  ewe.summary[i, StartCatch := ewe.catch[[2]][1, i]]
  ewe.summary[i, EndCatch   := ewe.catch[[2]][nrow(ewe.catch[[2]]) - 1, i]]
}

ecosim.diff <- merge(ocean.summary, ewe.summary, by = 'Group')
ecosim.diff <- ecosim.diff[!Group %in% c('Outside', 'Detritus')]
ecosim.diff[, StartBio := ((StartBio.x - StartBio.y) / StartBio.y) * 100]
ecosim.diff[, EndBio   := ((EndBio.x   - EndBio.y)   / EndBio.y)   * 100]
ecosim.diff[, StartCatch := (((StartCatch.x * 12) - StartCatch.y) / StartCatch.y) * 100]
ecosim.diff[, EndCatch   := (((EndCatch.x   * 12) - EndCatch.y)   / EndCatch.y)   * 100]
ecosim.diff[, c(paste(c('StartBio', 'EndBio', 'StartCatch', 'EndCatch'), '.x', sep = ''),
                  paste(c('StartBio', 'EndBio', 'StartCatch', 'EndCatch'), '.y', sep = ''),
                  'X') := NULL]

#png(file = paste(out.dir, 'Ecosim_differences_s1.png', sep = ''), height = 800, width = 1063, res = 300)
opar <- par(mar = c(3, 3, 1, 1))
boxplot(ecosim.diff[, 2:5, with = F], axes = F)
axis(1, at = axTicks(1), labels = c('StartBio', 'EndBio', 'StartCatch', 'EndCatch'), cex.axis = 0.7, padj = -1.5)
axis(2, las = T, cex.axis = 0.8, hadj = .7)
box(lwd = 2)
mtext(1, text = 'Ecosim outputs',  line = 1.3)
mtext(2, text = 'Percent difference', line = 2.1)
#dev.off()

