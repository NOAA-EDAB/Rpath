#R Ecosystem
#Tutorial for using Rpath

#User parameters - file locations
#I have a windows machine and a linux machine, hence the windows toggle
windows <- F
if(windows == T){
  r.dir    <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\code\\"
  data.dir <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\data\\"
  out.dir  <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\outputs\\"
  stable   <- "L:\\PhD\\Rpath\\code\\"
}
if(windows == F){
  r.dir    <- "/home/slucey/slucey/Rpath/code/"
  data.dir <- "/home/slucey/slucey/Rpath/data/"
  out.dir  <- "/home/slucey/slucey/Rpath/outputs/"
  stable   <- "/home/slucey/slucey/PhD/Rpath/code/"
}

#To download Rpath
#This only needs to be done the first time you run the script
library(devtools)
devtools::install_github('slucey/Rpath/Rpath', 
                          auth_token = 'd95526d2fb3f6e34f9c8481b1740e0033ac1d623')

library(Rpath); library(data.table)

#Identify the location of your parameter files
#For future work you can use create.rpath.param() to generate the parameter 
#file skeletons
modfile  <- paste(data.dir, 'REco_mod.csv',  sep = '')
dietfile <- paste(data.dir, 'REco_diet.csv', sep = '')
pedfile  <- paste(data.dir, 'REco_ped.csv',  sep = '')
juvfile  <- paste(data.dir, 'REco_juv.csv',  sep = '')

#Run ecosim on the R Ecosystem parameter files
REco <- ecopath(modfile, dietfile, pedfile, 'R Ecosystem')

#There is an rpath method for generic functions print and summary
REco
summary(REco)

#print also includes the mortalities toggle
print(REco, morts = T)

#Summary table of parameters or mortalities can be output
write.Rpath(REco, file = paste(out.dir, 'R_Ecosystem_Parameters.csv', sep = ''))
write.Rpath(REco, file = paste(out.dir, 'R_Ecosystem_Mortalities.csv', sep = ''), morts = T)

#Webplot plots the resultant food web
set.seed(34)
#Set the plot order to avoid overlapping titles/points (optional)
my.order <- c(c(22, 20, 21), c(17, 19, 18), c(8, 4, 2, 14, 11, 13, 10, 15, 16, 9, 6),
              c(25, 1, 7, 12), c(3, 5, 24), 23)

tiff(file = paste(out.dir, 'R_Ecosystem.tif'), height = 1500, width = 1700, res = 300)
webplot(REco, labels = T, fleets = T, box.order = my.order, label.cex = 0.65)
dev.off()

#Highlight example
set.seed(34)
tiff(file = paste(out.dir, 'Highlight_AduRoundfish1.tif'), height = 1800, width = 2000, res = 300)
webplot(REco, fleets = T, highlight = 5, box.order = my.order, label.cex = 0.65)
dev.off()



#Use the Rpath object from above to run a 100 year Ecosim
#Scenario 1 - Equilibrium
REco.init <- ecosim.init(REco, YEARS = 100, juvfile)
REco.i1   <- copy(REco.init)
REco.s1   <- ecosim.run(REco.i1, 0, 100)

#Write out the basic outputs from ecosim
write.Rpath.sim(REco.s1, file = paste(out.dir, 'R_Ecosystem_Ecosim_s1.csv', sep = ''))

#Plot Relative biomass for sim run
biomass <- as.data.table(REco.s1$out_BB[, 2:23])
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
tiff(file = paste(out.dir, 'REco_relative_bio_Base_run.tif', sep = ''), 
     height = 1800, width = 2000, res = 300)
layout(matrix(c(1, 2), 1, 2), widths = c(4, 1))
opar <- par(mar = c(4, 6, 2, 0))
plot(0, 0, ylim = c(ymin, ymax), xlim = c(0, xmax), 
     axes = F, xlab = '', ylab = '', type = 'n')
axis(1)
axis(2, las = T)
box(lwd = 2)
mtext(1, text = 'Months', line = 2.5, cex = 1.8)
mtext(2, text = 'Relative Biomass', line = 3, cex = 1.8)

line.col <- rainbow(22)
for(i in 1:(ncol(rel.bio) - 1)){
  lines(rel.bio[, c(1, i + 1), with = F], col = line.col[i], lwd = 3)
}

opar <- par(mar = c(0, 0, 0, 0))
plot(0, 0, xlab = '', ylab = '', axes = F)
legend('center', legend = REco.s1$spname[2:23], fill = line.col, cex = .6)
dev.off()


#Scenario 2 - Increase F on Adult Roundfish 1
REco.i2  <- copy(REco.init)
REco.s2a <- ecosim.run(REco.i2, 0, 24)

REco.s2a$FORCED_FRATE$'AduRoundfish1'[24:100] <- 0.025 
REco.s2b <- ecosim.run(REco.s2a, 24, 100)

#Write out the basic outputs from ecosim
write.Rpath.sim(REco.s2b, file = paste(out.dir, 'R_Ecosystem_Ecosim_s2.csv', sep = ''))

#Plot Relative biomass for sim run
biomass <- as.data.table(REco.s2b$out_BB[, 2:23])
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
tiff(file = paste(out.dir, 'REco_relative_bio_scenario_2.tif', sep = ''), 
     height = 1800, width = 2000, res = 300)
layout(matrix(c(1, 2), 1, 2), widths = c(4, 1))
opar <- par(mar = c(4, 6, 2, 0))
plot(0, 0, ylim = c(ymin, ymax), xlim = c(0, xmax), 
     axes = F, xlab = '', ylab = '', type = 'n')
axis(1)
axis(2, las = T)
box(lwd = 2)
mtext(1, text = 'Months', line = 2.5, cex = 1.8)
mtext(2, text = 'Relative Biomass', line = 3, cex = 1.8)

line.col <- rainbow(22)
for(i in 1:(ncol(rel.bio) - 1)){
  lines(rel.bio[, c(1, i + 1), with = F], col = line.col[i], lwd = 3)
}

opar <- par(mar = c(0, 0, 0, 0))
plot(0, 0, xlab = '', ylab = '', axes = F)
legend('center', legend = REco.s1$spname[2:23], fill = line.col, cex = .6)
dev.off()




