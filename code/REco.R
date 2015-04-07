#R Ecosystem
#Tutorial for using Rpath

#User parameters - file locations
#I have a windows machine and a linux machine, hence the windows toggle
windows <- F
if(windows == T){
  data.dir <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\data\\"
  out.dir  <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\outputs\\"
}
if(windows == F){
  data.dir <- "/home/slucey/slucey/Rpath/data/"
  out.dir  <- "/home/slucey/slucey/Rpath/outputs/"
}

#To download Rpath
#This only needs to be done the first time you run the script
library(devtools)
devtools::install_github('slucey/Rpath/Rpath', 
                          auth_token = 'd95526d2fb3f6e34f9c8481b1740e0033ac1d623')

library(Rpath); library(data.table)

#Identify the location of your parameter files
#For future work you can use create.rpath.param() to generate the parameter 
groups <- c('Whale', 'Fish1', 'Fish2', 'Zooplankton', 'Phytoplankton', 
            'Detritus', 'Discards', 'Fleet1', 'Fleet2')
types <- c(0, 0, 0, 0, 1, 2, 2, 3, 3)
model <- create.rpath.param(group = groups, type = types, parameter = 'model')
model
diet <- create.rpath.param(group = groups, type = types, parameter = 'diet')
diet
juv <- create.rpath.param(group = 'Fish1', parameter = 'juvenile')
juv
ped <- create.rpath.param(group = groups, type = types, parameter = 'pedigree')
ped

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
ecosim.plot(REco.s1)
ecosim.plot(ewe.s1.list)

#Scenario 2 - Increase F on Adult Roundfish 1
REco.i2  <- copy(REco.init)
REco.s2a <- ecosim.run(REco.i2, 0, 25)

REco.s2a$FORCED_FRATE$'AduRoundfish1'[25:100] <- 0.025 
REco.s2b <- ecosim.run(REco.s2a, 25, 100)

#Write out the basic outputs from ecosim
write.Rpath.sim(REco.s2b, file = paste(out.dir, 'R_Ecosystem_Ecosim_s2.csv', sep = ''))

#Plot Relative biomass for sim run
ecosim.plot(REco.s2b)
ecosim.plot(ewe.s2.list)

#Testing
REco.i2  <- copy(REco.init)
REco.s2a <- ecosim.run(REco.i2, 0, 25)

fish <- data.table(Group = REco.s2a$FishFrom,
                   Gear  = REco.s2a$FishThrough,
                   Q     = REco.s2a$FishQ)
fish[Group == 5, Q := c(.169, .007, .019, .105)]
REco.s2a$FishQ <- fish[, Q]

REco.s2b <- ecosim.run(REco.s2a, 25, 100, init_run = F)
ecosim.plot(REco.s2b)
ecosim.plot(ewe.s2.list)

#Scenario 3 - Increase F on Foragefish 1
REco.i3  <- copy(REco.init)
REco.s3a <- ecosim.run(REco.i3, 0, 25)

REco.s3a$FORCED_FRATE$'Foragefish1'[25:100] <- 0.0172 
REco.s3b <- ecosim.run(REco.s3a, 25, 100)

ecosim.plot(REco.s3b)
ecosim.plot(ewe.s3.list)

#Write out the basic outputs from ecosim
write.Rpath.sim(REco.s3b, file = paste(out.dir, 'R_Ecosystem_Ecosim_s3.csv', sep = ''))

REco.i4  <- copy(REco.init)
REco.s4a <- ecosim.run(REco.i4, 0, 25)

fish <- data.table(Group = REco.s4a$FishFrom,
                   Gear  = REco.s4a$FishThrough,
                   Q     = REco.s4a$FishQ)
fish[Group == 13, Q := c(.1764, .000588, .0294)]
REco.s4a$FishQ <- fish[, Q]

REco.s4b <- ecosim.run(REco.s4a, 25, 100, init_run = F)
ecosim.plot(REco.s4b)

