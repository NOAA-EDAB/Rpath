#EMAX
#SML

#User parameters
if(Sys.info()['sysname']=="Windows"){
    data.dir <- "L:\\Rpath\\data\\"
    out.dir  <- "L:\\Rpath\\outputs\\"
    memory.limit(4000)
}
if(Sys.info()['sysname']=="Linux"){
    data.dir <- "/home/slucey/slucey/Rpath/data/"
    out.dir  <- "/home/slucey/slucey/Rpath/outputs/"
}

#-------------------------------------------------------------------------------
#Required packages
library(data.table); library(Rpath)

#-------------------------------------------------------------------------------
#User created functions

#-------------------------------------------------------------------------------
#Georges Bank
modfile  <- paste(data.dir, 'GB_EMAX_mod.csv', sep = '')
dietfile <- paste(data.dir, 'GB_EMAX_diet.csv', sep = '')

GB <- ecopath(modfile, dietfile, 'Georges Bank')

#Webplot plots the resultant food web
png(file = paste(out.dir, 'Georges_Bank_EMAX_Foodweb.png', sep = ''),
    height = 1700, width = 2200, res = 200)
my.groups <- c(c(30, 1, 31), c(2, 3, 4, 7, 10, 13, 5, 9, 15, 11, 12), 
               c(6, 8, 14, 16, 18), c(27, 23, 20, 22, 19, 17),
               c(26, 24, 29, 32), c(28, 25, 21))
webplot(GB, labels = T, box.order = my.groups)
dev.off()

