#Georges Bank Rpath model
#Expanded model of Georges Bank
#SML

#User parameters
if(Sys.info()['sysname']=="Windows"){
  data.dir <- "L:\\Rpath\\data\\"
  out.dir  <- "L:\\Rpath\\outputs\\"
  memory.limit(4000)
}
if(Sys.info()['sysname']=="Linux"){
  data.dir <- "/home/slucey/slucey/Rpath/data"
  out.dir  <- "/home/slucey/slucey/Rpath/outputs"
}

#-------------------------------------------------------------------------------
#Required packages
library(data.table); library(Rpath)

#-------------------------------------------------------------------------------
#User created functions

#-------------------------------------------------------------------------------
#Georges Bank
#This model will have a much less aggregated box structure.
groups <- c('Seabirds', 'Seals', 'Turtles', 'BalWhale', 'ToothWhale', 'Dolphins',
            'HMS', 'Sharks', 'Cod', 'AtlHalibut', 'Bluefish', 'Fourspot', 
            'Goosefish', 'OffHake', 'SilverHake', 'SpinyDogfish', 'Striped Bass',
            'SummerFlounder', 'Weakfish', 'WinterSkate', 'Redfish', 'Pollock', 
            'RedHake','Windowpane', 'SmoothDogfish', 'LittleSkate', 'SmallSkates',
            'AmLobster', 'AmPlaice', 'Barndoor', 'BlackSeaBass', 'Megabenthos', 
            'Red Crab', 'Haddock', 'Ocean Pout', 'Scup', 'WinterFlounder', 
            'WitchFlounder', 'YTFlounder', 'OtherDemersals', 'OtherCephalopods',
            'Illex', 'Loligo', 'WhiteHake', 'AmShad', 'AtlHerring', 'AtlMackerel',
            'AtlMenhaden', 'RiverHerring', 'Butterfish', 'SmallPelagics', 
            'Macrobenthos', 'AtlScallops', 'Clams', 'NShrimp', 'OtherShrimp', 
            'Micronekton', 'GelZooplankton', 'Mesozooplankton', 'Microzooplankton', 
            'Microplankton', 'Nano-picoplankton', 'Bacteria', 'Detritus', 
            'Discards', 'DredgeScallop', 'DredgeClam', 'Gillnet', 'Longline',
            'Seine', 'PotTrap', 'OttertrawlSm', 'OttertrawlLg', 'Midwater', 
            'OtherFisheries')

types <- c(rep(0, 60), 1, 1, 0, 2, 2, rep(3, 10))

GB.params <- create.rpath.params(groups, types)
