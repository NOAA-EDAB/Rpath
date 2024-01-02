#Anchovy Bay
#SML

#Required packages
library(data.table); library(Rpath); library(here)

#Anchovy Bay----
groups <- c('whales', 'seals', 'cod', 'whiting', 'mackerel', 'anchovy', 'shrimp',
            'benthos', 'zooplankton', 'phytoplankton', 'detritus', 'sealers', 
            'trawlers', 'seiners', 'bait boats', 'shrimpers')

types <- c(rep(0, 9), 1, 2, rep(3, 5))

AB.params <- create.rpath.params(groups, types)

#Biomass, Production, consumption
biomass <- c(0.08, 0.0609, 3, 1.8, 1.2, 7, 0.8, NA, 14.8, 9, 10, rep(NA, 5))

pb <- c(0.05, 0.164, 0.340, 0.581, 0.723, 1.140, 3, 3, 35, 240, rep(NA, 6))

qb <- c(9, 15, 2.58, 3.3, 4.4, 9.13, rep(NA, 10))

AB.params$model[, Biomass := biomass]
AB.params$model[, PB      := pb]
AB.params$model[, QB      := qb]

AB.params$model[Group == 'shrimp',      ProdCons := 0.25]
AB.params$model[Group == 'benthos',     ProdCons := 0.25]
AB.params$model[Group == 'zooplankton', ProdCons := 0.25]

#Add EE's for unknown biomasses
AB.params$model[Group == 'benthos', EE := 0.6]
#Biomass accumulation and unassimilated production
AB.params$model[, BioAcc  := c(rep(0, 11), rep(NA, 5))]
AB.params$model[, Unassim := c(rep(0.2, 9), rep(0, 2), rep(NA, 5))]

#Detrital fate
AB.params$model[, detritus := c(rep(1, 10), rep(0, 6))]

#Landings/Discards
AB.params$model[Group == 'seals',    sealers      := .0045]
AB.params$model[Group == 'cod',      trawlers     := 0.45]
AB.params$model[Group == 'whiting',  trawlers     := 0.2]
AB.params$model[Group == 'mackerel', seiners      := 0.4]
AB.params$model[Group == 'anchovy',  seiners      := 1.2]
AB.params$model[Group == 'anchovy',  "bait boats" := 0.2]
AB.params$model[Group == 'shrimp',   shrimpers    := 0.05]

#Diet
AB.params$diet[, whales      := c(rep(NA, 2), 0.1, 0.1, 0.2, 0.5, NA, 0.1, rep(NA, 4))]
AB.params$diet[, seals       := c(NA, NA, 0.04, 0.05, NA, NA, 0.01, 0.9, rep(NA, 4))]
AB.params$diet[, cod         := c(NA, NA, NA, 0.05, NA, 0.1, 0.01, 0.84, rep(NA, 4))]
AB.params$diet[, whiting     := c(NA, NA, 0.05, 0.05, NA, 0.45, 0.01, 0.44, rep(NA, 4))]
AB.params$diet[, mackerel    := c(rep(NA, 4), 0.05, 0.5, NA, NA, 0.45, rep(NA, 3))]
AB.params$diet[, anchovy     := c(rep(NA, 8), 1, rep(NA, 3))]
AB.params$diet[, shrimp      := c(rep(NA, 7), 1, rep(NA, 4))]
AB.params$diet[, benthos     := c(rep(NA, 7), 0.1, 0.1, 0.1, 0.7, NA)]
AB.params$diet[, zooplankton := c(rep(NA, 9), 0.9, 0.1, NA)]

usethis::use_data(AB.params, overwrite = T)
