## ----rpath load, echo = F-----------------------------------------------------
knitr::opts_chunk$set(
  comment = '#>',
  collapse = T)
library(Rpath); library(data.table)

## ----creating Rpath object, echo = F------------------------------------------
#Groups and types for the R Ecosystem

groups <- c('Seabirds', 'Whales', 'Seals', 'JuvRoundfish1', 'AduRoundfish1', 
            'JuvRoundfish2', 'AduRoundfish2', 'JuvFlatfish1', 'AduFlatfish1',
            'JuvFlatfish2', 'AduFlatfish2', 'OtherGroundfish', 'Foragefish1',
            'Foragefish2', 'OtherForagefish', 'Megabenthos', 'Shellfish',
            'Macrobenthos', 'Zooplankton', 'Phytoplankton', 'Detritus', 
            'Discards', 'Trawlers', 'Midwater', 'Dredgers')

types  <- c(rep(0, 19), 1, rep(2, 2), rep(3, 3))

stgroups <- c(rep(NA, 3), rep('Roundfish1', 2), rep('Roundfish2', 2), 
              rep('Flatfish1', 2), rep('Flatfish2', 2), rep(NA, 14))

REco.params <- create.rpath.params(group = groups, type = types, stgroup = stgroups)

#Model
biomass <- c(0.0149, 0.454, NA, NA, 1.39, NA, 5.553, NA, 5.766, NA,
             0.739, 7.4, 5.1, 4.7, 5.1, NA, 7, 17.4, 23, 10, rep(NA, 5))

pb <- c(0.098, 0.031, 0.100, 2.026, 0.42, 2.1, 0.425, 1.5, 0.26, 1.1, 0.18, 0.6,
        0.61, 0.65, 1.5, 0.9, 1.3, 7, 39, 240, rep(NA, 5))

qb <- c(76.750, 6.976, 34.455, NA, 2.19, NA, 3.78, NA, 1.44, NA, 1.69,
        1.764, 3.52, 5.65, 3.6, 2.984, rep (NA, 9))

REco.params$model[, Biomass := biomass]
REco.params$model[, PB := pb]
REco.params$model[, QB := qb]

#EE for groups w/o biomass
REco.params$model[Group %in% c('Seals', 'Megabenthos'), EE := 0.8]

#Production to Consumption for those groups without a QB
REco.params$model[Group %in% c('Shellfish', 'Zooplankton'), ProdCons:= 0.25]
REco.params$model[Group == 'Macrobenthos', ProdCons := 0.35]

#Biomass accumulation and unassimilated production
REco.params$model[, BioAcc  := c(rep(0, 22), rep(NA, 3))]
REco.params$model[, Unassim := c(rep(0.2, 18), 0.4, rep(0, 3), rep(NA, 3))]

#Detrital Fate
REco.params$model[, Detritus := c(rep(1, 20), rep(0, 5))]
REco.params$model[, Discards := c(rep(0, 22), rep(1, 3))]

#Fisheries
#Landings
trawl  <- c(rep(0, 4), 0.08, 0, 0.32, 0, 0.09, 0, 0.05, 0.2, rep(0, 10), rep(NA, 3))
mid    <- c(rep(0, 12), 0.3, 0.08, 0.02, rep(0, 7), rep(NA, 3))
dredge <- c(rep(0, 15), 0.1, 0.5, rep(0, 5), rep(NA, 3))
REco.params$model[, Trawlers := trawl]
REco.params$model[, Midwater := mid]
REco.params$model[, Dredgers := dredge]

#Discards
trawl.d  <- c(1e-5, 1e-7, 0.001, 0.001, 0.005, 0.001, 0.009, 0.001, 0.04, 0.001,
              0.01, 0.08, 0.001, 0.001, 0.001, rep(0, 7), rep(NA, 3))
mid.d    <- c(rep(0, 2), 0.001, 0.001, 0.01, 0.001, 0.01, rep(0, 4), 0.05, 0.05,
              0.01, 0.01, rep(0, 7), rep(NA, 3))
dredge.d <- c(rep(0, 3), 0.001, 0.05, 0.001, 0.05, 0.001, 0.05, 0.001, 0.01, 0.05,
              rep(0, 3), 0.09, 0.01, 1e-4, rep(0, 4), rep(NA, 3))
REco.params$model[, Trawlers.disc := trawl.d]
REco.params$model[, Midwater.disc := mid.d]
REco.params$model[, Dredgers.disc := dredge.d]

#Group parameters
REco.params$stanzas$stgroups[, VBGF_Ksp := c(0.145, 0.295, 0.0761, 0.112)]
REco.params$stanzas$stgroups[, Wmat     := c(0.0769, 0.561, 0.117,  0.321)]

#Individual stanza parameters
REco.params$stanzas$stindiv[, First   := c(rep(c(0, 24), 3), 0, 48)]
REco.params$stanzas$stindiv[, Last    := c(rep(c(23, 400), 3), 47, 400)]
REco.params$stanzas$stindiv[, Z       := c(2.026, 0.42, 2.1, 0.425, 1.5, 
                                           0.26, 1.1, 0.18)]
REco.params$stanzas$stindiv[, Leading := rep(c(F, T), 4)]

REco.params <- rpath.stanzas(REco.params)

#Diets
REco.params$diet[, Seabirds        := c(rep(NA, 11), 0.1, 0.25, 0.2, 0.15, 
                                         rep(NA, 6), 0.3, NA)]
REco.params$diet[, Whales          := c(rep(NA, 3), 0.01, NA, 0.01, NA, 0.01, 
                                         NA, 0.01, rep(NA, 4), 0.1, rep(NA, 3), 
                                         0.86, rep(NA, 4))]
REco.params$diet[, Seals           := c(rep(NA, 3), 0.05, 0.1, 0.05, 0.2, 0.005, 
                                         0.05, 0.005, 0.01, 0.24, rep(0.05, 4), 
                                         0.09, rep(NA, 6))]
REco.params$diet[, JuvRoundfish1   := c(rep(NA, 3), rep(c(1e-4, NA), 4), 1e-3, 
                                         rep(NA, 2), 0.05, 1e-4, NA, .02, 0.7785, 
                                         0.1, 0.05, NA, NA)]
REco.params$diet[, AduRoundfish1   := c(rep(NA, 5), 1e-3, 0.01, 1e-3, 0.05, 1e-3, 
                                         0.01, 0.29, 0.1, 0.1, 0.347, 0.03, NA, 
                                         0.05, 0.01, rep(NA, 4))]
REco.params$diet[, JuvRoundfish2   := c(rep(NA, 3), rep(c(1e-4, NA), 4), 1e-3, 
                                         rep(NA, 2), 0.05, 1e-4, NA, .02, 0.7785, 
                                         0.1, .05, NA, NA)]
REco.params$diet[, AduRoundfish2   := c(rep(NA, 3), 1e-4, NA, 1e-4, NA, rep(1e-4, 4), 
                                         0.1, rep(0.05, 3), 0.2684, 0.01, 0.37, 0.001, 
                                         NA, 0.1, NA, NA)]
REco.params$diet[, JuvFlatfish1    := c(rep(NA, 3), rep(c(1e-4, NA), 4), rep(NA, 3), 
                                         rep(1e-4, 2), NA, 0.416, 0.4334, 0.1, 0.05, 
                                         NA, NA)]
REco.params$diet[, AduFlatfish1    := c(rep(NA, 7), rep(1e-4, 5), rep(NA, 2), 0.001, 
                                         0.05, 0.001, 0.6, 0.2475, NA, 0.1, NA, NA)]
REco.params$diet[, JuvFlatfish2    := c(rep(NA, 3), rep(c(1e-4, NA), 4), rep(NA, 3),
                                         rep(1e-4, 2), NA, 0.416, 0.4334, 0.1, 0.05, 
                                         NA, NA)]
REco.params$diet[, AduFlatfish2    := c(rep(NA, 7), 1e-4, NA, 1e-4, rep(NA, 4), 
                                         rep(1e-4, 3), 0.44, 0.3895, NA, 0.17, NA, NA)]
REco.params$diet[, OtherGroundfish := c(rep(NA, 3), rep(1e-4, 8), 0.05, 0.08, 0.0992, 
                                         0.3, 0.15, 0.01, 0.3, 0.01, rep(NA, 4))]
REco.params$diet[, Foragefish1     := c(rep(NA, 3), rep(c(1e-4, NA), 4), rep(NA, 7), 
                                         0.8196, 0.06, 0.12, NA, NA)]
REco.params$diet[, Foragefish2     := c(rep(NA, 3), rep(c(1e-4, NA), 4), rep(NA, 7), 
                                         0.8196, 0.06, 0.12, NA, NA)]
REco.params$diet[, OtherForagefish := c(rep(NA, 3), rep(c(1e-4, NA), 4), rep(NA, 7), 
                                         0.8196, 0.06, 0.12, NA, NA)]
REco.params$diet[, Megabenthos     := c(rep(NA, 15), 0.1, 0.03, 0.55, rep(NA, 2), 0.32,
                                         NA, NA)]
REco.params$diet[, Shellfish       := c(rep(NA, 18), 0.3, 0.5, 0.2, NA, NA)]
REco.params$diet[, Macrobenthos    := c(rep(NA, 16), 0.01, rep(0.2, 2), NA, 0.59, NA, NA)]
REco.params$diet[, Zooplankton     := c(rep(NA, 18), 0.2, 0.6, 0.2, NA, NA)]

REco <- rpath(REco.params, eco.name = 'R Ecosystem')

## ----rsim.scenario------------------------------------------------------------
REco.sim <- rsim.scenario(REco, REco.params, years = 1:100)

## ----Change effort------------------------------------------------------------
REco.sim <- adjust.fishing(REco.sim, 'ForcedEffort', group = 'Trawlers', sim.year = 25:100, 
                           value = 2)

## ----rsim.run-----------------------------------------------------------------
REco.run1 <- rsim.run(REco.sim, method = 'RK4', years = 1:100)

## ----ecosim equilibrium plot, fig.height = 7, fig.width = 9-------------------
rsim.plot(REco.run1, groups[1:22])

