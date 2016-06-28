#Generic - 37
#SML

#User parameters
if(Sys.info()['sysname']=="Windows"){
  data.dir <- "L:\\Rpath\\data"
  out.dir  <- "L:\\Rpath\\outputs"
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
groups <- c('Baleen_whales', 'Toothed_whales', 'Seals', 'Birds', 'Sharks_L',
            'Sharks_SM', 'Rays_L', 'Rays_SM', 'Pelagics_L', 'Pelagics_M', 
            'Pelagics_S_carniv', 'Pelagics_S_herbiv', 'Benthopelagics_L', 
            'Benthopelagics_SM', 'Demersals_L', 'Demersals_M', 'Demersals_S', 
            'Reeffish_L', 'Reeffish_M', 'Flatfish_L', 'Flatfish_SM', 
            'Reeffish_S', 'Bathypelagics', 'Bathydemersals', 'Jellyfish', 
            'Cephalopods', 'Shrimps', 'Lobsters_crabs', 'Macrobenthos', 
            'Meiobenthos', 'Corals', 'Softcorals_sponges_etc', 'Krill', 
            'Zooplankton', 'Phytoplankton', 'Benthic plants', 'Detritus', 'Fleet1')

types <- c(rep(0, 30), 0.63, rep(0, 3), rep(1, 2), 2, 3)
#Note - corals are partial primary producers....

G37.params <- create.rpath.param(groups, types)

#Biomass, Production, consumption
biomass <- c(0.001, 0.002, 0.003, 0.001, 0.1, 0.3, 0.1, 0.3, 0.4, 1.2, 5, 2.5,
             0.2, 0.4, 0.5, 2, 5, 0.1, 0.5, 0.1, 1, 1, 0.5, 0.2, 1, 0.5, 2, 0.5,
             10, 15, 0.1, 2, 1, 10, 15, 2, 100, NA)

pb <- c(0.03, 0.05, 0.07, 0.1, 0.3, 0.6, 0.3, 0.6, 0.3, 0.6, 0.6, 0.6, 0.3,
        0.6, 0.3, 0.6, 1, 0.3, 0.6, 0.3, 0.8, 1, 0.5, 0.2, 10, 1, 2.5, 2, 2, 10,
        1, 0.2, 5, 30, 150, 10, NA, NA)

G37.params$model[, Biomass := biomass]
G37.params$model[, PB      := pb]

G37.params$model[Group == 'Baleen_whales',  QB := 30]
G37.params$model[Group == 'Toothed_whales', QB := 40]
G37.params$model[Group == 'Seals',          QB := 50]
G37.params$model[Group == 'Birds',          QB := 100]
G37.params$model[Group == 'Corals',         QB := 1.5]

pq <- c(rep(NA, 4), 0.150, 0.2, 0.15, 0.2, 0.2, rep(0.25, 5), 0.15, 0.2, 0.25,
        0.15, 0.2, 0.15, rep(0.25, 3), 0.3, 0.25, rep(0.3, 5), NA, 0.3, 0.25, 
        0.25, rep(NA, 4))

G37.params$model[, ProdCons := pq]

#Biomass accumulation and unassimilated production
G37.params$model[, BioAcc  := c(rep(0, 37), NA)]
G37.params$model[, Unassim := c(rep(0.2, 33), 0.4, rep(0, 3), NA)]

#Detrital fate
G37.params$model[, Detritus := c(rep(1, 12), 0, 0, rep(1, 22), 0, 0)]

#Landings/Discards
G37.params$model[, Fleet1 := c(rep(0, 4), 0.01, .04, 0.012, 0.02, 0.05, 0.4, 0.1,
                               0.1, 0.05, 0.2, 0.13, 0.8, 0, 0.02, 0.2, 0.025,
                               0.5, rep(0, 4), 0.1, 0.2, 0.4, 0.05, rep(0, 8),
                               NA)]

#Diet
G37.params$diet[, Baleen_whales := c(rep(NA, 9), 0.03, 0.1, 0.1, rep(NA, 20),
                                     0.4, 0.37, rep(NA, 3))]

G37.params$diet[, Toothed_whales := c(NA, NA, 0.001, rep(NA, 5), 0.09, 0.2, 0.309,
                                      0.1, rep(NA, 13), 0.1, rep(NA, 6), 0.2,
                                      rep(NA, 4))]

G37.params$diet[, Seals := c(rep(NA, 9), 0.1, 0.2, 0.1, rep(NA, 3), 0.2, 0.2, 
                             rep(NA, 3), 0.1, rep(NA, 4), 0.1, rep(NA, 11))]

G37.params$diet[, Birds := c(rep(NA, 10), 0.4, 0.2, rep(NA, 4), 0.1, rep(NA, 11),
                             0.1, rep(NA, 3), 0.1, 0.1, rep(NA, 3))]

G37.params$diet[, Sharks_L := c(rep(NA, 5), 0.1, rep(NA, 3), 0.2, 0.1, 0.1,
                                rep(NA, 3), 0.1, 0.1, rep(NA, 4), 0.1, rep(NA, 6),
                                0.2, rep(NA, 8))]

G37.params$diet[, Sharks_SM := c(rep(NA, 10), 0.1, 0.1, rep(NA, 3), 0.1, 0.1, 
                                 rep(NA, 9), 0.1, 0.05, 0.45, rep(NA, 8))]

G37.params$diet[, Rays_L := c(rep(NA, 11), 0.1, rep(NA, 4), 0.1, rep(NA, 9), 0.2,
                              0.1, 0.5, rep(NA, 8))]

G37.params$diet[, Rays_SM := c(rep(NA, 26), 0.1, 0.05, 0.8, rep(NA, 4), 0.05,
                               rep(NA, 3))]

G37.params$diet[, Pelagics_L := c(rep(NA, 9), 0.3, 0.3, 0.1, rep(NA, 3), 0.1, 0.1,
                                  rep(NA, 15), 0.1, rep(NA, 4))]

G37.params$diet[, Pelagics_M := c(rep(NA, 10), 0.5, 0.3, rep(NA, 13), 0.1,
                                  rep(NA, 6), 0.1, rep(NA, 4))]
                
G37.params$diet[, Pelagics_S_carniv := c(rep(NA, 32), 0.1, 0.9, rep(NA, 3))]

G37.params$diet[, Pelagics_S_herbiv := c(rep(NA, 33), 0.222, 0.778, NA, NA)]

G37.params$diet[, Benthopelagics_L := c(rep(NA, 10), 0.3, 0.1, rep(NA, 4), 0.1, 
                                        rep(NA, 9), 0.2, NA, 0.1, rep(NA, 3),
                                        0.1, 0.1, rep(NA, 3))] 

G37.params$diet[, Benthopelagics_SM := c(rep(NA, 10), 0.2, 0.1, rep(NA, 4), 0.2,
                                         rep(NA, 9), 0.2, NA, 0.2, rep(NA, 4), 
                                         0.1, rep(NA, 3))]

G37.params$diet[, Demersals_L := c(rep(NA, 10), 0.05, 0.05, rep(NA, 3), 0.1, 0.5,
                                   rep(NA, 9), 0.2, NA, 0.1, rep(NA, 8))]

G37.params$diet[, Demersals_M := c(rep(NA, 16), 0.6, rep(NA, 9), 0.2, 0.05, 0.15,
                                   rep(NA, 8))]

G37.params$diet[, Demersals_S := c(rep(NA, 26), 0.1, NA, 0.5, 0.2, rep(NA, 3),
                                   0.2, rep(NA, 3))]

G37.params$diet[, Reeffish_L := c(rep(NA, 10), 0.05, 0.05, rep(NA, 4), 0.1, NA, 
                                  0.2, NA, NA, 0.4, rep(NA, 6), 0.2, rep(NA, 8))]

G37.params$diet[, Reeffish_M := c(rep(NA, 16), 0.1, rep(NA, 4), 0.4, rep(NA, 4),
                                  0.1, NA, 0.15, 0.1, 0.03, 0.02, NA, 0.1, 
                                  rep(NA, 3))] 

G37.params$diet[, Flatfish_L := c(rep(NA, 15), 0.05, 0.25, rep(NA, 3), 0.2, 
                                  rep(NA, 4), 0.1, 0.1, NA, 0.2, rep(NA, 4), 0.1,
                                  rep(NA, 3))]

G37.params$diet[, Flatfish_SM := c(rep(NA, 16), 0.05, rep(NA, 9), 0.1, NA, 0.6,
                                   0.1, rep(NA, 3), 0.15, rep(NA, 3))]

G37.params$diet[, Reeffish_S := c(rep(NA, 26), 0.1, rep(NA, 5), 0.1, 0.8, 
                                  rep(NA, 3))]

G37.params$diet[, Bathypelagics := c(rep(NA, 22), 0.05, rep(NA, 9), 0.2, 0.65, 
                                     NA, NA, 0.1)]

G37.params$diet[, Bathydemersals := c(rep(NA, 23), 0.05, rep(NA, 4), 0.5, 0.2,
                                      rep(NA, 3), 0.05, NA, NA, 0.2)]

G37.params$diet[, Jellyfish := c(rep(NA, 24), 0.02, rep(NA, 8), 0.8, NA, NA, 0.18)]

G37.params$diet[, Cephalopods := c(rep(NA, 10), 0.01, 0.01, rep(NA, 20), 0.1, 
                                   0.78, NA, NA, 0.1)]

G37.params$diet[, Shrimps := c(rep(NA, 29), 0.1, rep(NA, 6), 0.9)]

G37.params$diet[, Lobsters_crabs := c(rep(NA, 28), 0.2, 0.5, rep(NA, 6), 0.3)]

G37.params$diet[, Macrobenthos := c(rep(NA, 28), 0.02, 0.6, rep(NA, 3), 0.02, NA,
                                    0.01, 0.35)]

G37.params$diet[, Meiobenthos := c(rep(NA, 29), 0.1, rep(NA, 3), 0.4, 0.2, NA, 
                                   0.3)]

#G37.params$diet[, Corals := c(rep(NA, 33), 0.5, 0.1, NA, 0.4)]
G37.params$diet[, Corals := c(rep(NA, 33), 0.185, 0.037, NA, 0.148)]

G37.params$diet[, Softcorals_sponges_etc := c(rep(NA, 33), 0.5, 0.1, NA, 0.4)]

G37.params$diet[, Krill := c(rep(NA, 32), 0.05, 0.7, 0.1, NA, 0.15)]

G37.params$diet[, Zooplankton := c(rep(NA, 34), 0.9, NA, 0.1)]

check.rpath.param(G37.params)

save(G37.params, file = file.path(data.dir, "Generic_37_params.RData"))

G37.params$model[Group == 'Demersals_S', BioAcc := -.2]
#Ecopath
G37 <- rpath(G37.params, 'Generic 37')

webplot(G37, labels = T)

#Ecosim
G37.sim <- rsim.scenario(G37, G37.params, 10)
G37.run1 <- rsim.run(G37.sim, method = 'AB', years = 10)
rsim.plot(G37.run1, groups[1:37])
write.Rsim(REco.run1, file = paste(out.dir, 'REco_baserun_AB.csv', sep = ''))

G37.sim <- rsim.scenario(G37, G37.params, 9)
G37.sim <- adjust.fishing(G37.sim, 'EFFORT', gear = 'Fleet1', 
                           year = 0:9, value = c(1, 2.279, 2.7, 2.846, 3.084,
                                                 3.084, 3.084, 3.03, 2.993, 2.993))
G37.run2 <- rsim.run(G37.sim, method = 'AB', years = 9)
rsim.plot(G37.run2, groups[1:37])


