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
groups <- c('Phytoplankton- Primary Producers', 'Bacteria', 'Microzooplankton',
            'Small copepods', 'Large Copepods', 'Gelatinous Zooplankton', 
            'Micronekton', 'Mesopelagics', 'Macrobenthos- polychaetes', 
            'Macrobenthos- crustaceans', 'Macrobenthos- molluscs', 
            'Macrobenthos- other', 'Megabenthos- filterers', 'Megabenthos- other',
            'Shrimp et al.', 'Larval-juv fish- all', 'Small Pelagics- commercial',
            'Small Pelagics- other', 'Small Pelagics- squid', 
            'Small Pelagics- anadromous', 'Medium Pelagics- (piscivores & other)',
            'Demersals- benthivores', 'Demersals- omnivores', 
            'Demersals- piscivores', 'Sharks- pelagics', 'HMS', 'Baleen Whales',
            'Odontocetes', 'Sea Birds', 'Discard', 'Detritus-POC', 'Fishery')
types <- c(1, rep(0, 28), rep(2, 2), 3)

GB.params <- create.rpath.param(groups, types)

GB.params$model$Biomass <- c(25.70472, 6.517908, 5.587981, 12.98514, 6.980794,
                             1.319463, 3.805126, 0.045, 11.40272, 10.87353,
                             9.8865, 40.02257, 3.613779, 3.965064, 0.09, 
                             0.6293996, 14.97737, 1.0737, 1.262162, 0.25,
                             0.2915301, 4.576176, 3.438957, 2.244675, 0.04334066,
                             0.00439539, 0.4167178, 0.1127281, 0.003496863,
                             0.4836398, 50, NA)

GB.params$model$PB <- c(166.1342, 91.24998, 72.00002, 41.66504, 54.63586, 
                        40, 14.25, 0.9503762, 2.5, 3, 2, 2, 5, 2, 2, 15, 
                        0.3452712, 0.9571092, 0.9503762, 0.4249809, 0.459,
                        0.45, 0.45, 0.486, 0.102, 0.682623, 0.03802086,
                        0.04, 0.275, NA, NA, NA)

GB.params$model$QB <- c(NA, 380.2082, 242.4243, 127.75, 109.5, 143.08, 36.5,
                        1.825, 17.5, 21, 14, 17.64, 18, 18, 5, 45, 2, 2,
                        2.75, 2, 2.3814, 0.92, 0.83, 2.205567, 0.5328019,
                        2.053014, 4.5, 13.82976, 4.379231, NA, NA, NA)

GB.params$model$BioAcc <- c(rep(0, 29), rep(NA, 3))

GB.params$model$Unassim <- c(0, 0.2, 0.1, 0.25, 0.25, 0.35, 0.25, 0.15, 0.5, 
                             0.5, 0.6, 0.5, 0.7, 0.3, 0.3, 0.15, 0.15,  0.35, 
                             0.15, 0.15, 0.15, 0.3, 0.35, 0.15, 0.15, 0.15, 0.2, 
                             0.2, 0.15, NA, NA, NA)

GB.params$model$Discard <- c(rep(0.01, 7), 0, rep(0.01, 23), NA)

GB.params$model$'Detritus-POC' <- c(rep(0.99, 7), 0, rep(0.99, 21), 0.001, 0, NA)

GB.params$model$Fishery <- c(rep(0, 11), 0.6478018, 0.03305462, 0, 0.00022901,
                             0.2990043, 0, 0.003389353, 0.03275818, 0.01015249,
                             0.1, 0.005314156, 0.5314448, 0, 0.003406343, 
                             rep(0, 6), NA)

GB.params$model$Fishery.disc <- c(rep(0, 5), 6.36e-7, 0, 3.7e-12, rep(0, 4), .2,
                                 0, 0, 4.8e-9, 0.0299, 0, 3.39e-4, 3.28e-3, 
                                 3.05e-3, 7.92e-2, 1.59e-3, .159, 2.38e-4, 
                                 1.02e-3, 1.25e-8, 1.02e-8, 1.44e-4, 0, 0, NA)

#Diet
GB.params$diet[, Bacteria := c(0.24, rep(NA, 29), 0.76)]

GB.params$diet[, Microzooplankton := c(0.216, 0.16, 0.12, rep(NA, 27), 0.504)]

GB.params$diet[, 'Small copepods' := c(0.724, NA, 0.08, 0.065, rep(NA, 26), 0.131)]

GB.params$diet[, 'Large Copepods' := c(0.546, NA, 0.0439, 0.174, 0.122, 0.0531, 
                                       rep(NA, 3), 1.92e-4, NA, 1.02e-4, 
                                       rep(NA, 18), 0.0594)]

GB.params$diet[, 'Gelatinous Zooplankton' := 
                 c(0.087, 0.02, 0.051, 0.335, 0.366, 0.021, rep(NA, 9), 0.01, 
                   0.005, 0.002, 0.00043, 0.000016, rep(NA, 10), 0.102)]

GB.params$diet[, Micronekton := c(0.162, NA, NA, .308, .325, NA, .041, rep(NA, 23),
                                   0.163)]

GB.params$diet[, Mesopelagics := c(0.028, 0.017, 0.072, 0.34, 0.524, NA, 0.014, 
                                   rep(NA, 23), 0.004)]

GB.params$diet[, 'Macrobenthos- polychaetes' := 
                 c(0.128, 0.308, rep(NA, 6), 0.015, 9.8e-4, 5.79e-4, 3.54e-3, 
                   3.56e-3, 1.9e-4, rep(NA, 15), 0.00593, 0.534)]

GB.params$diet[, 'Macrobenthos- crustaceans' :=
                 c(0.212, 0.187, NA, 0.019, 0.037, rep(NA, 3), 0.015, 0.008, 
                   0.008, 0.026, 0.018, 2.9e-4, rep(NA, 7), 8.5e-4, 5.1e-4, 
                   rep(NA, 6), 0.009, 0.459)]

GB.params$diet[, 'Macrobenthos- molluscs' :=
                 c(0.432, 0.199, rep(NA, 8), 0.004, 0.003, 0.013, 4.3e-4, 
                   rep(NA, 15), 0.006, 0.342)]

GB.params$diet[, 'Macrobenthos- other' := 
                 c(0.215, 0.231, rep(NA, 6), 0.017, 0.025, 0.016, 0.057, 0.006,
                   0.003, rep(NA, 7), 8.5e-4, 4.2e-4, 2.4e-7, rep(NA, 5), 0.01,
                   0.418)]

GB.params$diet[, 'Megabenthos- filterers' := c(0.69, 0.08, rep(NA, 28), 0.23)]

GB.params$diet[, 'Megabenthos- other' := c(NA, 0.176, rep(NA, 6), 0.078, 0.098,
                                           0.032, 0.294, 0.048, 0.045, rep(NA, 7),
                                           0.003, 0.002, 2.4e-7, rep(NA, 5), 0.048,
                                           0.176)]

GB.params$diet[, 'Shrimp et al.' := c(0.062, 0.365, rep(NA, 4), 0.123, NA, NA, 
                                      0.009, NA, 0.013, NA, NA, 5.4e-4, 
                                      rep(NA, 14), 0.062, 0.365)]

GB.params$diet[, 'Larval-juv fish- all' := c(0.062, NA, NA, 0.456, 0.264, NA, 
                                             0.072, NA, 0.01, 0.007, 0.005, 0.005,
                                             rep(NA, 3), 0.055, rep(NA, 14), 0.062)]

GB.params$diet[, 'Small Pelagics- commercial' :=
                 c(0.0113, NA, NA, 0.15, 0.436, 0.0818, 0.165, NA, 0.00892, 
                   0.0247, 0.00892, 0.0113, rep(NA, 3), 0.0963, rep(NA, 5), 
                   9.06e-4, 9.06e-4, 0.00413, rep(NA, 7))]

GB.params$diet[, 'Small Pelagics- other' :=
                 c(0.158, NA, NA, 0.115, 0.573, 0.102, 0.042, NA, NA, 6.9e-4, 
                   4.8e-4, 2.4e-4, rep(NA, 3), 0.007, rep(NA, 14), 0.0007)]

GB.params$diet[, 'Small Pelagics- squid' := 
                 c(rep(NA, 4), 0.129, NA, 0.456, NA, NA, 0.099, NA, 0.018, NA, 
                   NA, 0.011, 0.176, 0.016, 0.018, 0.077, 1e-4, rep(NA, 11))]

GB.params$diet[, 'Small Pelagics- anadromous' := 
                 c(0.012, NA, NA, 0.056, 0.9, NA, 0.02, NA, 4.9e-4, 0.002, 
                   rep(NA, 5), 0.008, rep(NA, 14), 0.001)]

GB.params$diet[, 'Medium Pelagics- (piscivores & other)' :=
                 c(rep(NA, 5), 0.001, NA, 0.003, NA, 0.013, NA, 0.011, 0.003, 
                   0.018, 0.001, 0.002, 0.576, 0.044, 0.114, 0.01, 0.011, 0.098,
                   0.014, 0.079, rep(NA, 6), 0.001)]

GB.params$diet[, 'Demersals- benthivores' :=
                 c(rep(NA, 5), 0.005, 0.001, NA, 0.111, 0.13, 0.111, 0.133, 0.104,
                   0.133, 0.006, NA, 0.11, 0.001, 0.01, NA, NA, 0.076, 0.029,
                   0.018, rep(NA, 5), 0.011, 0.011)]
               
GB.params$diet[, 'Demersals- omnivores' :=
                 c(rep(NA, 5), 0.006, 0.028, NA, 0.132, 0.066, 0.066, 0.066, 0.072,
                   0.286, 0.01, 0.013, 0.132, 0.002, 0.022, NA, NA, 0.048, 0.004,
                   0.025, rep(NA, 5), 0.011, 0.011)]

GB.params$diet[, 'Demersals- piscivores' :=
                 c(rep(NA, 5), 0.026, 0.001, 0.007, 0.013, 0.013, 0.02, 0.118,
                   0.246, 0.019, 0.012, 0.013, 0.322, 0.032, 0.015, 0.009, 0.001,
                   0.038, 0.007, 0.086, rep(NA, 6), 0.001)]

GB.params$diet[, 'Sharks- pelagics' :=
                 c(rep(NA, 4), 0.031, 0.01, NA, 0.007, NA, 0.01, NA, 0.01, 
                   rep(NA, 4), 0.216, 0.082, 0.165, 0.001, 0.124, 0.051, 0.082,
                   0.072, 0.014, 0.01, 0.01, 0.021, 0.031, NA, 0.051)]

GB.params$diet[, HMS := c(rep(NA, 5), 0.103, rep(NA, 10), 0.135, 0.741, 0.022, 
                          rep(NA, 12))]

GB.params$diet[, 'Baleen Whales' := c(rep(NA, 3), 0.059, 0.473, 0.001, 0.296, NA, 
                                      NA, 0.059, 0.006, 0.024, NA, 0.004, NA, NA,
                                      0.059, 0.004, 0.003, 1.3e-4, rep(NA, 10),
                                      0.012)]

GB.params$diet[, Odontocetes := c(rep(NA, 5), 0.003, 0.027, rep(NA, 9), 0.379,
                                  0.205, 0.274, 0.002, 3.2e-4, NA, 0.068, 0.04,
                                  rep(NA, 3), 0.002, rep(NA, 3))]

GB.params$diet[, 'Sea Birds' := c(rep(NA, 4), 0.034, NA, 0.137, rep(NA, 7), 0.01,
                                  NA, 0.305, 0.263, 0.068, 0.001, 5.4e-4, NA, 
                                  0.041, 0.013, rep(NA, 5), 0.126, NA)]
               
GB <- rpath(GB.params, 'Georges Bank')

#Webplot plots the resultant food web
png(file = paste(out.dir, 'Georges_Bank_EMAX_Foodweb.png', sep = ''),
    height = 1700, width = 2200, res = 200)
my.groups <- c(c(30, 1, 31), c(2, 3, 4, 7, 10, 13, 5, 9, 15, 11, 12), 
               c(6, 8, 14, 16, 18), c(27, 23, 20, 22, 19, 17),
               c(26, 24, 29, 32), c(28, 25, 21))
webplot(GB, labels = T, box.order = my.groups)
dev.off()

