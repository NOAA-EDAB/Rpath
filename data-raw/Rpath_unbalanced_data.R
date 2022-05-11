## Code to create unbalanced ecopath objects for Rpath examples
## Doesn't run on package install, needs to be run in repo and then Rdata
## files created in data and inst/extdata need to be pushed/committed.
##
## Only re-run if changing input data changes, or if structure of unbalanced
## ecosystem object changes.
##
## Also create entry in R/data.R file for each ecosystem.

## This must be run from the Rpath package's root directory (the
## directory with the /data, /R, /data-raw and /inst folders).

library(usethis)
library(Rpath)

Ecosense.EBS <- read.rpath.params("data-raw/Ecosense_EBS_base.csv",
                                  "data-raw/Ecosense_EBS_diet.csv",
                                  "data-raw/Ecosense_EBS_pedigree.csv"
)
Ecosense.GOA <- read.rpath.params("data-raw/Ecosense_GOA_base.csv",
                                  "data-raw/Ecosense_GOA_diet.csv",
                                  "data-raw/Ecosense_GOA_pedigree.csv"
)
Ecosense.ECS <- read.rpath.params("data-raw/Ecosense_ECS_base.csv",
                                  "data-raw/Ecosense_ECS_diet.csv",
                                  "data-raw/Ecosense_ECS_pedigree.csv"
)
usethis::use_data(Ecosense.EBS, overwrite = TRUE)
usethis::use_data(Ecosense.GOA, overwrite = TRUE)
usethis::use_data(Ecosense.ECS, overwrite = TRUE)
