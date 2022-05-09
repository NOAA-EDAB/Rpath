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

EBS.sensemodel <- read.rpath.params("")

usethis::use_data(DATASET, overwrite = TRUE)
