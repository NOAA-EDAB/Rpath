# ---- Modify this toggle to TRUE to generate the baseline files. ----
# ---- Reset it back to FALSE to run the tests. ----------------------
# CREATE_BASELINE_FILES <- TRUE
CREATE_BASELINE_FILES <- FALSE

NUMBER_OF_STEPS <- 5 # should be an odd multiple of nrows=600 (i.e., 5,15,30)
FACTOR_VALUE <- 5
SEED_VALUE   <- 1 
SEED_OFFSET <- 2000
TOLERANCE_VALUE <- 1e-5
RUN_QUIET <- TRUE
YLIMIT_DIFFERENCE_PLOTS <- 0.05
PLOT_TYPE <- 1 # 1 = Baseline and Current superimposed, 2 = difference of (Current-Baseline)
PLOT_SHOW <- 1 # 1 - All Plots, 2 = Only plots reflecting test errors # Not sure if can be implemented
# INPUT_DATA_DIR_BASELINE  <- 'data/input/baseline'
INPUT_DATA_DIR_CURRENT   <- here::here('tests/testthat/data/input')
INPUT_DATA_DIR_BASELINE  <- INPUT_DATA_DIR_CURRENT
OUTPUT_DATA_DIR          <- here::here('tests/testthat/data/output')
JITTER_AMOUNT_PCT <- 1 # This is 1% jitter
