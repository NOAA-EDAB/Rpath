# Initial set up for Rsim dynamic simulation

Converts the outputs from Rpath into rates for use in Rsim.

## Usage

``` r
rsim.params(
  Rpath,
  mscramble = 2,
  mhandle = 1000,
  preyswitch = 1,
  scrambleselfwt = 0,
  handleselfwt = 0,
  steps_yr = 12,
  steps_m = 1
)
```

## Arguments

- Rpath:

  R object containing a static `Rpath` model.

- mscramble:

  Base value for vulnerability in functional response; default = 2
  (mixed response).

- mhandle:

  Base value for handling time in functional response; default = 1000
  (off).

- preyswitch:

  Exponent for prey switching in functional response; default = 1 (off).

- scrambleselfwt:

  Value of 1 indicates all predators overlap in the foraging arena while
  0 treats predators individually.

- handleselfwt:

  Value of 1 indicates all prey overlap in the arena and contribute to
  saturation while 0 treats prey individually.

- steps_yr:

  Number of time steps per year; default = 12.

- steps_m:

  Number of time steps per month; default = 1.

## Value

Returns an object of class `Rsim.params`, a list of 39 objects that is
passed to the `rsim.run` function via the `rsim.scenario` function.

- `NUM_GROUPS`, number of total model groups

- `NUM_LIVING`, number of living model groups

- `NUM_DEAD`, number of detritus model groups

- `NUM_GEARS`, number of fishery model groups

- `NUM_BIO`, number of living + detritus model groups

- `spname`, named `NUM_GROUPS`+1 length character vector beginning with
  "Outside" containing names of all model groups

- `spnum`, named `NUM_GROUPS`+1 length character vector beginning with
  "Outside" numbered 0, containing numbers of all model groups

- `B_BaseRef`, named `NUM_GROUPS`+1 length character vector beginning
  with "Outside" containing Rpath base biomass of all model groups

- `MzeroMort`, named `NUM_GROUPS`+1 length character vector beginning
  with "Outside" containing unaccounted mortality, calculated as PB \*
  (1.0 - EE), of all model groups

- `UnassimRespFrac`, named `NUM_GROUPS`+1 length character vector
  beginning with "Outside" containing the proportion of consumption that
  goes to detritus of all model groups

- `ActiveRespFrac`, named `NUM_GROUPS`+1 length character vector
  beginning with "Outside" containing the proportion of consumption that
  is "lost to heat" for all model groups

- `FtimeAdj`, named `NUM_GROUPS`+1 length character vector beginning
  with "Outside" containing rate of change of feeding time, currently
  set to 0 for all model groups

- `FtimeQBOpt`, named `NUM_GROUPS`+1 length character vector beginning
  with "Outside" containing base QB for all consumer model groups, or
  base PB for primary producers

- `PBopt`, named `NUM_GROUPS`+1 length character vector beginning with
  "Outside" containing base PB for all model groups

- `NoIntegrate`, named `NUM_GROUPS`+1 length character vector beginning
  with "Outside" containing flag set to 0 for high turnover model groups
  and set to `spnum` for all others

- `Handleself`, named `NUM_GROUPS`+1 length character vector beginning
  with "Outside" containing flag for handling time influence, set for
  all model groups with function argument `handleselfwt`, default 0 for
  individual prey handling time

- `Scrambleself`, named `NUM_GROUPS`+1 length character vector beginning
  with "Outside" containing flag for predator density influence, set for
  all model groups with function argument `scrambleselfwt`, default 0
  for individual predator density dependent predation

- `PreyFrom`, numeric vector length `NumPredPreyLinks`+1, spnum of prey
  for each predator prey interaction pair in the model

- `PreyTo`, numeric vector length `NumPredPreyLinks`+1, spnum of
  predator for each predator prey interaction pair in the model

- `QQ`, numeric vector length `NumPredPreyLinks`+1, base consumption
  rate for each predator prey interaction pair in the model

- `DD`, numeric vector length `NumPredPreyLinks`+1, handling time effect
  on functional response for each predator prey pair, set for all model
  groups with function argument `mhandle`, default = 1000

- `VV`, numeric vector length `NumPredPreyLinks`+1, vulnerability effect
  on functional response for each predator prey pair, set for all model
  groups with function argument `mscramble`, default = 2

- `HandleSwitch`, numeric vector length `NumPredPreyLinks`+1, prey
  density dependence effect on functional response for each predator
  prey pair, set for all model groups with function argument
  `preyswitch`, default = 1

- `PredPreyWeight`, numeric vector length `NumPredPreyLinks`+1, relative
  weight of individual predator to total predators for each predator
  prey pair, used if `scrambleselfwt`\>0

- `PreyPreyWeight`, numeric vector length `NumPredPreyLinks`+1, relative
  weight of individual prey to total prey for each predator prey pair,
  used if `handleselfwt`\>0

- `NumPredPreyLinks`, number of predator to prey linkages over all
  groups in the model

- `FishFrom`, numeric vector length `NumFishingLinks`+1, spnum of
  landing and discard for each fishery interaction in the model

- `FishThrough`, numeric vector length `NumFishingLinks`+1, spnum of
  gear type for each fishery interaction in the model

- `FishQ`, numeric vector length `NumFishingLinks`+1, landings or
  discards relative to base fished group biomass for each fishery
  interaction in the model

- `FishTo`, numeric vector length `NumFishingLinks`+1, spnum of sink for
  each fishery interaction in the model ("Outside" for landings or
  detritus group for discards)

- `NumFishingLinks`, number of model group landings and discards to
  fishery links over all groups in the model

- `DetFrac`, numeric vector length `NumDetLinks`+1, fraction of detritus
  going to DetTo for each living and detritus group in the model

- `DetFrom`, numeric vector length `NumDetLinks`+1, spnum flowing to
  detritus for each living and detritus group in the model

- `DetTo`, numeric vector length `NumDetLinks`+1, spnum of detritus sink
  for each living and detritus group in the model

- `NumDetLinks`, number of model group links to detritus over all groups
  in the model

- `BURN_YEARS`, number of model run burn-in (spin up) years, default
  value -1

- `COUPLED`, number to control species interactions, value of 0 allows
  density dependent non-interacting species, default value 1

- `RK4_STEPS`, number of RK4 integration steps per month, default value
  4

- `SENSE_LIMIT`, numeric vector of multipliers on biomass determining
  acceptable range for continuing a model run within `BURN_YEARS`

## See also

Other Rsim functions:
[`extract.node()`](https://noaa-edab.github.io/Rpath/reference/extract.node.md),
[`rsim.plot()`](https://noaa-edab.github.io/Rpath/reference/rsim.plot.md),
[`rsim.run()`](https://noaa-edab.github.io/Rpath/reference/rsim.run.md),
[`rsim.scenario()`](https://noaa-edab.github.io/Rpath/reference/rsim.scenario.md),
[`rsim.step()`](https://noaa-edab.github.io/Rpath/reference/rsim.step.md),
[`write.Rsim()`](https://noaa-edab.github.io/Rpath/reference/write.Rsim.md)

## Examples

``` r
# Read in Rpath parameter file and generate model object
Rpath <- rpath(AB.params)
# Create default dynamic parameters from Rpath model
Rsim.params <- rsim.params(Rpath) 
```
