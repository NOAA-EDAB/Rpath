# Convert EwE to Rpath

``` r

library(Rpath)
```

If you have developed your model in the Ecopath with Ecosim (EwE)
software package and you now want to bring it into Rpath for analysis,
there are tools available in `Rpath` to do just that. First, you need to
export your model from EwE to XML format, using EwE’s
`File -> Export Model -> To XML` menu option to create a .eiixml file
(all EwE functionality described here was tested using EwE version
6.7.0).

For this example, the Western Bering Sea EwE model available for
download on EcoBase (Aydin et al. 2002) was first imported into EwE,
using EwE’s Import Model menu option, then exported to the file
`Western_Bering_Sea.eiixml`. We can read in this file and convert it to
an `Rpath` object using the function
[`create.rpath.from.eiixml()`](https://noaa-edab.github.io/Rpath/reference/create.rpath.from.eiixml.md).

``` r

# Read in the model
eiixml_file <- system.file("extdata/xml", "Western_Bering_Sea.eiixml", package = "Rpath")
model <- create.rpath.from.eiixml(eiixml_file)
```

Currently, the
[`create.rpath.from.eiixml()`](https://noaa-edab.github.io/Rpath/reference/create.rpath.from.eiixml.md)
will only import input parameters related to the unbalanced Ecopath
model, and will not import Ecosim information or other supplementary
tables. However, the full set of tables in an .eiixml file can be
examined (or manually imported) using the
[`import.eiixml()`](https://noaa-edab.github.io/Rpath/reference/import.eiixml.md)
helper function.

The variable `model` is a list object (an unbalanced Ecopath model)
containing the data objects required to balance your model and run
simulations. There are 4 objects that comprise an `Rpath` object:

- model - including basic inputs, detritus fate, other production, and
  fishery landings/discards.
- diet - a diet matrix, including diet import.
- stanzas - parameters for multi-stanza functional groups.
- pedigree - data pedigree (quality) for the model.

Prepare your imported model for balancing by first calculating
multistanza parameters that are based on age structure (such as B and
Q/B), using the
[`rpath.stanzas()`](https://noaa-edab.github.io/Rpath/reference/rpath.stanzas.md)
function, then using the
[`check.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/check.rpath.params.md)
function to ensure all of the model’s parameters are entered (If your
model has no multistanza groups,
[`rpath.stanzas()`](https://noaa-edab.github.io/Rpath/reference/rpath.stanzas.md)
can be called but will not change the model). The
[`check.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/check.rpath.params.md)
function will make sure the data was read in correctly, and produce
warnings in situations that might prevent the model from balancing; for
example, if a functional group is missing too many parameters to be
balanced, or if diet compositions for a predator do not sum to 1.

``` r

model <- rpath.stanzas(model)
check.rpath.params(model)
#> Rpath parameter file is functional.
```

## Model parameters

To compare Rpath balance results with EwE, first balance the model in
EwE and use EwE’s “Save to a csv” icon on EwE’s Outputs -\> Basic
estimates tab (upper right corner). It is recommended that you first
increase the precision (digits displayed) in EwE to 8 digits. The
following script, using the dplyr library, can then compare the
resulting balances.

``` r


# dplyr is used to match the Rpath and EwE balance results
  library(dplyr, warn.conflicts=F)

# EwE output files - eiixml file and Basic Estimates csv
  eiifile <- system.file("extdata/xml", "Western_Bering_Sea.eiixml", package = "Rpath")
  
  zip_path <- system.file("extdata", "xml", "western_bering.zip", package = "Rpath")
  csvfile <- read.csv(unz(zip_path, "Western Bering Sea-Basic estimates.csv"))

# Load and balance the model in Rpath  
  unbal  <- create.rpath.from.eiixml(eiifile)
  unbal <- rpath.stanzas(unbal)
  check.rpath.params(unbal)
#> Rpath parameter file is functional.
  bal <- rpath(unbal)
  
# Load the csv file and clean up the format
  csv.out <- csvfile
  # If there's no stanzas, the total mortality column will be missing from the CSV - add an extra column
  if(length(csv.out)<13){csv.out <- data.frame(append(csv.out, list(tm=NA), after=6))}
  # cleaned column names matching EwE output
  names(csv.out)<-  c("X", "Group.name", "Trophic.level", "Hab.area", "Biomass.in.habitat.area",
                    "Biomass", "Total.mortality", "Production.biomass", "Consumption.biomass",   
                    "Ecotrophic.Efficiency", "Production.consumption", "Biomass.accumulation",
                    "BA.rate")  
  # Drop placeholder lines that head stanza groups
  csv.out <- csv.out[!is.na(csv.out$X),]
  # Clean names to rpath standard
  csv.out$group <- janitor::make_clean_names(csv.out$Group.name)
  
# Create a data frame of rpath (balanced model) output 
  rpath.dat <- data.frame(group=bal$Group,type=bal$type, tl=bal$TL, biomass=bal$Biomass,
                          pb=bal$PB,qb=bal$QB, EE=bal$EE, pc=bal$GE, 
                          ba=bal$BA)[c(rpath.living(bal), rpath.detrital(bal)),]
  
# Join rpath and EwE balance outputs, using dplyr
  rpath_ewe_table <- rpath.dat %>% 
    dplyr::left_join(csv.out, by="group") %>%
    dplyr::mutate( biomass_test = abs(biomass-Biomass)/Biomass, 
                   ee_test      = abs(EE-Ecotrophic.Efficiency)) %>%
    dplyr::select(group, rpath.biomass=biomass, ewe.biomass=Biomass, prop.bio.diff=biomass_test,
                  rpath.EE=EE, ewe.EE=Ecotrophic.Efficiency, EE.diff=ee_test)

  rpath_ewe_table 
#>                   group rpath.biomass ewe.biomass prop.bio.diff     rpath.EE
#> 1         baleen_whales      0.391000    0.391000  0.000000e+00 0.9891304348
#> 2        toothed_whales      0.042000    0.042000  2.380952e-08 0.0000000000
#> 3          sperm_whales      0.020000    0.020000  0.000000e+00 0.0000000000
#> 4  walrus_bearded_seals      0.262000    0.262000  3.816794e-08 0.4284351145
#> 5                 seals      0.097000    0.097000  3.092783e-08 0.5962199313
#> 6     steller_sea_lions      0.035000    0.035000  0.000000e+00 0.3500000000
#> 7              seabirds      0.010000    0.010000  2.000000e-08 0.0000000000
#> 8         adult_pollock     15.000000   15.000000  0.000000e+00 0.9497486121
#> 9      juvenile_pollock      3.757000    3.757000  7.985095e-09 0.4800152741
#> 10          pacific_cod      3.187000    3.187000  1.255099e-08 0.3541614076
#> 11            p_halibut      0.083000    0.083000  3.614458e-08 0.9878137841
#> 12     greenland_turbot      0.058000    0.058000  3.448276e-08 0.9465517241
#> 13  arrowtooth_flounder      0.052000    0.052000  1.923077e-08 0.8895514435
#> 14       small_flatfish      0.992000    0.992000  2.016129e-08 0.7397461896
#> 15               skates      0.271000    0.271000  0.000000e+00 0.3782287823
#> 16    sculpins_rockfish      0.677000    0.677000  1.477105e-08 0.9444742460
#> 17          macrouridae      1.156000    1.156000  1.730104e-08 0.6748010381
#> 18            zoarcidae      0.900000    0.900000  2.222222e-08 0.9880651111
#> 19          tanner_crab      0.083000    0.083000  3.614458e-08 0.3979768578
#> 20            snow_crab      0.249000    0.249000  0.000000e+00 0.8712332284
#> 21            king_crab      0.119000    0.119000  0.000000e+00 0.8268279357
#> 22               shrimp      2.105753    2.105753  4.777329e-09 0.9000000000
#> 23             epifauna    114.962000  114.961998  1.739705e-08 0.2450068922
#> 24              infauna    125.687000  125.686996  3.182509e-08 0.2740845025
#> 25    benthic_amphipods     13.812000   13.812000  2.172024e-08 0.9944898735
#> 26      pacific_herring      0.787000    0.787000  0.000000e+00 0.8228223326
#> 27       pacific_salmon      0.039000    0.039000  2.564102e-08 0.3907632487
#> 28          cephalopods      4.830000    4.830000  1.656315e-08 0.9109046695
#> 29          forage_fish     19.419371   19.419371  3.826047e-08 0.9000000000
#> 30            jellyfish      1.400000    1.400000  1.428571e-08 0.0164664762
#> 31    large_zooplankton    120.740000  120.739998  1.656452e-08 0.5520143674
#> 32             copepods    122.620000  122.620003  2.446583e-08 0.9742557243
#> 33        phytoplankton     15.000000   15.000000  0.000000e+00 0.7382177445
#> 34     pelagic_detritus      0.000000          NA            NA          Inf
#> 35     benthic_detritus      0.000000          NA            NA          Inf
#> 36        detrital_pool   6930.982701          NA            NA 0.0001562706
#>        ewe.EE      EE.diff
#> 1  0.98913050 6.521739e-08
#> 2  0.00000000 0.000000e+00
#> 3  0.00000000 0.000000e+00
#> 4  0.42843521 9.549618e-08
#> 5  0.59622002 8.872852e-08
#> 6  0.35000011 1.100000e-07
#> 7  0.00000000 0.000000e+00
#> 8  0.94974858 3.209067e-08
#> 9  0.48001531 3.589832e-08
#> 10 0.35416141 2.425322e-09
#> 11 0.98781377 1.409639e-08
#> 12 0.94655168 4.413793e-08
#> 13 0.88955140 4.350962e-08
#> 14 0.73974621 2.036136e-08
#> 15 0.37822881 2.771218e-08
#> 16 0.94447428 3.403988e-08
#> 17 0.67480099 4.806228e-08
#> 18 0.98806512 8.888889e-09
#> 19 0.39797691 5.216867e-08
#> 20 0.87123322 8.413654e-09
#> 21 0.82682788 5.571429e-08
#> 22 0.90000021 2.100000e-07
#> 23 0.24500690 7.818708e-09
#> 24 0.27408451 7.461256e-09
#> 25 0.99449003 1.565455e-07
#> 26 0.82282239 5.739020e-08
#> 27 0.39076331 6.128171e-08
#> 28 0.91090471 4.047966e-08
#> 29 0.89999992 8.000000e-08
#> 30 0.01646648 3.809524e-09
#> 31 0.55201441 4.260628e-08
#> 32 0.97425568 4.431013e-08
#> 33 0.73821777 2.547578e-08
#> 34         NA           NA
#> 35         NA           NA
#> 36 0.00015627 6.081843e-10
```

## Differences between Rpath and EwE balances

Rpath was tested using the 150+ models available on EcoBase. In general,
differences between Rpath and EwE were within 1.0e-7 of each other for
any given estimated EE, B, or P/B, with differences due to numerical
precision of inputs. Stanza groups may show differences of up to 1.0e-3
due to accumulated differences summing across age groups. However, there
are some specific differences in functionality in some models:

- Q/B - Rpath does not currently estimate Q/B as part of the main
  balancing (i.e. to estimate Q/B, both P/B and PC must be supplied).
- Interdetrital flows - For models with three or more stages of
  interdetrital flows (that is, if detritus fate specifies thart
  Detritus A -\> Detritus B -\> Detritus C), results may differ between
  Rpath and EwE. Two-stage detrital flows (Detritus A -\> Detritus B)
  give matching results.
- Rpath will not currently estimate BA. In cases where all of B, EE, and
  P/B are supplied, Rpath will not estimate BA or make changes to the
  model, and may produce inconsistent results if one of those paramters
  isn’t set to NA.
- Rpath does not currently support BA rate, only total BA. If your model
  has a BA rate, it is recommended that you calculate the total BA in
  EwE (as an output), and enter that value into Rpath.
