# Rpath 1.1.0

- Added ability to import .eiixml files via the `create.rpath.from.eiixml()` function
- Added balance functionality to estimate P/B from input Biomass and EE, also
allowing a missing Q/B to be estimated from input PC and the estimated P/B
- Corrected accounting for single-stage interdetrital flows during balance
- Behavior change: if all of P/B, Q/B, and PC are supplied as inputs for a group, recalculate PC during balance to ensure internal model consistency (instead of leaving PC unchanged and possibly inaccurate)
- Improved error messages in `rpath()` balance routine to help diagnose models that are missing parameters
- Improved warning message content in `check.rpath.params()` to aid diagnosis
- Modified `rpath.stanzas()` so it no longer produces errors when called with a model that has no multistanza groups (instead it returns the model unchanged)
- Added Western Bering Sea model .eiixml and EwE output csv files as an import example
- Added `Convert_EwE_to_Rpath` vignette to describe .eiixml import process
- Tested eiixml input routines on over 150 models available via EcoBase with most producing consistent results between Rpath and EwE; some remaining balancing differences between Rpath and EwE are noted in the `Convert_EwE_to_Rpath` vignette
- Fixed some table ordering issues with stanza inputs
- Documentation fixes


# Rpath 1.0.0

- Fix group parameter reference in `adjust.scenario`
- Fix Force bymort double-counting bug `ecosim.cpp`
- rsim step dyt fix in `ecosim_multistep`
- Removed old datasets
- Switch from single to double operators in `ecosim.cpp`
- Added examples to all functions
- Updated documentation/website
- Added package hex
- Added a bibliography file
- Created issue templates
- Added contributors guideline and code of conduct

# Rpath 0.9.1

- New function documentation added

# Rpath 0.9.0

- Added Ecosense routines and vignettes
- Added example ecosystems
- Added Unit Tests
- Fixed bug in stanza month calculations that affected some species

# Rpath 0.8.0

- Added bioenergetic code
- adjust functions can now supply values by individual months rather than by years

# Rpath 0.7.0

- Added CI workflows
- Added/updated documentation (pkgdown)
- Stanza changes

# Rpath 0.6.0

- Added BAB to multistanzas
- adjust.scenario can be applied to multiple groups at once rather than one at a time.
- Changed variable names to align with manuscript
  - BB to Biomass
  - CC to Catch
  - GS to Unassim
  - Catch to Landings
  - Recruits from r to R
  - WWa to QageS
  - surv to survive_L
  - B to Biomass in Pedigree param list
  - All "byx" functions to "ForcedX"
  - EFFORT to ForcedEffort
  - CATCH and FRATE to ForcedCatch and ForcedFRate
- Change default values of scramble flags from 1 to 0. Zero indicating the flag is off, 1 indicating it is on.

# Rpath 0.5.1

- Fixed mixotroph diet calculations and MTI calculation.

# Rpath 0.5.0

- Initial Release
- Acknowledgement of package submission to CRAN along with manuscript (see release notes)
- Fixed bug with adjust.fishing so one value can be applied over a range of sim.years
- Fixed bug with original scenario file and stanza files being overwritten by C code.
- Fixed bug with egg stanza calculation where the Wmat matrix was indexed wrong inthe calculation.
- Updated adjustment functions.
- Renamed parameter "year" to "sim.year".
- Removed the "gear" parameter from adjust fishing as it is a group and there is never a group/fleet combo that needs to be specified.
- Fixed indexing issue with rsim.step and issue with rsim.run modifying initial rsim.scenario data table.
- ecosim scenario objects now parameterized with a vector of years for labeling, e.g. years=c(1970,2017), not a single number of years.
- rsim.run will run for the years specified (e.g. years=c(2013,2014) will save outputs in those two year slots). Haven't carried this change over to the multistep function.
- Labeling in the fishing/forcing matrices (columns labelled by species, rows labeled by years (for annual matrices, no labeling for monthly rows)
- The addition of force_bybio (biomass forced to stay at input value)
- Added gear specific catch tracking and new extract.node function.
- Fixed issue with consumption calculations due to the addition of diet import in the diet matrix.

