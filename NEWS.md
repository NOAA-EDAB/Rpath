# Rpath 0.0.1.2

## Major Changes
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

# Rpath 0.0.1.1

- Fixed mixotroph diet calculations and MTI calculation.

# Rpath 0.0.1.0

- Acknowledgement of package submission along with manuscript

# Rpath 0.0.0.917

- Fixed bug with adjust.fishing so one value can be applied over a range of sim.years

# Rpath 0.0.0.916

- Fixed bug with original scenario file and stanza files being overwritten by C code.

# Rpath 0.0.0.915

- Fixed bug with egg stanza calculation where the Wmat matrix was indexed wrong inthe calculation.

# Rpath 0.0.0.9014

- Updated adjustment functions.
- Renamed parameter "year" to "sim.year".
- Removed the "gear" parameter from adjust fishing as it is a group and there is never a group/fleet combo that needs to be specified.

# Rpath 0.0.0.9013

- Fixed indexing issue with rsim.step and issue with rsim.run modifying initial rsim.scenario data table.

# Rpath 0.0.0.911

- ecosim scenario objects now parameterized with a vector of years for labeling, e.g. years=c(1970,2017), not a single number of years.
- rsim.run will run for the years specified (e.g. years=c(2013,2014) will save outputs in those two year slots). Haven't carried this change over to the multistep function.
- labeling in the fishing/forcing matrices (columns labelled by species, rows labeled by years (for annual matrices, no labeling for monthly rows)
- the addition of force_bybio (biomass forced to stay at input value)

# Rpath 0.0.0.906

- Added gear specific catch tracking and new extract.node function.

# Rpath 0.0.0.905

- Caught up with Kerim's edits.

# Rpath 0.0.0.904

- Data.table package updated and created a few bugs that have been solved.

# Rpath 0.0.0.903

- Fixed issue with consumption calculations due to the addition of diet import in the diet matrix.

Did not keep good track of notes prior to version 0.0.0.9003 - My Bad

The initial files were the code sent by Kerim to Sean on March 31, 2014 (ecopath_sim_r_v0_04_toSean_31_03_14.zip) Files have been renamed to drop the version number and take advantage of this VCS.