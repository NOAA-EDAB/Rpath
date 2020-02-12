Rpath
=====
This is the repository for Rpath - the R version of the Ecopath with Ecosim model

The new public facing repository (NOAA-EDAB/Rpath) was established for release with 
the Rpath manuscript (Lucey et al. submitted).  Previous repository is no longer maintained.

The initial files are the code sent by Kerim to Sean on March 31, 2014 (ecopath_sim_r_v0_04_toSean_31_03_14.zip)
Files have been renamed to drop the version number and take advantage of this VCS

Beta version notes:

v0.0.1.1 - Fixed mixotroph diet calculations and MTI calculation.

v0.0.1.0 - Acknowledgement of package submission along with manuscript

v0.0.0.917 - Fixed bug with adjust.fishing so one value can be applied over a range of sim.years

v0.0.0.916 - Fixed bug with original scenario file and stanza files being overwritten by C code.

v0.0.0.915 - Fixed bug with egg stanza calculation where the Wmat matrix was indexed wrong inthe calculation.

v0.0.0.9014 - Updated adjustment functions.  Renamed parameter "year" to "sim.year".  Removed the "gear" parameter from adjust fishing as it is a group and there is never a group/fleet combo that needs to be specified.

v0.0.0.9013 - Fixed indexing issue with rsim.step and issue with rsim.run modifying initial rsim.scenario data table.

v0.0.0.911 - ecosim scenario objects now parameterized with a vector of years for labeling, e.g. years=c(1970,2017), not a single number of years.
- rsim.run will run for the years specified (e.g. years=c(2013,2014) will save outputs in those two year slots).  Haven't carried this change over to the multistep function.
- labeling in the fishing/forcing matrices (columns labelled by species, rows labeled by years (for annual matrices, no labeling for monthly rows)
- the addition of force_bybio (biomass forced to stay at input value)

v0.0.0.906 - Added gear specific catch tracking and new extract.node function.

v0.0.0.905 - Caught up with Kerim's edits.

v0.0.0.904 - Data.table package updated and created a few bugs that have been solved.

v0.0.0.903 - Fixed issue with consumption calculations due to the addition of
diet import in the diet matrix.

Did not keep good track of notes prior to version 0.0.0.9003 - My Bad


*This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.*
