# Package index

## Rpath

Functions for creating a mass balance food web

- [`create.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/create.rpath.params.md)
  : Creates a shell for the Rpath.params object
- [`check.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/check.rpath.params.md)
  : Check Rpath parameter files
- [`read.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/read.rpath.params.md)
  : Read Rpath parameters from .csv files
- [`rpath()`](https://noaa-edab.github.io/Rpath/reference/rpath.md) :
  Ecopath module of Rpath
- [`rpath.stanzas()`](https://noaa-edab.github.io/Rpath/reference/rpath.stanzas.md)
  : Calculate biomass and consumption for multistanza groups
- [`rpath.consumers()`](https://noaa-edab.github.io/Rpath/reference/rpath.consumers.md)
  : List of Rpath consumer groups
- [`rpath.detrital()`](https://noaa-edab.github.io/Rpath/reference/rpath.detrital.md)
  : List of Rpath detrital groups
- [`rpath.gears()`](https://noaa-edab.github.io/Rpath/reference/rpath.gears.md)
  : List of Rpath gears groups
- [`rpath.groups()`](https://noaa-edab.github.io/Rpath/reference/rpath.groups.md)
  : Rpath functional group names
- [`rpath.living()`](https://noaa-edab.github.io/Rpath/reference/rpath.living.md)
  : List of Rpath living groups
- [`rpath.mixotrophs()`](https://noaa-edab.github.io/Rpath/reference/rpath.mixotrophs.md)
  : List of Rpath mixotroph groups
- [`rpath.producers()`](https://noaa-edab.github.io/Rpath/reference/rpath.producers.md)
  : List of Rpath producer groups

## Rsim

Functions for creating dynamic simulations of the mass balance model

- [`rsim.scenario()`](https://noaa-edab.github.io/Rpath/reference/rsim.scenario.md)
  : Specify a dynamic simulation scenario
- [`rsim.run()`](https://noaa-edab.github.io/Rpath/reference/rsim.run.md)
  : Run a dynamic simulation
- [`rsim.step()`](https://noaa-edab.github.io/Rpath/reference/rsim.step.md)
  : Function to add a time step to an Rsim run
- [`rsim.params()`](https://noaa-edab.github.io/Rpath/reference/rsim.params.md)
  : Initial set up for Rsim dynamic simulation
- [`rsim.fishing()`](https://noaa-edab.github.io/Rpath/reference/rsim.fishing.md)
  : Initialize fishing matrix for Rsim scenario
- [`rsim.forcing()`](https://noaa-edab.github.io/Rpath/reference/rsim.forcing.md)
  : Initialize forcing matrix for Rsim scenario
- [`rsim.deriv()`](https://noaa-edab.github.io/Rpath/reference/rsim.deriv.md)
  : Calculate the derivatives for a time step
- [`rsim.diet()`](https://noaa-edab.github.io/Rpath/reference/rsim.diet.md)
  : Output consumption of each prey by an individual predator group from
  Rsim run
- [`rsim.mort()`](https://noaa-edab.github.io/Rpath/reference/rsim.mort.md)
  : Output biomass removed by each predator for an individual prey group
  from Rsim run
- [`rsim.stanzas()`](https://noaa-edab.github.io/Rpath/reference/rsim.stanzas.md)
  : Initialize list of age structured group (stanza) parameters for
  dynamic simulation
- [`rsim.state()`](https://noaa-edab.github.io/Rpath/reference/rsim.state.md)
  : Initialize starting state matrix for Rsim scenario
- [`rsim.sense()`](https://noaa-edab.github.io/Rpath/reference/rsim.sense.md)
  : Ecosense function for Rpath
- [`get.rsim.fishing()`](https://noaa-edab.github.io/Rpath/reference/get.rsim.fishing.md)
  : Retrieve fishing forcing parameters from an Rsim scenario
- [`get.rsim.forcing()`](https://noaa-edab.github.io/Rpath/reference/get.rsim.forcing.md)
  : Retrieve forcing parameters from an Rsim scenario
- [`get.rsim.params()`](https://noaa-edab.github.io/Rpath/reference/get.rsim.params.md)
  : Retrieve parameters from an Rsim scenario
- [`get.rsim.stanzas()`](https://noaa-edab.github.io/Rpath/reference/get.rsim.stanzas.md)
  : Retrieve stanza parameters from an Rsim scenario
- [`get.rsim.start_state()`](https://noaa-edab.github.io/Rpath/reference/get.rsim.start_state.md)
  : Retrieve starting state values from an Rsim scenario
- [`set.rsim.scene()`](https://noaa-edab.github.io/Rpath/reference/set.rsim.scene.md)
  : Set Rsim.scenario parameters

## Models

Example ecosysem models

- [`Ecosense.EBS`](https://noaa-edab.github.io/Rpath/reference/Ecosense.EBS.md)
  : Eastern Bering Sea 1990s Ecopath model
- [`Ecosense.ECS`](https://noaa-edab.github.io/Rpath/reference/Ecosense.ECS.md)
  : Eastern Chukchi Sea Ecopath model
- [`Ecosense.GOA`](https://noaa-edab.github.io/Rpath/reference/Ecosense.GOA.md)
  : Gulf of Alaska (west/central) 1990s Ecopath model
- [`AB.params`](https://noaa-edab.github.io/Rpath/reference/AB.params.md)
  : AB.params
- [`REco.params`](https://noaa-edab.github.io/Rpath/reference/REco.params.md)
  : REco.params

## Adjustment Tools

Functions for adjusting Rsim scenarios

- [`adjust.scenario()`](https://noaa-edab.github.io/Rpath/reference/adjust.scenario.md)
  : Adjust Rsim.scenario parameters
- [`adjust.fishing()`](https://noaa-edab.github.io/Rpath/reference/adjust.fishing.md)
  : Adjust Fishing Mortality
- [`adjust.forcing()`](https://noaa-edab.github.io/Rpath/reference/adjust.forcing.md)
  : Adjust Forcing Parameters

## Visualization Tools

Functions for visualizing rpath and rsim outputs

- [`webplot()`](https://noaa-edab.github.io/Rpath/reference/webplot.md)
  : Plot routine for Ecopath food web
- [`ggwebplot()`](https://noaa-edab.github.io/Rpath/reference/ggwebplot.md)
  : ggplot routine for Rpath food web
- [`rsim.plot()`](https://noaa-edab.github.io/Rpath/reference/rsim.plot.md)
  : Plot routine for Rsim runs
- [`stanzaplot()`](https://noaa-edab.github.io/Rpath/reference/stanzaplot.md)
  : Plot routine for Ecopath multistanzas
- [`summarize.for.webplot()`](https://noaa-edab.github.io/Rpath/reference/summarize.for.webplot.md)
  : Prepare food web model to plot

## Output Tools

Functions for writing and exporting outputs from Rpath

- [`write.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/write.rpath.params.md)
  : Write Rpath parameters to .csv files
- [`write.Rpath()`](https://noaa-edab.github.io/Rpath/reference/write.Rpath.md)
  : Write function for Rpath food web object
- [`write.Rsim()`](https://noaa-edab.github.io/Rpath/reference/write.Rsim.md)
  : Write function for Rsim run
- [`frate.table()`](https://noaa-edab.github.io/Rpath/reference/frate.table.md)
  : Fishing Mortality Table
- [`extract.node()`](https://noaa-edab.github.io/Rpath/reference/extract.node.md)
  : Extract node/group data from an Rsim run

## Visualization Tools

Functions for visualizing rpath and rsim outputs

## Network Tools

Functions for trophic network analysis

- [`MTI()`](https://noaa-edab.github.io/Rpath/reference/MTI.md) :
  Quantify mixed trophic impacts

## EwE Compatibility Tools

Functions for importing EwE models

- [`import.eiixml()`](https://noaa-edab.github.io/Rpath/reference/import.eiixml.md)
  : Reads in EwE exported XML file and parses into data frames
- [`create.rpath.from.eiixml()`](https://noaa-edab.github.io/Rpath/reference/create.rpath.from.eiixml.md)
  : Creates an Rpath object from an EwE exported model (EIIXML format)

## Model Fitting Tools

Functions for fitting Ecosim runs to data

- [`fitcatch.to.forcecatch()`](https://noaa-edab.github.io/Rpath/reference/fitcatch.to.forcecatch.md)
  : Convert Fitting Catch to Forced Catch in Rsim Scenario
- [`get.pp.link()`](https://noaa-edab.github.io/Rpath/reference/get.pp.link.md)
  : Get Predator-Prey Link Index
- [`read.fitting.biomass()`](https://noaa-edab.github.io/Rpath/reference/read.fitting.biomass.md)
  : Read and Process Biomass Fitting Data
- [`read.fitting.catch()`](https://noaa-edab.github.io/Rpath/reference/read.fitting.catch.md)
  : Read and Process Catch Fitting Data
- [`rsim.fit.apply()`](https://noaa-edab.github.io/Rpath/reference/rsim.fit.apply.md)
  : Apply Fitted Parameters to Rpath Simulation Scene Helper function
  that maps a vector of optimized parameter modifications back into the
  \`Rpath\` simulation parameters object (\`scene.params\`). Since this
  returns a scenario params object instead of a full scenario, it is
  primarily intended for internal use.
- [`rsim.fit.list.bio.bysurvey()`](https://noaa-edab.github.io/Rpath/reference/rsim.fit.list.bio.bysurvey.md)
  : List Biomass Timeseries Groups by Survey Source
- [`rsim.fit.list.bio.series()`](https://noaa-edab.github.io/Rpath/reference/rsim.fit.list.bio.series.md)
  : List Biomass Timeseries Groups and Sources
- [`rsim.fit.list.catch.series()`](https://noaa-edab.github.io/Rpath/reference/rsim.fit.list.catch.series.md)
  : List Catch Timeseries Groups
- [`rsim.fit.obj()`](https://noaa-edab.github.io/Rpath/reference/rsim.fit.obj.md)
  : Calculate Objective Function for Rsim Model Fit
- [`rsim.fit.obj.species()`](https://noaa-edab.github.io/Rpath/reference/rsim.fit.obj.species.md)
  : Extract Objective Function Fits by Species
- [`rsim.fit.remove.bio.timeseries()`](https://noaa-edab.github.io/Rpath/reference/rsim.fit.remove.bio.timeseries.md)
  : Remove Biomass Timeseries from Fitting Data
- [`rsim.fit.run()`](https://noaa-edab.github.io/Rpath/reference/rsim.fit.run.md)
  : Run Fitted Rpath Simulation and Evaluate Objective Function
- [`rsim.fit.set.bio.wt()`](https://noaa-edab.github.io/Rpath/reference/rsim.fit.set.bio.wt.md)
  : Set Biomass Fitting Data Weights
- [`rsim.fit.set.catch.wt()`](https://noaa-edab.github.io/Rpath/reference/rsim.fit.set.catch.wt.md)
  : Set Catch Fitting Data Weights
- [`rsim.fit.set.q()`](https://noaa-edab.github.io/Rpath/reference/rsim.fit.set.q.md)
  : Set Catchability (q) for Biomass Fitting Data
- [`rsim.fit.table()`](https://noaa-edab.github.io/Rpath/reference/rsim.fit.table.md)
  : Extract Summary Table of Model Fits
- [`rsim.fit.update()`](https://noaa-edab.github.io/Rpath/reference/rsim.fit.update.md)
  : Update Rpath Scene with Fitted Parameters
- [`rsim.plot.fit()`](https://noaa-edab.github.io/Rpath/reference/rsim.plot.fit.md)
  : Plot Biomass and Catch Fitting Results
- [`rsim.plot.fitspecies()`](https://noaa-edab.github.io/Rpath/reference/rsim.plot.fitspecies.md)
  : Plot Full Biomass and Catch Fits for a Single Species Generates a
  multi-panel figure displaying the simulation fits for a single species
  across all available biomass observation sources, alongside its catch
  fitting data.
- [`rsim.predprey.table()`](https://noaa-edab.github.io/Rpath/reference/rsim.predprey.table.md)
  : Extract Predator-Prey Interaction Table
