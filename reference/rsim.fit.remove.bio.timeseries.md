# Remove Biomass Timeseries from Fitting Data

Removes specific biomass observation timeseries from an \`Rsim\`
scenario object based on their unique identification strings.

## Usage

``` r
rsim.fit.remove.bio.timeseries(Rsim.scenario, timeseries)
```

## Arguments

- Rsim.scenario:

  An \`Rsim\` scenario object containing the fitting data setup.

- timeseries:

  A character vector of \`ID\` strings representing the timeseries to be
  removed (e.g., \`"SpeciesName:SurveySource"\`). Existing timeseries
  IDs can be listed using the \`rsim.fit.list.bio.series\` function.

## Value

Modifies and returns the \`Rsim.scenario\` object, with the specified
timeseries removed from the \`Rsim.scenario\$fitting\$Biomass\` data
frame.

## Details

This is a utility function used to drop specific observation datasets
(e.g., a specific survey for a specific species) prior to running a
model fit. It filters the \`scene\$fitting\$Biomass\` data frame,
keeping only the rows where the \`ID\` does NOT match the values
provided in the \`timeseries\` argument.
