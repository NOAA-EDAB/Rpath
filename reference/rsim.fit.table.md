# Extract Summary Table of Model Fits

Runs the objective function and aggregates the negative log-likelihood
fit scores by functional group for both Biomass and Catch data,
returning them in a standardized data frame.

## Usage

``` r
rsim.fit.table(Rsim.scenario, Rsim.output)
```

## Arguments

- Rsim.scenario:

  An Rsim scenario object containing the fitting data.

- Rsim.output:

  An Rsim simulation output to be evaluated.

## Value

A \`data.frame\` where the row names correspond to the model's full
species list. It contains two columns (\`Biomass\` and \`Catch\`)
representing the summed negative log-likelihood fit scores for each
group.

## Details

This function acts as a convenient wrapper around \`rsim.fit.obj()\`. It
runs the objective function in verbose mode to extract individual
observation fit scores, sums those scores by species/group, and then
aligns the results against the complete species list defined in the
model.
