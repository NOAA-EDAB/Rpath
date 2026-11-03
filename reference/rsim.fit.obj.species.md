# Extract Objective Function Fits by Species

Runs the objective function and extracts the detailed, row-by-row fit
scores strictly for specific species or functional groups from both the
Biomass and Catch data.

## Usage

``` r
rsim.fit.obj.species(Rsim.scenario, Rsim.output, species = NULL)
```

## Arguments

- Rsim.scenario:

  An \`Rsim\` scenario object containing the fitting data.

- Rsim.output:

  An Rsim simulation output to be evaluated.

- species:

  A character vector of species or group names to extract. Defaults to
  \`NULL\`.

## Value

A list containing two data frames (\`Biomass\` and \`Catch\`)
representing the detailed observation-level fit data filtered to the
requested species.

## Details

This function acts as a filter on top of \`rsim.fit.obj()\`. It runs the
base objective function in verbose mode to calculate the individual
negative log-likelihoods for all observation data. It then subsets the
resulting \`Biomass\` and \`Catch\` data frames, returning only the rows
where the \`Group\` matches the provided \`species\` vector.
