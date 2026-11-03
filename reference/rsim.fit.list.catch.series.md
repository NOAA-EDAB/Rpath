# List Catch Timeseries Groups

Extracts the unique species or functional groups that have observation
data within the catch fitting matrix of an \`Rsim\` simulation scenario.

## Usage

``` r
rsim.fit.list.catch.series(Rsim.scenario)
```

## Arguments

- Rsim.scenario:

  An \`Rsim\` scenario list object containing fitting data setup.

## Value

A list containing a single element, \`groups\`, which holds a vector of
the unique group names present in the catch data.

## Details

This is a helper function designed to quickly identify which groups have
empirical catch data loaded into the scenario. It reads the \`Group\`
column from the \`Rsim.scenario\$fitting\$Catch\` data frame and returns
the unique values.
