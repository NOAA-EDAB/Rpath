# Set Biomass Fitting Data Weights

Updates the statistical weighting (\`wt\`) for specific biomass
observation timeseries prior to running a model fit in an \`Rsim\`
simulation scene. Matches are based on both the functional group and the
survey source.

## Usage

``` r
rsim.fit.set.bio.wt(Rsim.scenario, group, source, wt)
```

## Arguments

- Rsim.scenario:

  An Rsim scenario object containing the fitting data.

- group:

  A character string or vector of group names to update.

- source:

  A character string or vector of survey source names to update.

- wt:

  A numeric value to assign as the new weight.

## Value

Modifies and returns the \`Rsim.scenario\`object with the updated
weights applied to \`Rsim.scenario\$fitting\$Biomass\`.

## Details

This function locates all observation rows in the
\`Rsim.scenario\$fitting\$Biomass\` data frame that match BOTH the
specified \`group\` and \`source\` names. It replaces their \`"wt"\`
column with the provided \`wt\` value. If the specific combination is
not found in the biomass data, it issues a warning and returns the scene
unmodified.
