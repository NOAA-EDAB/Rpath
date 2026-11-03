# Set Catch Fitting Data Weights

Updates the statistical weighting (\`wt\`) for specific catch
observation timeseries prior to running a model fit in an \`Rsim\`
simulation scene. Matches are based on both the functional group and the
survey source.

## Usage

``` r
rsim.fit.set.catch.wt(Rsim.scenario, group, wt)
```

## Arguments

- Rsim.scenario:

  An \`Rsim\` scenario object containing the fitting data.

- group:

  A character string or vector of group names to update.

- wt:

  A numeric value to assign as the new weight.

## Value

Modifies and returns the \`Rsim.scenario\`object with the updated
weights applied to \`Rsim.scenario\$fitting\$Biomass\`.

## Details

This function locates all observation rows in the
\`Rsim.scenario\$fitting\$Catch\` data frame that match BOTH the
specified \`group\` and \`source\` names. It replaces their \`"wt"\`
column with the provided \`wt\` value. If the specific combination is
not found in the catch data, it issues a warning and returns the scene
unmodified.
