# Set Catchability (q) for Biomass Fitting Data

Updates the catchability coefficient (\`initial_q\`) and timeseries
\`Type\` for specific biomass observations within an \`Rpath\`
simulation scene.

## Usage

``` r
rsim.fit.set.q(
  Rsim.scenario,
  group,
  source,
  q = NULL,
  years = NULL,
  type = NULL
)
```

## Arguments

- Rsim.scenario:

  An Rsim scenario object containing the fitting data.

- group:

  A character string (or vector) of the functional group(s) to update.

- source:

  A character string (or vector) of the survey source(s) to update.

- q:

  An optional numeric value to directly set the catchability
  coefficient.

- years:

  An optional numeric or character vector of years used to calculate
  \`q\` empirically.

- type:

  An optional character string to set the timeseries type (defaults to
  \`"fixed"\`).

## Value

Modifies and returns the \`Rsim.scenario\` (scene) object with updated
\`initial_q\` and \`Type\` columns in the \`scene\$fitting\$Biomass\`
data frame.

## Details

This function operates in three distinct modes depending on the
arguments provided: 1. \*\*Default Reset (\`q\` and \`years\` are
\`NULL\`):\*\* Resets \`initial_q\` to \`1.0\` and sets \`Type\` to
\`"index"\`. 2. \*\*Direct Assignment (\`q\` is provided):\*\* Assigns
the provided numeric \`q\` to \`initial_q\` and assigns the specified
\`type\` (which defaults to \`"fixed"\`). 3. \*\*Empirical Calculation
(\`years\` is provided, \`q\` is \`NULL\`):\*\* Calculates \`q\` by
taking the mean observed biomass (\`Value \* Scale\`) over the specified
\`years\` and dividing it by the model's base reference biomass
(\`B_BaseRef\`) for the specified \`group\`.
