# Read and Process Biomass Fitting Data

Reads, cleans, and formats biomass observation data to prepare it for
model fitting within an \`Rsim\` scenario object.

## Usage

``` r
read.fitting.biomass(Rsim.scenario, biomass_data)
```

## Arguments

- Rsim.scenario:

  An \`Rsim\` scenario object containing simulation parameters and
  setup.

- biomass_data:

  A \`data.frame\` containing biomass observation data.

## Value

Modifies and returns a new Rsim scenario object, attaching the cleaned
data frame to \`Rsim.scenario\$fitting\$Biomass\` with additional
calculated columns (\`obs\`, \`sd\`, \`initial_q\`, \`wt\`).

## Details

Observation data must be provided in the \`biomass_data\` data frame.
\`biomass_data\` must contain the following columns: \`Group\`,
\`Year\`, \`Value\`, \`Scale\`, \`Source\`, one of either \`Stdev\` or
\`CV\` (or both), and optionally \`Type\`. If both \`Stdev\` and \`CV\`
are provided, \`Stdev\` is used for each data point unless it is missing
(NA) in which case \`CV\` is used.
