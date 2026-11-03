# List Biomass Timeseries Groups and Sources

Extracts all unique groups, survey sources, and combined identifiers
from the biomass fitting data within an \`Rsim\` simulation scene.

## Usage

``` r
rsim.fit.list.bio.series(Rsim.scenario)
```

## Arguments

- Rsim.scenario:

  An \`Rsim\` simulation scene list object containing the fitting data
  setup.

## Value

A list containing three elements: \* \`groups\`: A character vector of
unique species/groups. \* \`sources\`: A character vector of unique
survey sources. \* \`all\`: A character vector of unique combinations
(\`"Group:Source"\`).

## Details

This helper function reads the \`Group\` and \`Source\` columns from the
\`scene\$fitting\$Biomass\` data frame. It returns a summary list
containing the unique functional groups, the unique survey sources, and
the unique combinations of the two formatted as \`"Group:Source"\`.
