# List Biomass Timeseries Groups by Survey Source

Extracts and groups the unique species or functional groups associated
with each survey source in the biomass fitting data of an \`Rsim\`
simulation scene.

## Usage

``` r
rsim.fit.list.bio.bysurvey(Rsim.scenario)
```

## Arguments

- Rsim.scenario:

  An \`Rpath\` simulation scene list object containing the fitting data
  setup.

## Value

A named list where each element name corresponds to a unique survey
\`Source\`, and the contents are character vectors of the \`Group\`
names associated with that source.

## Details

This helper function reads the \`Source\` and \`Group\` columns from the
\`scene\$fitting\$Biomass\` data frame. It iterates through every unique
survey source and compiles a list of all groups (species) that have
observation data attributed to that specific source.
