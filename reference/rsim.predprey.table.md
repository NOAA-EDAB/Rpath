# Extract Predator-Prey Interaction Table

Generates a data frame of predator-prey linkages and their associated
foraging parameters from an \`Rpath\` simulation scenario.

## Usage

``` r
rsim.predprey.table(Rsim.scenario, pred = NULL, prey = NULL)
```

## Arguments

- Rsim.scenario:

  An \`Rsim\` scenario object containing the fitting data.

- pred:

  An optional character string or vector of predator names. If not
  provided, the function will include all model predators in its output.

- prey:

  An optional character string or vector of prey names. If not provided,
  the function will include all model prey in its output.

## Value

A \`data.frame\` containing the columns:

- `from`, prey name,

- `to` predator name,

- `QQ` base consumption rate of the predator/prey link,

- `VV` vulnerability parameter of the predator/prey link,

- `DD` handling time parameter of the predator/prey link

## Details

This function outputs a data frame of predator/prey link parameters
contained in an Rsim scenario for a given set of predators and prey.
