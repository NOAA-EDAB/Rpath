# Output biomass removed by each predator for an individual prey group from Rsim run

Extracts a matrix of biomass removed (in model biomass units) for a
particular prey by each of its predators from the Rsim model output of
[`rsim.run()`](https://noaa-edab.github.io/Rpath/reference/rsim.run.md).

## Usage

``` r
rsim.mort(Rsim.output, group)
```

## Arguments

- Rsim.output:

  R object containing the output from `rsim.run`.

- group:

  Predator group from the `Rpath` model that is of interest

## Value

a numeric matrix of annual (rows) biomass removed by each predator group
(columns)

## Examples

``` r
# Read in Rpath parameter file and generate model object
Rpath <- rpath(AB.params)
# Create a 50 yr Rsim scenario
Rsim.scenario <- rsim.scenario(Rpath, AB.params, years = 1:50)
# Run the Rsim simulation
Rsim.output <- rsim.run(Rsim.scenario, method = "RK4", years = 1:50)
# Extract a prey's biomass loss from each predator over the model run
Rsim.mort.anchovy <- rsim.mort(Rsim.output, group = "anchovy")
```
