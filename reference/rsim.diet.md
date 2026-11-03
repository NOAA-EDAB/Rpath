# Output consumption of each prey by an individual predator group from Rsim run

Extracts a matrix of consumption (in model biomass units) of each prey
over time by a particular predator from the Rsim model output of
[`rsim.run()`](https://noaa-edab.github.io/Rpath/reference/rsim.run.md).

## Usage

``` r
rsim.diet(Rsim.output, group)
```

## Arguments

- Rsim.output:

  R object containing the output from `rsim.run`.

- group:

  Predator group from the `Rpath` model that is of interest

## Value

a numeric matrix of annual (rows) consumption of each prey group
(columns)

## Examples

``` r
# Read in Rpath parameter file and generate model object
Rpath <- rpath(AB.params)
# Create a 50 yr Rsim scenario
Rsim.scenario <- rsim.scenario(Rpath, AB.params, years = 1:50)
# Run the Rsim simulation
Rsim.output <- rsim.run(Rsim.scenario, method = "RK4", years = 1:50)
# Extract a preadator's consumption of each prey over the model run
Rsim.diet.whales <- rsim.diet(Rsim.output, group = "whales")
```
