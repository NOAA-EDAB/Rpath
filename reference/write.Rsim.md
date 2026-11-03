# Write function for Rsim run

Outputs starting biomass, ending biomass, and catch to a .csv or .Rdata
file.

## Usage

``` r
write.Rsim(Rsim.output, file = NA)
```

## Arguments

- Rsim.output:

  object created by the
  [`rsim.run()`](https://noaa-edab.github.io/Rpath/reference/rsim.run.md)
  function.

- file:

  file name for resultant file. Need to specify ".csv" or ".RData".

## Value

Writes a list object or file with the start and end biomass and catch
per group from an Rpath.sim object.

- Group:

  Vector of group names

- StartBio:

  Vector of Biomass values at simulation start

- EndBio:

  Vector of Biomass values at simulation end

- BioES:

  Vector of ratio end:start Biomass values

- StartCatch:

  Vector of Catch values at simulation start

- EndCatch:

  Vector of Catch values at simulation end

- CatchES:

  Vector of ratio end:start Catch values

## See also

Other Rsim functions:
[`extract.node()`](https://noaa-edab.github.io/Rpath/reference/extract.node.md),
[`rsim.params()`](https://noaa-edab.github.io/Rpath/reference/rsim.params.md),
[`rsim.plot()`](https://noaa-edab.github.io/Rpath/reference/rsim.plot.md),
[`rsim.run()`](https://noaa-edab.github.io/Rpath/reference/rsim.run.md),
[`rsim.scenario()`](https://noaa-edab.github.io/Rpath/reference/rsim.scenario.md),
[`rsim.step()`](https://noaa-edab.github.io/Rpath/reference/rsim.step.md)

## Examples

``` r
# Read in Rpath parameter file and generate model object
Rpath <- rpath(AB.params)
# Create a 50 yr Rsim scenario
Rsim.scenario <- rsim.scenario(Rpath, AB.params, years = 1:50)
# Run the Rsim simulation
Rsim.output <- rsim.run(Rsim.scenario, method = "RK4", years = 1:50)
# Save summary output to a .csv file
write.Rsim(Rsim.output, file = "output.csv")

```
