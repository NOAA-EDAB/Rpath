# Plot routine for Rsim runs

Plots the relative biomass of each group from an `Rsim.output` object.

## Usage

``` r
rsim.plot(Rsim.output, spname = "all", indplot = F)
```

## Arguments

- Rsim.output:

  R object containing the output from
  [`rsim.run()`](https://noaa-edab.github.io/Rpath/reference/rsim.run.md).

- spname:

  Vector of species names to be displayed. The default "all" will
  display all group names

- indplot:

  Logical value of whether to plot a single group or multiple groups

## Value

Creates a figure of relative biomass

## Details

For a detailed example of how to use this plot command please see
[`vignette("RunRsim")`](https://noaa-edab.github.io/Rpath/articles/RunRsim.md)

## See also

Other Rsim functions:
[`extract.node()`](https://noaa-edab.github.io/Rpath/reference/extract.node.md),
[`rsim.params()`](https://noaa-edab.github.io/Rpath/reference/rsim.params.md),
[`rsim.run()`](https://noaa-edab.github.io/Rpath/reference/rsim.run.md),
[`rsim.scenario()`](https://noaa-edab.github.io/Rpath/reference/rsim.scenario.md),
[`rsim.step()`](https://noaa-edab.github.io/Rpath/reference/rsim.step.md),
[`write.Rsim()`](https://noaa-edab.github.io/Rpath/reference/write.Rsim.md)
