# Write Rpath parameters to .csv files

Creates a series of .csv files from an `Rpath.params` object.

## Usage

``` r
write.rpath.params(Rpath.params, eco.name, path = "")
```

## Arguments

- Rpath.params:

  R object containing the Rpath parameters. Most likely this was created
  using
  [`create.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/create.rpath.params.md)
  or
  [`read.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/read.rpath.params.md).

- eco.name:

  ecosystem name that will be included in all the file names.

- path:

  location for the output files.

## Value

Outputs a series of .csv files named by the provided eco.name and the
parameters they represent. For example the model parameters will be
named "eco.name_model.csv".

## See also

Other Rpath functions:
[`MTI()`](https://noaa-edab.github.io/Rpath/reference/MTI.md),
[`check.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/check.rpath.params.md),
[`create.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/create.rpath.params.md),
[`frate.table()`](https://noaa-edab.github.io/Rpath/reference/frate.table.md),
[`read.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/read.rpath.params.md),
[`rpath()`](https://noaa-edab.github.io/Rpath/reference/rpath.md),
[`rpath.consumers()`](https://noaa-edab.github.io/Rpath/reference/rpath.consumers.md),
[`rpath.detrital()`](https://noaa-edab.github.io/Rpath/reference/rpath.detrital.md),
[`rpath.gears()`](https://noaa-edab.github.io/Rpath/reference/rpath.gears.md),
[`rpath.groups()`](https://noaa-edab.github.io/Rpath/reference/rpath.groups.md),
[`rpath.living()`](https://noaa-edab.github.io/Rpath/reference/rpath.living.md),
[`rpath.mixotrophs()`](https://noaa-edab.github.io/Rpath/reference/rpath.mixotrophs.md),
[`rpath.producers()`](https://noaa-edab.github.io/Rpath/reference/rpath.producers.md),
[`rpath.stanzas()`](https://noaa-edab.github.io/Rpath/reference/rpath.stanzas.md),
[`stanzaplot()`](https://noaa-edab.github.io/Rpath/reference/stanzaplot.md),
[`webplot()`](https://noaa-edab.github.io/Rpath/reference/webplot.md),
[`write.Rpath()`](https://noaa-edab.github.io/Rpath/reference/write.Rpath.md)
