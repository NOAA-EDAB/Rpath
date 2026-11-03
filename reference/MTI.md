# Quantify mixed trophic impacts

Uses a static `Rpath` model and creates a matrix of mixed trophic
impacts.

## Usage

``` r
MTI(Rpath, Rpath.params)
```

## Arguments

- Rpath:

  Rpath R object containing a static `Rpath` model.

- Rpath.params:

  R object containing the `Rpath` model parameters. This is generated
  either by the
  [`create.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/create.rpath.params.md)
  or
  [`read.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/read.rpath.params.md)
  functions.

## Value

Returns a matrix of mixed trophic impacts.

## See also

Other Rpath functions:
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
[`write.Rpath()`](https://noaa-edab.github.io/Rpath/reference/write.Rpath.md),
[`write.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/write.rpath.params.md)
