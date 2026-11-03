# Ecopath module of Rpath

Performs initial mass balance using a `Rpath.params()` file

## Usage

``` r
rpath(Rpath.params, eco.name = NA, eco.area = 1)
```

## Arguments

- Rpath.params:

  R object containing the parameters needed to create a `Rpath` model.
  This is generated either by the
  [`create.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/create.rpath.params.md)
  or
  [`read.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/read.rpath.params.md)
  functions.

- eco.name:

  Optional name of the ecosystem which becomes an attribute of the
  `Rpath` model.

- eco.area:

  Optional area of the ecosystem which becomes an attribute of the
  `Rpath` model.

## Value

Returns a static `Rpath` model that can be supplied to the
`rsim.scenario` function.

## See also

Other Rpath functions:
[`MTI()`](https://noaa-edab.github.io/Rpath/reference/MTI.md),
[`check.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/check.rpath.params.md),
[`create.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/create.rpath.params.md),
[`frate.table()`](https://noaa-edab.github.io/Rpath/reference/frate.table.md),
[`read.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/read.rpath.params.md),
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
