# Check Rpath parameter files

Logical check that the parameter files are filled out correctly, i.e.
data is entered where it is expected.

## Usage

``` r
check.rpath.params(Rpath.params)
```

## Arguments

- Rpath.params:

  R object containing the parameters needed to create a `Rpath` model.
  This is generated either by the
  [`create.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/create.rpath.params.md)
  or
  [`read.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/read.rpath.params.md)
  functions.

## Value

Checks Rpath parameter files for consistency. An error message will be
produced if one of the logical checks fails. Checks include:
verification that all types are represented (e.g. consumer, producer,
detrital, and fleet); check if input parameters were entered correctly;
check if Diet columns sum to 1. (NOTE: This does not ensure data is
correct just that it is in the right places).

## See also

Other Rpath functions:
[`MTI()`](https://noaa-edab.github.io/Rpath/reference/MTI.md),
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
