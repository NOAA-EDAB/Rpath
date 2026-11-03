# Creates a shell for the Rpath.params object

Creates a shell of the `Rpath.params` list object which contains the
*model*, *diet*, *multistanza*, and *pedigree* parameters.

## Usage

``` r
create.rpath.params(group, type, stgroup = NA)
```

## Arguments

- group:

  Vector of group names.

- type:

  Numeric vector of group type. Living = 0, Producer = 1, Detritus = 2,
  Fleet = 3.

- stgroup:

  Vector of multistanza group names. Include NA for non-stanza groups.

## Value

Outputs a list object of `Rpath.params` which are populated with values
of NA or logical default values. Values can then be filled in using R.
Use
[`check.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/check.rpath.params.md)
to ensure parameter files are filled out correctly (NOTE: This does not
ensure data is correct just that it is in the right places).

## See also

Other Rpath functions:
[`MTI()`](https://noaa-edab.github.io/Rpath/reference/MTI.md),
[`check.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/check.rpath.params.md),
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
