# Read Rpath parameters from .csv files

Creates an `Rpath.param` object from a series of .csv files.

## Usage

``` r
read.rpath.params(
  modfile,
  dietfile,
  pedfile = NA,
  stanzagroupfile = NA,
  stanzafile = NA
)
```

## Arguments

- modfile:

  file location of the flat file containing the *model* parameters.

- dietfile:

  file location of the flat file containing the *diet* parameters.

- pedfile:

  file location of the flat file containing the *pedgigree* parameters.

- stanzagroupfile:

  file location of the flat file containing the *group* parameters for
  multistanza groups. If not specified a blank stanza list will be
  created.

- stanzafile:

  file location of the flat file containing the individual stanza
  parameters for multistanza groups. If not specified a blank stanza
  list will be created.

## Value

Outputs an `Rpath.param` object that can be used for Rpath and
subsequently Rsim. (NOTE: This does function does not ensure data is
correct or in the correct locations.Run
[`check.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/check.rpath.params.md)
to ensure the appropriate columns are present).

## See also

Other Rpath functions:
[`MTI()`](https://noaa-edab.github.io/Rpath/reference/MTI.md),
[`check.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/check.rpath.params.md),
[`create.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/create.rpath.params.md),
[`frate.table()`](https://noaa-edab.github.io/Rpath/reference/frate.table.md),
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
