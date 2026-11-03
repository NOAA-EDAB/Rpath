# Plot routine for Ecopath food web

Plots the food web associated with an Rpath object.

## Usage

``` r
webplot(
  Rpath.obj,
  eco.name = attr(Rpath.obj, "eco.name"),
  line.col = "grey",
  highlight = NULL,
  highlight.col = c("black", "red", "orange"),
  labels = FALSE,
  label.pos = NULL,
  label.num = FALSE,
  label.cex = 1,
  fleets = FALSE,
  type.col = "black",
  box.order = NULL
)
```

## Arguments

- Rpath.obj:

  Rpath model created by
  [`rpath()`](https://noaa-edab.github.io/Rpath/reference/rpath.md).

- eco.name:

  Optional name of the ecosystem. Default is the \`eco.name\` attribute
  from the rpath object created from running
  [`rpath()`](https://noaa-edab.github.io/Rpath/reference/rpath.md).

- line.col:

  The color of the lines between nodes of the food web.

- highlight:

  Group number or name to highlight the connections of that group. Valid
  values are found in the \`Group\` field of the object created from
  running
  [`rpath()`](https://noaa-edab.github.io/Rpath/reference/rpath.md).

- highlight.col:

  Color of the connections to the highlighted group, vector of length 3.
  Defaults to black for prey, red for predator, orange for fleet.

- labels:

  Logical whether or not to display group names. If `TRUE` and
  `label.pos` = NULL, no points will be plotted, just label names.

- label.pos:

  A position specifier for the labels. Values of 1, 2, 3, 4 indicate
  positions below, to the left of, above, and to the right of the
  points, respectively. A null value will cause the labels to be plotted
  without the points, assuming that `labels` = TRUE.

- label.num:

  Logical value indication whether group numbers should be used for
  labels instead of names.

- label.cex:

  Numeric value of the relative size of the labels within the plot.

- fleets:

  Logical value indicating whether or not to include fishing fleets in
  the food web plot.

- type.col:

  The color of the points corresponding to the types of the group. Must
  be of length 1 or 4. Color order will be consumers, primary producers,
  detrital, and fleet groups.

- box.order:

  Vector of box numbers to change the default plot order. Must include
  all box numbers

## Value

Returns a plot visualization of the food web.

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
[`write.Rpath()`](https://noaa-edab.github.io/Rpath/reference/write.Rpath.md),
[`write.rpath.params()`](https://noaa-edab.github.io/Rpath/reference/write.rpath.params.md)

## Examples

``` r
# Read in Rpath parameter file, generate and name model object
Rpath.obj <- rpath(AB.params, eco.name = "Anchovy Bay")
# Plot food web diagram with all groups labeled, including fleets
webplot(Rpath.obj, labels = TRUE, fleets = TRUE)

# Plot food web diagram without labels, highlighting connections of cod group
webplot(Rpath.obj, highlight = "cod",fleets = TRUE)

```
