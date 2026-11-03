# Prepare food web model to plot

Prepare and summarize food web model for plotting by generating x-values
for node locations and creating a data frame of all trophic connections.
This function could also be used to build custom plotting functions.

## Usage

``` r
summarize.for.webplot(Rpath.obj, fleets = FALSE, box.order = NULL)
```

## Arguments

- Rpath.obj:

  Rpath model created by the
  [`rpath()`](https://noaa-edab.github.io/Rpath/reference/rpath.md)
  function.

- fleets:

  Logical value indicating whether or not to include fishing fleets in
  the food web.

- box.order:

  Vector of box numbers to change the default plot order. Must include
  all box numbers

## Value

Returns a list of length 2. The first element, `pointmap` is a data
table of the functional group names, group numbers, type (0- living, 1-
producer, 2- detritus, 3- fleet), trophic level, biomass, and `x.pos`
which determines the x coordinate of the node when plotting. The second
element is a data table with one row for each trophic connection, and
columns for predator and prey trophic levels, predator and prey `x.pos`
from `pointmap`, and predator and prey group numbers (`id`).
