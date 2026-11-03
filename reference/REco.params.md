# REco.params

Rpath inputs for a tutorial Rpath ecosystem.

## Usage

``` r
REco.params
```

## Format

An unbalanced Rpath model object that contains:

- model:

  A data.table containing Ecopath unbalanced ecosystem parameters (base
  parameters and fisheries).

- diet:

  A data.table containing the Ecopath model's diet matrix.

- stanzas:

  Parameters for multistanza groups.

- pedigree:

  A data.table containing the data quality (pedigree) for the Ecopath
  model.

## References

Lucey SM, Gaichas SK, Aydin KY (2020) Conducting reproducible ecosystem
modeling using the open source mass balance model Rpath. Ecol Model
427:11. doi:10.1016/j.ecolmodel.2020.109057

## See also

[`vignette("ModelSetup", package = "Rpath")`](https://noaa-edab.github.io/Rpath/articles/ModelSetup.md)
for an example of how to build this tutorial ecosystem.

## Examples

``` r
# Balance the model
REco_bal <- rpath(REco.params)
```
