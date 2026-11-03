# AB.params

Anchovy Bay (`AB.params`) is a sample ecosystem frequently used as an
example ecosystem with EwE (see Christensen and Walters, 2024).

## Usage

``` r
AB.params
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

## Source

Christensen V, Walters CJ (2024) Ecosystem Modelling with EwE. The
University of British Columbia, Vancouver, BC. doi:10.14288/24d7-ab68

## Examples

``` r
# Balance the model
AB_bal <- rpath(AB.params)
 
```
