# Gulf of Alaska (west/central) 1990s Ecopath model

Rpath inputs for the Gulf of Alaska (west/central) 1990s Ecopath model
(49 biological groups and 1 fleet).

## Usage

``` r
Ecosense.GOA
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

Whitehouse and Aydin 2020. Assessing the sensitivity of three Alaska
marine food webs to perturbations: an example of Ecosim simulations
using Rpath. https://doi.org/10.1016/j.ecolmodel.2020.109074

## Note

This is an aggregated version of the 1990s Gulf of Alaska Ecopath model
of Aydin et al. (2007) and does not include stanzas.

## References

Aydin KY, Gaichas S, Ortiz I, Kinzey D, Friday N (2007) A comparison of
the Bering Sea, Gulf of Alaska, and Aleutian Islands large marine
ecosystems through food web modeling. U.S. Dep. Commer., NOAA Tech.
Memo. NMFS-AFSC-178.

## Examples

``` r
# Balance the model
goa_bal <- rpath(Ecosense.GOA)
```
