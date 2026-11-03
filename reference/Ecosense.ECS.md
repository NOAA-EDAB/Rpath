# Eastern Chukchi Sea Ecopath model

Rpath inputs for the eastern Chukchi Sea Ecopath model (52 biological
groups and 1 fleet).

## Usage

``` r
Ecosense.ECS
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

This is an aggregated version of the eastern Chukchi Sea Ecopath model
of Whitehouse and Aydin (2016).

## References

Whitehouse GA, Aydin KY (2016) Trophic structure of the eastern Chukchi
Sea: An updated mass balance food web model. U.S. Dep Commer, NOAA Tech.
Memo. NMFS-AFSC-318. doi:10.7289/V5/TM-AFSC-318

## Examples

``` r
# Balance the model
ecs_bal <- rpath(Ecosense.ECS)
```
