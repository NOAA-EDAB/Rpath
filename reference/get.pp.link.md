# Get Predator-Prey Link Index

Retrieves the numeric index of a specific predator-prey interaction from
the model's internal diet network structure.

## Usage

``` r
get.pp.link(Rsim.scenario, predator, prey)
```

## Arguments

- Rsim.scenario:

  An \`Rsim\` scenario object containing the fitting data..

- predator:

  A character string specifying the name of the predator group.

- prey:

  A character string specifying the name of the prey group.

## Value

A numeric scalar (or vector) representing the internal link index
associated with the specified predator-prey interaction.

## Details

This helper function scans the internal C++ 0-indexed arrays
(\`PreyFrom\` and \`PreyTo\`), converts them to 1-based R indices to
look up the corresponding species names in \`spname\`, and returns the
link index (or indices) where both the predator and prey match the
requested strings.
