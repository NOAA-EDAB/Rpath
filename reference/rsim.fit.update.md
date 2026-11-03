# Update Rpath Scene with Fitted Parameters

A function that applies fit variables to an Rsim scenario, returning a
new scenario object.

## Usage

``` r
rsim.fit.update(values, species, vartype, Rsim.scenario)
```

## Arguments

- values:

  A numeric vector of parameter values or modifiers. For a parameter of
  vartype \`"mzero"\`, the unit of the value is a mortality rate. Values
  for \`"predvul"\`, \`"preyvul"\`, and \`"ppvul"\` are log-scaled
  vulnerabilities.

- species:

  A character vector indicating the species/groups to modify (or string
  representations of numeric link indices for \`"ppvul"\`).

- vartype:

  A character vector indicating the types of parameters to modify.
  Currently supported parmaters are \`"mzero"\`, \`"predvul"\`,
  \`"preyvul"\`, and \`"ppvul"\`.

- Rsim.scenario:

  An \`Rsim\` scenario object containing the fitting data.

## Value

Modifies and returns the complete \`Rsim.scenario\` object with updated
vulnerability and mortality parameters.

## Details

This function translates the defined fitting variables into the
appropriate variables of an Rsim.scenario object. In particular,
applying an \`"mzero"\` parameter to a species replaces the
Rsim.scenario\$params\$MzeroMort paramter for that species.
\`"predvul"\`, \`"preyvul"\`, and \`"ppvul"\` are applied using the
formula for each predator/prey link
`Rsim.scenario$params$VV(predator, prey) = (1 + exp(log(Rsim.scenario$params$VV(predator, prey)-1) + preddiff + preydiff + ppdiff))`
