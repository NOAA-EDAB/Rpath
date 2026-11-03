# Apply Fitted Parameters to Rpath Simulation Scene Helper function that maps a vector of optimized parameter modifications back into the \`Rpath\` simulation parameters object (\`scene.params\`). Since this returns a scenario params object instead of a full scenario, it is primarily intended for internal use.

Apply Fitted Parameters to Rpath Simulation Scene Helper function that
maps a vector of optimized parameter modifications back into the
\`Rpath\` simulation parameters object (\`scene.params\`). Since this
returns a scenario params object instead of a full scenario, it is
primarily intended for internal use.

## Usage

``` r
rsim.fit.apply(values, species, vartype, scene.params)
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

- scene.params:

  An \`Rpath\` simulation parameters list (usually \`scene\$params\`).

## Value

Modifies and returns the \`scene.params\` object with updated
\`MzeroMort\` and \`VV\` matrices.

## Details

This function processes a vector of parameter modifications (\`values\`)
based on their assigned \`vartype\`. It handles four parameter types: \*
\`"mzero"\`: Modifies non-predation mortality directly. \*
\`"predvul"\`: Modifies vulnerability from the predator's perspective.
\* \`"preyvul"\`: Modifies vulnerability from the prey's perspective. \*
\`"ppvul"\`: Modifies vulnerability for a single, specific predator-prey
link. For vulnerabilities, the modifications are applied to the base
vulnerability (\`VV\`) using a logarithmic transformation: \`VV_new =
1 + exp(log(VV_base - 1) + modifications)\`.
