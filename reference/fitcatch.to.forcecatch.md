# Convert Fitting Catch to Forced Catch in Rsim Scenario

Performs 3 steps: Sets effort to 0, sets default F, sets catch.
Integrates observed catch data into an \`Rpath\` simulation scene by
zeroing out effort, setting a baseline fishing mortality rate, and
forcing specific catches for observation years.

## Usage

``` r
fitcatch.to.forcecatch(Rsim.scenario, Rpath)
```

## Arguments

- Rsim.scenario:

  An \`Rsim\` scenario object containing the observed fishing data.

- Rpath:

  A balanced \`Rpath\` model object containing base biomass and
  landings.

## Value

Modifies and returns the \`Rsim.scenario\` object with updated
\`ForcedEffort\`, \`ForcedFRate\`, and \`ForcedCatch\` matrices.

## Details

The function modifies the fishing parameters of the simulation by: 1.
Setting all gear effort (\`ForcedEffort\`) to zero. 2. Calculating a
baseline Fishing Mortality rate (F) for all living and detrital groups
based on the balanced \`Rpath\` model (\`(Landings + Discards) /
Biomass\`). 3. Assigning this calculated F to the \`ForcedFRate\` matrix
for all years. 4. Overriding specific years and groups with empirical
catch data (\`fitting\$Catch\$obs\`). For these targeted overrides, the
\`ForcedFRate\` is set to 0 to allow the forced catch to drive the
simulation.
