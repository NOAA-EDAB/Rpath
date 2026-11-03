# Plot Full Biomass and Catch Fits for a Single Species Generates a multi-panel figure displaying the simulation fits for a single species across all available biomass observation sources, alongside its catch fitting data.

Plot Full Biomass and Catch Fits for a Single Species Generates a
multi-panel figure displaying the simulation fits for a single species
across all available biomass observation sources, alongside its catch
fitting data.

## Usage

``` r
rsim.plot.fitspecies(Rsim.scenario, Rsim.output, species)
```

## Arguments

- Rsim.scenario:

  An \`Rsim\` scenario object containing the fitting data.

- Rsim.output:

  An \`Rsim\` simulation output object resulting from a fitted run.

- species:

  A single character string representing the species name.

## Value

Invoked for its side effect of drawing a multi-panel plot to the active
graphics device. Returns \`NULL\`.

## Details

This function aggregates the biomass and catch visualization into a
single row of plots. It dynamically determines the number of biomass
data sources for the specified species and configures the graphics
layout (\`par(mfrow)\`) to fit all biomass plots plus one additional
panel for the catch plot. It calculates log-normal confidence intervals
for observations and overlays base reference lines.
