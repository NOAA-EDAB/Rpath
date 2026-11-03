# Plot Biomass and Catch Fitting Results

Generates multi-panel diagnostic plots comparing simulation run outputs
against observed biomass and catch data across a specified list of
species.

## Usage

``` r
rsim.plot.fit(Rsim.scenario, Rsim.output, species.list, scene.name = NULL)
```

## Arguments

- Rsim.scenario:

  An \`Rsim\` scenario object containing the fitting data.

- Rsim.output:

  An \`Rsim\` simulation output object resulting from a fitted run.

- species.list:

  A character vector of species names to include in the plots.

- scene.name:

  An optional character string used to label the main title of the
  plots. Defaults to \`NULL\` (resulting in a blank space).

## Value

Invoked for its side effect of drawing plots to external graphics
devices. Returns \`NULL\`.

## Details

This function opens two separate graphic windows. The first displays
biomass fits, iterating through all available data sources for each
specified species. The second window displays catch fits. Both plots
include an overarching title that dynamically calculates and reports the
negative log-likelihood (NLL) of the fit.
