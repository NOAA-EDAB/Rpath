# Creates an Rpath object from an EwE exported model (EIIXML format)

Parses an eiixml file exported using the Ecopath with Ecosim (EwE) GUI
into an unbalanced Rpath model object, which can then be balanced. This
function was tested on files exported by EwE version 6.7.

## Usage

``` r
create.rpath.from.eiixml(eiifile, verbose = FALSE)
```

## Arguments

- eiifile:

  Path to exported EwE XML file

- verbose:

  Logical. Use for debugging. If TRUE, prints out a list of parsed
  tables and number of rows read from the XML file

## Value

An Rpath object (list) with the following components:

- stanzas:

  Stanza parameters

- pedigree:

  Pedigree parameters

- diet:

  Diet parameters

- model:

  Model parameters

## See also

Other xml:
[`import.eiixml()`](https://noaa-edab.github.io/Rpath/reference/import.eiixml.md)

## Examples

``` r
# Import an eiixml file previously exported from the EwE GUI, check that the
# unbalanced model is functional, recalculate age-structured groups based on
# imported stanza parameters (if the imported model has stanzas), and balance
# the model.
eiixml <- system.file("extdata/xml/Western_Bering_Sea.eiixml", package = "Rpath")
rpath_object <- create.rpath.from.eiixml(eiixml)
check.rpath.params(rpath_object)
#> Rpath parameter file is functional. 
rpath_object <- rpath.stanzas(rpath_object)
rpath_balanced <- rpath(rpath_object)
```
