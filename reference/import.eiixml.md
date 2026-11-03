# Reads in EwE exported XML file and parses into data frames

Parses an eiixml file exported using the Ecopath with Ecosim (EwE) GUI
into a list of data frames, one frame for each table in the exported XML
file. This function is usually called by the function
[`create.rpath.from.eiixml()`](https://noaa-edab.github.io/Rpath/reference/create.rpath.from.eiixml.md)
that converts these tables into an unbalanced rpath model object.
However import.xml can be used on its own to examine the full set of
tables exported by EwE, including tables not currently imported into
Rpath objects, such as Ecosim runs or model metadata. This function was
tested on files exported by EwE version 6.7.

## Usage

``` r
import.eiixml(eiifile, verbose = FALSE)
```

## Arguments

- eiifile:

  Path to exported EwE XML file

- verbose:

  Logical. Use for debugging. If TRUE, prints out a list of parsed
  tables and number of rows read from the XML file

## Value

A list of data frames, one data frame for each node (exported EwE table)
in the XML file. Each table has the naming convention ewe\_\[table
name\] where \[table name\] is the name of the table provided by EwE.

## See also

Other xml:
[`create.rpath.from.eiixml()`](https://noaa-edab.github.io/Rpath/reference/create.rpath.from.eiixml.md)

## Examples

``` r
# Import an eiixml file previously exported from the EwE GUI into a list of
# data frames containing the model data
eiixml <- system.file("extdata/xml","Western_Bering_Sea.eiixml", package = "Rpath")
xml_data <- import.eiixml(eiixml)

```
