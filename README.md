

Rpath
=====

<!-- badges: start -->
[![gitleaks](https://github.com/NOAA-EDAB/Rpath/actions/workflows/secretScan.yml/badge.svg)](https://github.com/NOAA-EDAB/Rpath/actions/workflows/secretScan.yml)
[![pkgdown](https://github.com/NOAA-EDAB/Rpath/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/NOAA-EDAB/Rpath/actions/workflows/pkgdown.yaml)
[![R-CMD-check](https://github.com/NOAA-EDAB/Rpath/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/NOAA-EDAB/Rpath/actions/workflows/R-CMD-check.yaml)
[![tests](https://github.com/NOAA-EDAB/Rpath/actions/workflows/tests.yml/badge.svg)](https://github.com/NOAA-EDAB/Rpath/actions/workflows/R-CMD-check.yaml)
<!-- [![tests](https://github.com/NOAA-EDAB/Rpath/actions/workflows/tests.yml/badge.svg)](https://github.com/NOAA-EDAB/Rpath/actions/workflows/tests.yml) -->
<!-- badges: end -->

Rpath is an implementation (in R) of the ecosystem model Ecopath with Ecosim.

Rpath was developed and is maintained through a collaboration between the Alaska Fisheries Science Center and the Northeast Fisheries Science Center.  

## Installation

To install the package and build all of the vignettes locally

```
remotes::install_github("noaa-edab/Rpath",build_vignettes=TRUE)`
```

If you experience issues installing the package using `remotes` or don't need the vignettes locally then please use this alternative

```
pak::pak("noaa-edab/Rpath")
```





*This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.*

<img src="https://raw.githubusercontent.com/nmfs-general-modeling-tools/nmfspalette/main/man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" height="75" alt="NOAA Fisheries">

[U.S. Department of Commerce](https://www.commerce.gov/) | [National Oceanographic and Atmospheric Administration](https://www.noaa.gov) | [NOAA Fisheries](https://www.fisheries.noaa.gov/)
