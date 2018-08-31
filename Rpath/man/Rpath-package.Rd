\name{Rpath-package}
\alias{Rpath-package}
\alias{Rpath}
\docType{package}
\title{
R implementation of the Ecopath with Ecosim model
}
\description{
This package implements the core functionality of the popular Ecopath with Ecosim
mass balance model.  Ecopath is a static snap shot of the energy flow within an
ecosystem.  Ecosim is a dynamic simulation using the Ecopath values as a starting
point.
}
\details{
\tabular{ll}{
Package: \tab Rpath\cr
Type: \tab Package\cr
Version: \tab 0.0.0.9001\cr
Date: \tab 2016-07-15\cr
License: \tab GPL (>= 2)\cr
}
There are two main components to the package: rpath and rsim.  The function rpath
uses a parameter object create by create.rpath.params to conduct the mass balance
of the system.  The rsim component is carried out using a series of functions.
The first function is rsim.scenario which converts the rpath output to rates 
appropriate for running the simulation.  You can adjust the base scenario through
a series of adjust functions (adjust.scenario, adjust.fishing, adjust.forcing).
Finally you run the simulation with the function rsim.run.  Visualizations of the
foodweb and rsim results can be made using the functions webplot and rsim.plot 
respectively.
}
\author{
Sean M. Lucey, Sarah K. Gaichas, and Kerim Y. Aydin

Maintainer: Sean M. Lucey <Sean.Lucey@NOAA.gov>
}
\references{
~~ Literature or other references for background information ~~
}

\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
%% ~~ simple examples of the most important functions ~~
}
