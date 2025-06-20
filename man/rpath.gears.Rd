% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rpath_support.R
\name{rpath.gears}
\alias{rpath.gears}
\title{List of Rpath gears groups}
\usage{
rpath.gears(Rpath)
}
\arguments{
\item{Rpath}{Balanced Rpath model generated by \code{rpath()}}
}
\value{
Returns a string vector containing the names of Rpath functional groups
that are gear types
}
\description{
List of gears groups from Rpath object with species type = 3
}
\examples{
# Build the balanced Rpath model and parameter file by calling `rpath`
Rpath <- rpath(Ecosense.EBS)
# Run the groups command
gearGroups = rpath.gears(Rpath)
# Print out the first few gear type group names
head(gearGroups)

}
\seealso{
Other Rpath functions: 
\code{\link{MTI}()},
\code{\link{check.rpath.params}()},
\code{\link{create.rpath.params}()},
\code{\link{frate.table}()},
\code{\link{read.rpath.params}()},
\code{\link{rpath}()},
\code{\link{rpath.consumers}()},
\code{\link{rpath.detrital}()},
\code{\link{rpath.groups}()},
\code{\link{rpath.living}()},
\code{\link{rpath.mixotrophs}()},
\code{\link{rpath.producers}()},
\code{\link{rpath.stanzas}()},
\code{\link{stanzaplot}()},
\code{\link{webplot}()},
\code{\link{write.Rpath}()},
\code{\link{write.rpath.params}()}
}
\concept{Rpath functions}
