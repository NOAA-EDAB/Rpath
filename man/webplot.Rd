% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecopath_plot.R
\name{webplot}
\alias{webplot}
\title{Plot routine for Ecopath food web}
\usage{
webplot(
  Rpath.obj,
  eco.name = attr(Rpath.obj, "eco.name"),
  line.col = "grey",
  highlight = NULL,
  highlight.col = c("black", "red", "orange"),
  labels = FALSE,
  label.pos = NULL,
  label.num = FALSE,
  label.cex = 1,
  fleets = FALSE,
  type.col = "black",
  box.order = NULL
)
}
\arguments{
\item{Rpath.obj}{Rpath model created by \code{rpath()}.}

\item{eco.name}{Optional name of the ecosystem.  Default is the `eco.name` attribute from the
rpath object created from running \code{rpath()}.}

\item{line.col}{The color of the lines between nodes of the food web.}

\item{highlight}{Group number or name to highlight the connections of that group. Valid values are found in the 
`Group` field of the object created from running \code{rpath()}.}

\item{highlight.col}{Color of the connections to the highlighted group, vector of length 3.
Defaults to black for prey, red for predator, orange for fleet.}

\item{labels}{Logical whether or not to display group names.  If \code{TRUE} and \code{label.pos} = NULL, no 
points will be plotted, just label names.}

\item{label.pos}{A position specifier for the labels.  Values of 1, 2, 3, 4
indicate positions below, to the left of, above, and to the right of the points, respectively. A null 
value will cause the labels to be plotted without the points, assuming that \code{labels} = TRUE.}

\item{label.num}{Logical value indication whether group numbers should be used for labels 
instead of names.}

\item{label.cex}{Numeric value of the relative size of the labels within the plot.}

\item{fleets}{Logical value indicating whether or not to include fishing fleets in the food web plot.}

\item{type.col}{The color of the points corresponding to the types of the group.  Must be 
of length 1 or 4.  Color order will be consumers, primary producers, detrital, and fleet groups.}

\item{box.order}{Vector of box numbers to change the default plot order.  Must include all box numbers}
}
\value{
Returns a plot visualization of the food web.
}
\description{
Plots the food web associated with an Rpath object.
}
\examples{
# Read in Rpath parameter file, generate and name model object
Rpath.obj <- rpath(AB.params, eco.name = "Anchovy Bay")
# Plot food web diagram with all groups labeled, including fleets
webplot(Rpath.obj, labels = TRUE, fleets = TRUE)
# Plot food web diagram without labels, highlighting connections of cod group
webplot(Rpath.obj, highlight = "cod",fleets = TRUE)

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
\code{\link{rpath.gears}()},
\code{\link{rpath.groups}()},
\code{\link{rpath.living}()},
\code{\link{rpath.mixotrophs}()},
\code{\link{rpath.producers}()},
\code{\link{rpath.stanzas}()},
\code{\link{stanzaplot}()},
\code{\link{write.Rpath}()},
\code{\link{write.rpath.params}()}
}
\concept{Rpath functions}
