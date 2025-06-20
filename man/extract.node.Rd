% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Extract.R
\name{extract.node}
\alias{extract.node}
\title{Extract node/group data from an Rsim run}
\usage{
extract.node(Rsim.output, group)
}
\arguments{
\item{Rsim.output}{Object generated by \code{rsim.run()}.}

\item{group}{Name of the node/group for which you want to extract data. Valid values are found in the `Group` field of the object
created from running \code{rpath()}}
}
\value{
Returns a list object with monthly node specific data
\item{Biomass}{Vector of monthly Biomass values}
\item{AnnualBiomass}{Vector of annual Biomass}
\item{TotalCatch}{Vector of monthly Catch values}
\item{AnnualTotalCatch}{Vector of annual Catch values}
\item{Landings}{Dataframe (n x l) of monthly Landings. (Each column is a fishery that Lands the `group`) }
\item{Discards}{Dataframe (n x d) of monthly Discards (Each column is a fishery that Discards the `group`)}
}
\description{
Extracts group/node specific information from a scenario run using \code{rsim.run()}
}
\examples{
# Read in Rpath parameter file and balance model
Rpath <- rpath(Ecosense.EBS)
# Create a 50 yr Rpath scenario
Rsim.scenario <- rsim.scenario(Rpath, Ecosense.EBS, years = 1:50)
# Run a scenario
Rsim.output <- rsim.run(Rsim.scenario, method = "RK4", years = 1:50)
# Extracts group (node) specific data from the simulation run
group <- extract.node(Rsim.output,Rpath$Group['Squids'])


}
\seealso{
Other Rsim functions: 
\code{\link{rsim.params}()},
\code{\link{rsim.plot}()},
\code{\link{rsim.run}()},
\code{\link{rsim.scenario}()},
\code{\link{rsim.step}()},
\code{\link{write.Rsim}()}
}
\concept{Rsim functions}
