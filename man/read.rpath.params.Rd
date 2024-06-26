% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/param.R
\name{read.rpath.params}
\alias{read.rpath.params}
\title{Read Rpath parameters from .csv files}
\usage{
read.rpath.params(
  modfile,
  dietfile,
  pedfile = NA,
  stanzagroupfile = NA,
  stanzafile = NA
)
}
\arguments{
\item{modfile}{file location of the flat file containing the model parameters.}

\item{dietfile}{file location of the flat file containing the diet parameters.}

\item{pedfile}{file location of the flat file containg the pedgigree parameters.}

\item{stanzagroupfile}{file location of the flat file containing the group parameters
for multistanza groups.  If not specified a blank stanza list will be created.}

\item{stanzafile}{file location of the flat file containing the individual stanza 
parameters for multistanza groups.  If not specified a blank stanza list will 
be created.}
}
\value{
Outputs an Rpath.param object that can be used for Rpath and subsequently
 Rsim.  (NOTE: This does function does not ensure data is correct or in the 
 correct locations...run check.rpath.param to ensure the appropriate columns are
 present).
}
\description{
Creates an Rpath.param object from a series of .csv files.
}
\seealso{
Other Rpath functions: 
\code{\link{MTI}()},
\code{\link{adjust.forcing}()},
\code{\link{check.rpath.params}()},
\code{\link{create.rpath.params}()},
\code{\link{extract.node}()},
\code{\link{frate.table}()},
\code{\link{rpath}()},
\code{\link{rpath.groups}()},
\code{\link{rpath.stanzas}()},
\code{\link{stanzaplot}()},
\code{\link{webplot}()},
\code{\link{write.Rpath}()},
\code{\link{write.Rsim}()},
\code{\link{write.rpath.params}()}
}
\concept{Rpath functions}
