\name{lopimsSyn2}
\alias{lopimsSyn1}
\alias{lopimsSyn2}
\alias{lopimsSyn2_0frags}


\docType{data}

\title{
  LOPIMS data for the Synapter 2.0 paper
}

\description{
  TODO
}

\usage{
data("lopimsSyn1")
data("lopimsSyn2")
data("lopimsSyn2_0frags")
}

\format{

  These data are \code{MSnSet} instances, defined in the \code{MSnbase}
  package.

}

\examples{
data(lopimsSyn1)
data(lopimsSyn2)
data(lopimsSyn2_0frags)

## Visualisation
library("pRoloc")
par(mfrow = c(1, 3))
plot2D(lopimsSyn1, main = "Synapter 1", addLegend = "topleft")
plot2D(lopimsSyn2, main = "Synapter 2")
plot2D(lopimsSyn2_0frags, main = "Synapter 2 (0 fragments)")
}

\keyword{datasets}
