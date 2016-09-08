\name{trotter20010}
\alias{trotter2010}
\alias{trotter2010steep}
\alias{trotter2010shallow}

\docType{data}

\title{
  LOPIT data sets used in Trotter et al. (2010)
}

\description{
  The two Arabidobsis LOPIT data sets \code{trotter2010shallow} and
  \code{trotter2010steep} have been used in Trotter et al. (2010)
  to illustrate improvement of sub-cellular resolution upon data
  fusion. The data have originally been published in Dunkley et
  al. (2006) and Sadowski et al. (2008), respectively.

  The feature metadata contains the cellular compartment from TAIR8 and
  the pRoloc Arabidopsis markers (see
  \code{\link[pRoloc]{pRolocmarkers}}).
  
}

\usage{
data(trotter2010)
data(trotter2010shallow)
data(trotter2010steep)
}

\format{
  The data are instances of class \code{MSnSet} from package
  \code{MSnbase}. \code{trotter2010} corresponds to the combined steep
  and shallow data.
}

\source{
  Supporting information available on 
  \url{http://onlinelibrary.wiley.com/doi/10.1002/pmic.201000359/abstract}
}

\references{

  Trotter MWB, Sadowski PG, Dunkley TPJ, Groen AJ and Lilley KS.
  \emph{Improved sub-cellular resolution via simultaneous analysis of
  organelle proteomics data across varied experimental conditions.}
  Proteomics 2010 10(23):4213-4219. PMID 21058340.
 
  Sadowski PG, Groen AJ, Dupree P and Lilley KS.  \emph{Sub-cellular
  localization of membrane proteins.}  Proteomics 2008
  8(19):3991-4011. PMID 18780351.

}

\examples{
library(pRoloc)
## Replication of figure 4 from Trotter et al.
## individual data sets
data(trotter2010)
data(trotter2010steep)
data(trotter2010shallow)

par(mfrow = c(2,3))
plot2D(trotter2010shallow, fcol = "TAIR8", main = "Shallow (TAIR8)")
plot2D(trotter2010steep, fcol = "TAIR8", main = "Steep (TAIR8)")
plot2D(trotter2010, fcol = "TAIR8", main = "Combined (TAIR8)")
addLegend(trotter2010, where = "bottomleft", fcol = "TAIR8", ncol = 2)
plot2D(trotter2010shallow, main = "Shallow (markers)")
plot2D(trotter2010steep, main = "Steep (markers)")
plot2D(trotter2010, main = "Combined (markers)")
addLegend(trotter2010, where = "bottomleft", ncol = 2)

}

\keyword{datasets}
