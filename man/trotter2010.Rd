\name{trotter20010}
\alias{trotter20010}
\alias{trotter20010steep}
\alias{trotter20010shallow}

\docType{data}

\title{
  LOPIT data sets used in Trotter et al. 2010.
}

\description{
  The two Arabidobsis LOPIT data sets \code{trotter2010shallow} and
  \code{trotter2010steep} have been used in Trotter et al. (2010)
  to illustrate improvement of sub-cellular resolution upon data
  fusion. The data have originally been published in Dunkley et
  al. (2006) and Sadowski et al. (2008), respectively.
}

\usage{
data(trotter2010shallow)
data(trotter2010steep)
}

\format{
  The data are instances of class \code{MSnSet} from package \code{MSnbase}.
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
 
  Sadowski PG, Groen AJ, Dupree P and Lilley KS.
  \emph{Sub-cellular localization of membrane proteins.}
  Proteomics 2008 8(19):3991-4011. PMID 18780351.

  Dunkley TP, Hester S, Shadforth IP, Runions J, Weimar T, Hanton SL, Griffin
  JL, Bessant C, Brandizzi F, Hawes C, Watson RB, Dupree P, Lilley
  KS. \emph{Mapping the Arabidopsis organelle proteome.} Proc Natl Acad
  Sci U S A. 2006 Apr 25;103(17):6518-23. Epub 2006 Apr 17. PubMed PMID:
  16618929; PubMed Central PMCID: PMC1458916.
}

\examples{
library(pRoloc)
## Replication of figure 4 from Trotter et al.
## individual data sets
data(trotter2010steep)
data(trotter2010shallow)
## combined data
combined <- combine(shallow, steep)
par(mfrow = c(1,3))
plot2D(shallow, fcol = "TAIR8", main = "Shallow")
plot2D(steep, fcol = "TAIR8", main = "Steep")
plot2D(combined, fcol = "TAIR8", main = "Combined")
addLegend(combined, where = "bottomleft", fcol = "TAIR8", bty = "n", ncol = 2)
}

\keyword{datasets}
