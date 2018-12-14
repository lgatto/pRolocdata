\name{itzhak2017}
\alias{itzhak2017}
\alias{itzhak2017markers}

\docType{data}

\title{
  Data from Itzhak et al. 2017
}

\description{

  The data from Itzhal et al. 2017 defines a spatial map for mouse
  primary neurons. The data are composed of 5 spatial maps, each
  containing 6 differential centrifugation fraction (as described in
  Itzhak et al. 2016, see \code{\link{itzhak2016}}).

  The annotatied marker proteins are available in the
  \code{itzhak2017markers} dataset.

  This data corresponds to the \emph{Mouse neuron Intensity and LFQ
  data} data from supplemental table S3 (Proteomic Mouse Neuron Data
  Generated in This Study, Related to Figure 5).

}

\usage{
data(itzhak2017)
data(itzhak2017markers)
}

\format{
  The data is an instance of class \code{MSnSet} from package
  \code{MSnbase}.
}


\references{

  Itzhak DN, Davies C, Tyanova S, Mishra A, William son J, Antrobus R,
  Cox J, Weekes MP, Borner GHH. A Mass Spectrometry-Based Approach for
  Mapping Protein Subcellular Localization Reveals the Spatial Proteome
  of Mouse Primary Neurons.  Cell Rep. 2017 Sep
  12;20(11):2706-2718. doi: 10.1016/j.celrep.2017.08.063. PubMed PMID:
  28903049; PubMed Central PMCID: PMC5775508.

}

\examples{
data(itzhak2017)
itzhak2017

## experimental design
table(pData(itzhak2017))

## the expression data
exprs(itzhak2017)[1:5, 1:5]

## abstract
abstract(itzhak2017)

## split data by samples
x <- split(itzhak2017, "map")

## plotting
library("pRoloc")
par(mfrow = c(2, 3))
for (i in 1:5)
    plot2D(x[[i]], main = paste("Map", i))
plot2D(itzhak2017, main = "All maps")
addLegend(itzhak2017, where = "bottomleft")
}

\keyword{datasets}
