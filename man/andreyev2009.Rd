\name{andreyev2009}
\alias{andreyev2009}
\alias{andreyev2009rest}
\alias{andreyev2009activ}
\docType{data}

\title{

  Six sub-cellular fraction data from mouse macrophage-like RAW264.7
  cells from Andreyev et al. (2009)

}

\description{

  Data from Andreyev AY, Shen Z, Guan Z, Ryan A, Fahy E, Subramaniam S,
  Raetz CR, Briggs S, Dennis EA. Application of proteomic marker
  ensembles to subcellular organelle identification. Mol Cell
  Proteomics. 2010
  Feb;9(2):388-402. DOI:\url{http://dx.doi.org/10.1074/mcp.M900432-MCP200}. PubMed
  PMID:19884172; PubMed Central PMCID:PMC2830848.

  The 6 subcellular fractions are nuclei, mitochondria, cytoplasm,
  endoplasmic reticulum, plasma membrane and dense microsomal fraction.
  
}

\usage{
data("andreyev2009")
data("andreyev2009rest")
data("andreyev2009activ")
}


\details{

  \code{andreyev2009} is the full data where missing values were
  replaced by 0. \code{andreyev2009rest} and \code{andreyev2009activ}
  contain the resting (control) and Kdo2-lipid A-treated (activated)
  data respectively, which have been normalised (each reporter intensity
  was normalised by the sum over all replicates).
  
}

\source{

  These data were generated from supplemental tables S1 (quantitative
  data) and 2 (organelle markers)
  (\url{http://www.mcponline.org/content/9/2/388/suppl/DC1}). See
  \code{inst/scripts/andreyev2009.R} for details.
  
}

\examples{
data(andreyev2009rest)
data(andreyev2009activ)

library("pRoloc")
par(mfrow = c(1, 2))
plot2D(andreyev2009rest, main = "Resting (control)")
plot2D(andreyev2009activ, main = "Kdo2-lipid A-treated")
addLegend(andreyev2009activ)
}

\keyword{datasets}
