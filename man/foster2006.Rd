\name{foster2006}
\alias{foster2006}
\docType{data}

\title{
  PCP data from Foster et al. (2006)
}

\description{
  This is the data from Foster et al., \emph{A Mammalian Organelle Map
    by Protein Correlation Profiling}, Cell 2006 (PMID 16615899). See
  below for more details.
}

\usage{data(foster2006)}
\format{
  The data is an instance of class \code{MSnSet} from package \code{MSnbase}.
}

\details{
  This is a PCP experiment. Label-free quantification has been done on
  a totla of 26 high and low density-separated fractions (see
  \code{pData(foster2006)}). A total of 1555 proteins have been
  quantified in a subset of the fractions. The proteins are described in
  the \code{featureData} slot. Chi^2 calculations, as defined in the PCP
  experiment, has been performed using marker proteins for a total of 8
  organelles, as well as the authors' original assignment and notes are
  available in the \code{featureData} slot.
}

\source{
  Supplemental data on
  \url{http://www.cell.com/abstract/S0092-8674(06)00369-2}.
}

\references{
  Foster LJ, de Hoog CL, Zhang Y, Zhang Y, Xie X, Mootha VK, Mann
  M. \emph{A mammalian organelle map by protein correlation profiling.}
  Cell. 2006 Apr 7;125(1):187-99. PubMed PMID: 16615899.
}

\examples{
data(foster2006)
foster2006
phenoData(foster2006)
featureData(foster2006)
## organelle marker proteins
table(fData(foster2006)$train)
}
\keyword{datasets}
