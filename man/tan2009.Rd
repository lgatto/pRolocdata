\name{tan2009}
\alias{tan2009}
\alias{tan2009r1}
\alias{tan2009r2}
\alias{tan2009r3}

\docType{data}

\title{
  LOPIT data from Tan et al. 2009
}

\description{
  This is the data from Tan et al., \emph{Mapping organelle proteins and
    protein complexes in Drosophila melanogaster}, J Proteome Res. 2009
  Jun;8(6):2667-78. See below for more details. 
}

\usage{
data(tan2009r1)
data(tan2009r2)
data(tan2009r3)
}

\format{
  The data is an instance of class \code{MSnSet} from package \code{MSnbase}.
}

\details{
  This is a LOPIT experiment. Normalised intensities for proteins
  for four iTRAQ 4-plex labelled fractions are available for 3
  replicates (r1, r2 and r3 respectively). The markers used in the
  papers are available as \code{markers} feature meta-data and the
  partial least square discriminant analysis results as \code{PLSDA}
  feature meta-data.
}

\source{
  Supporting Information on
  \url{http://pubs.acs.org/doi/full/10.1021/pr800866n}
}

\references{
  Mapping organelle proteins and protein complexes in Drosophila 
  melanogaster. Tan DJ, Dvinge H, Christoforou A, Bertone P, Martinez 
  Arias A, Lilley KS. J Proteome Res. 2009 Jun;8(6):2667-78.
  PMID: 19317464
}

\examples{
data(tan2009r1)
tan2009r1
pData(tan2009r1)
head(exprs(tan2009r1))
table(fData(tan2009r1)$markers)
table(fData(tan2009r1)$PLSDA)
}

\keyword{datasets}
