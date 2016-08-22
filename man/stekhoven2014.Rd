\name{stekhoven2014}
\alias{stekhoven2014}

\docType{data}

\title{
  Data from Stekhoven et al. 2014
}

\description{

  Proteomics data provide unique insights into biological systems,
  including the predominant subcellular localization (SCL) of proteins,
  which can reveal important clues about their functions. Here we
  analyzed data of a complete prokaryotic proteome expressed under two
  conditions mimicking interaction of the emerging pathogen Bartonella
  henselae with its mammalian host. Normalized spectral count data from
  cytoplasmic, total membrane, inner and outer membrane fractions
  allowed us to identify the predominant SCL for 82% of the identified
  proteins. The spectral count proportion of total membrane versus
  cytoplasmic fractions indicated the propensity of cytoplasmic proteins
  to co-fractionate with the inner membrane, and enabled us to
  distinguish cytoplasmic, peripheral inner membrane and bona fide inner
  membrane proteins. Principal component analysis and k-nearest neighbor
  classification training on selected marker proteins or predominantly
  localized proteins, allowed us to determine an extensive catalog of at
  least 74 expressed outer membrane proteins, and to extend the SCL
  assignment to 94\% of the identified proteins, including 18% where in
  silico methods gave no prediction. Suitable experimental proteomics
  data combined with straightforward computational approaches can thus
  identify the predominant SCL on a proteome-wide scale. Finally, we
  present a conceptual approach to identify proteins potentially
  changing their SCL in a condition-dependent fashion.

}

\usage{data("stekhoven2014")}


\references{

  Stekhoven DJ, Omasits U, Quebatte M, Dehio C, Ahrens
  CH. \emph{Proteome-wide identification of predominant subcellular
  protein localizations in a bacterial model organism.} J
  Proteomics. 2014 Mar
  17;99:123-37. doi:10.1016/j.jprot.2014.01.015. Epub 2014 Jan
  28. PubMed PMID: 24486812.

}

\examples{
data(stekhoven2014)
library("pRoloc")
plot2D(stekhoven2014)
}

\keyword{datasets}
