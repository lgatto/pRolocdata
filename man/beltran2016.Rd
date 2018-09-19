\name{beltran2016}
\alias{beltran2016}

\alias{beltran2016HCMV24}
\alias{beltran2016HCMV48}
\alias{beltran2016HCMV72}
\alias{beltran2016HCMV96}
\alias{beltran2016HCMV120}
\alias{beltran2016MOCK24}
\alias{beltran2016MOCK48}
\alias{beltran2016MOCK72}
\alias{beltran2016MOCK96}
\alias{beltran2016MOCK120}


\docType{data}

\title{
  Data from Beltran et al. 2016
}

\description{

  The data contain the spatial proteomics data for 5 time points (24,
  48, 72, 96 and 120) and two condidiont (HMCV infection and MOCK),
  totalling 10 \code{MSnSet} object. Each contains expression data for
  1748 to 2220 proteins along 6 fractions quantified by TMT 6-plex.

}

\usage{
data(beltran2016HCMV24)
data(beltran2016HCMV48)
data(beltran2016HCMV72)
data(beltran2016HCMV96)
data(beltran2016HCMV120)
data(beltran2016MOCK24)
data(beltran2016MOCK48)
data(beltran2016MOCK72)
data(beltran2016MOCK96)
data(beltran2016MOCK120)
}

\format{
  The data is an instance of class \code{MSnSet} from package
  \code{MSnbase}.
}


\references{

  Jean Beltran PM, Mathias RA, Cristea IM. \emph{A Portrait of the Human
  Organelle Proteome In Space and Time during Cytomegalovirus
  Infection}. Cell Syst. 2016 Oct 26;3(4):361-373.e6. doi:
  10.1016/j.cels.2016.08.012. Epub 2016 Sep 15. PubMed PMID: 27641956;
  PubMed Central PMCID: PMC5083158.

}

\examples{
## load the two 24 hours datasets
data(beltran2016MOCK24)
data(beltran2016HCMV24)
beltran2016MOCK24
beltran2016HCMV24

## the expression data
head(exprs(beltran2016MOCK24))
head(exprs(beltran2016HCMV24))

## abstract
abstract(beltran2016HCMV24)

## plotting
library("pRoloc")
par(mfrow = c(1, 2))
plot2D(beltran2016HCMV24, main = "HCMV 24hpi")
plot2D(beltran2016MOCK24, main = "MOCK 24hpi")

## Combine the date as a list and keep only common features
ml <- MSnSetList(list(beltran2016HCMV24, beltran2016MOCK24))
ml <- commonFeatureNames(ml)

par(mfrow = c(1, 2))
plot2D(ml[[1]], main = "HCMV 24hpi")
plot2D(ml[[2]], main = "MOCK 24hpi")

}

\keyword{datasets}
