\name{at_chloro}
\alias{at_chloro}
\docType{data}

\title{
  The AT_CHLORO data base
}

\description{
  AT_CHLORO is a comprehensive chloroplast proteome database with
  subplastidial localization and curated information on envelope
  proteins.

  The \code{assayData} contains the raw spectral counts for 3
  chloroplastic fractions (the envelope, the stroma and the thylakoids)
  and for a complete chloroplast sample. The percentage of occurrence in
  each of the sub-chloroplast fraction as calculated in Ferro et
  al. (2010) are available as feature meta data (\code{Percent_ENV},
  \code{Percent_STR} and \code{Percent_THY}). 
}

\usage{data(at_chloro)}

\format{
  The data is an instance of class \code{MSnSet} from package \code{MSnbase}.
}

\source{
  Myriam Ferro
  Exploring the Dynamics of Proteomes (EDyP)
  Laboratoire Biologie a Grande Echelle (BGE)
  U1038 INSERM/CEA/UJF
  Institut de Recherches en Technologies et Sciences pour le Vivant (iRTSV)
  CEA/Grenoble
}

\references{
  Ferro M, Brugiere S, Salvi D, Seigneurin-Berny D, Court M, Moyet L,
  Ramus C, Miras S, Mellal M, Le Gall S, Kieffer-Jaquinod S, Bruley C,
  Garin J, Joyard J, Masselon C, Rolland  N. AT_CHLORO, a comprehensive
  chloroplast proteome database with subplastidial localization and
  curated information on envelope proteins. Mol Cell Proteomics. 2010
  Jun;9(6):1063-84. Epub 2010 Jan 10. PubMed PMID: 20061580; PubMed
  Central PMCID: PMC2877971 
}

\examples{
data(at_chloro)
dim(at_chloro)
pData(at_chloro)
head(exprs(at_chloro))
fvarLabels(at_chloro)
table(fData(at_chloro)$markers)
## check exprs data and 'TotalSpectralCount' feature meta data
all(fData(at_chloro)$TotalSpectralCount == rowSums(exprs(at_chloro)))
}

\keyword{datasets}
