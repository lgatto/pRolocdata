\name{fabre2015r1}
\alias{fabre2015r1}
\alias{fabre2015r2}
\alias{fabre2015}

\docType{data}

\title{
  Data from Fabre et al. 2015
}


\description{

  Duplicated experimental data from Fabre et al. 2015, \emph{Deciphering
  preferential interactions within supramolecular protein complexes: the
  proteasome case}. Protein complexes were separated by glycerol density
  gradient centrifugation. Proteins have been quantified by label-free
  (iBAQ) mass spectrometry.  
}

\usage{
data("fabre2015r1")
data("fabre2015r2")
}

\references{
  Fabre B, Lambour T, Garrigues L, Amalric F, Vigneron N, Menneteau T,
  Stella A, Monsarrat B, Van den Eynde B, Burlet-Schiltz O,
  Bousquet-Dubouch MP. Deciphering preferential interactions within
  supramolecular protein complexes: the proteasome case. Mol Syst
  Biol. 2015 Jan 5;11(1):771. doi: 10.15252/msb.20145497. PubMed PMID:
  25561571.
}



\examples{
data(fabre2015r1)
experimentData(fabre2015r1)
library("pRoloc")
plot2D(fabre2015r1)
addLegend(fabre2015r1, where = "topright")
}

\keyword{datasets}
