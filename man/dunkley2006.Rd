\name{dunkley2006}
\alias{dunkley2006}
\docType{data}
\title{
  LOPIT data from Dunkley et al. 2006
}

\description{
  This is the data from Dunkley et al., \emph{Mapping the Arabidopsis
    organelle proteome}, PNAS 2006 (PMID 16618929). See below for more
  details. 
}

\usage{data(dunkley2006)}

\format{
  The data is an instance of class \code{MSnSet} from package \code{MSnbase}.
}

\details{
  This is a LOPIT experiment. Normalised intensities for 689 proteins
  for four iTRAQ 4-plex labelled fraction and 2 membrane preparation in
  duplicate (16 samples, see \code{phenoData(dunkley2006)} for more
  details) are provided.
  Partial least square discriminant analysis (PLSDA) has originally been
  applied to the test data \code{fData(dunkley)$markers)}; assignment results are
  available with \code{fData(dunkley)$assigned)} for 5 organelles.
  This dataset was also used in testing the phenotype discovery algorithm from 
  Breckels et al., \emph{The Effect of Organelle Discovery upon Sub-Cellular 
  Protein Localisation}, J Proteomics, \emph{In Press.}, see \code{phenoDisco}. 
  New phenotype clusters identified from algorithm application are available as 
  \code{pd.2013} feature meta-data.
}

\source{
  Supporting Information on
  \url{http://www.pnas.org/content/103/17/6518.abstract}. 
}

\references{
  Dunkley TP, Hester S, Shadforth IP, Runions J, Weimar T, Hanton SL, Griffin
  JL, Bessant C, Brandizzi F, Hawes C, Watson RB, Dupree P, Lilley
  KS. \emph{Mapping the Arabidopsis organelle proteome.} Proc Natl Acad
  Sci U S A. 2006 Apr 25;103(17):6518-23. Epub 2006 Apr 17. PubMed PMID:
  16618929; PubMed Central PMCID: PMC1458916.

  Breckels LM, Gatto L, Christoforou A, Groen AJ, Lilley KS and Trotter MWB. 
  \emph{The Effect of Organelle Discovery upon Sub-Cellular Protein Localisation}
   J Proteomics. \emph{In Press.}
}

\examples{
data(dunkley2006)
dunkley2006 
phenoData(dunkley2006)
## Input training data (organelle markers)
table(fData(dunkley2006)$markers)
## PLSDA assignment results
table(fData(dunkley2006)$assigned)
## PhenoDisco results
table(fData(dunkley2006)$pd.2013)
}

\keyword{datasets}
