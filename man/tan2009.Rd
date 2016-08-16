\name{tan2009}
\alias{tan2009}
\alias{tan2009r1}
\alias{tan2009r2}
\alias{tan2009r3}
\alias{tan2009r1goCC}

\docType{data}

\title{
  LOPIT data from Tan et al. (2009)
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
data(tan2009r1goCC)
}

\format{
  The data is an instance of class \code{MSnSet} from package \code{MSnbase}.
}

\details{
  This is a LOPIT experiment. Normalised intensities for proteins
  for four iTRAQ 4-plex labelled fractions are available for 3
  replicates (r1, r2 and r3 respectively). The partial least square discriminant 
  analysis results from the paper are available as \code{PLSDA} feature 
  meta-data and the markers used in analysis are available as \code{markers} 
  feature meta-data (Note: the ER and Golgi organelle markers were combined in
  original PLSDA analysis). 
  
  Replicate 1 was also used in testing the phenotype discovery algorithm from 
  Breckels et al., \emph{The Effect of Organelle Discovery upon Sub-Cellular 
  Protein Localisation}, J Proteomics, \emph{In Press.}, see \code{phenoDisco}. 
  New phenotype clusters identified from algorithm application are available as 
  \code{pd.2013} feature meta-data.  
  
  The \code{tan2009r1goCC} instance contains binary assay data. Its
  columns represent GO CC terms that have been observed for the object's
  features. A 1 indicates that a GO term has been associated to a given
  feature (protein); a 0 means not such association was found in the GO
  ontology.
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

  Breckels LM, Gatto L, Christoforou A, Groen AJ, Lilley KS and Trotter MWB. 
  \emph{The Effect of Organelle Discovery upon Sub-Cellular Protein Localisation}
   J Proteomics. \emph{In Press.}
}

\examples{
data(tan2009r1)
tan2009r1
pData(tan2009r1)
head(exprs(tan2009r1))
# Organelle markers
table(fData(tan2009r1)$markers)
# PLSDA assignment results
table(fData(tan2009r1)$PLSDA)
}

\keyword{datasets}
