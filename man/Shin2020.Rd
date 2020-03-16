\name{Shin2010}
\alias{Shin2019MitoControl}
\alias{Shin2019MitoGCC88}
\alias{Shin2019MitoGol97}

\docType{data}

\title{
  Spatial proteomics defines the content of trafficking vesicles captured by golgin tethers
}

\description{
  
  Data from 'Spatial proteomics defines the content of trafficking vesicles captured by golgin tethers'
  
  Intracellular traffic between compartments of the secretory and endocytic pathways is mediated by
  vesicle-based carriers. The precise and complete proteomes of carriers destined for many
  organelles are ill-defined because the vesicular intermediates are transient, low-abundance and
  difficult to purify. Here, we combine vesicle relocalisation with organelle proteomics and Bayesian
  analysis to define the content of different endosome-derived vesicles destined for the trans-Golgi
  network (TGN). The golgin coiled-coil proteins golgin-97, golgin-245 and GCC88, shown previously
  to capture endosome-derived vesicles at the TGN, were individually relocalised to mitochondria
  and the content of subsequently re-routed vesicles was determined by organelle proteomics. Our
  findings revealed 45 integral and 51 peripheral membrane proteins re-routed by golgin-97,
  evidence for a distinct class of vesicles shared by golgin-97 and GCC88, and various cargoes
  specific to individual golgins. These results illustrate a general strategy for analysing intracellular
  sub-proteomes by combining acute cellular re-wiring with high-resolution spatial proteomics.
  
}

\usage{
  data("Shin2019MitoControl")
  data("Shin2019MitoGCC88")
  data("Shin2019MitoGol97")
}

\format{
  The data is an instance of class \code{MSnSet} from package
  \code{MSnbase}.
}

\references{
  
  \emph{Spatial proteomics defines the content of trafficking vesicles captured by golgin
    tethers. John J.H. Shin, Oliver M. Crook, Alicia C. Borgeaud, Jérôme Cattin-Ortolá, Sew-Yeu Peak-
      Chew, Jessica Chadwick, Kathryn S. Lilley, Sean Munro.
  
}

\examples{
  data(Shin2019MitoControl)
  Shin2019MitoControl
  pData(Shin2019MitoControl)
  exprs(Shin2019MitoControl)[1:3,1:3]
  
  library("pRoloc")
  plot2D(Shin2019MitoControl, main = "Shin 2019 HeK 293T WT")
}

\keyword{datasets}