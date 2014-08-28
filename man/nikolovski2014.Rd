\name{nikolovski2014}
\alias{nikolovski2014}

\docType{data}

\title{
  Data from Nikolovski et al. 2014
}


\description{
  This is the data used in Nikolovksi et al. (2014). See below for 
  details and references.  
}

\usage{
data(nikolovski2014)
}

\format{
  The data is an instance of class \code{MSnSet} from package \code{MSnbase}.
}

\details{

  Abstract: The proteomic composition of the Arabidopsis Golgi apparatus
  is currently reasonably well documented; however little is known about
  the relative abundances between different proteins within this
  compartment. Accurate quantitative information of Golgi resident
  proteins is of great importance: it facilitates a better understanding
  of the biochemical processes which take place within this organelle,
  especially those of different polysaccharide synthesis pathways. Golgi
  resident proteins are challenging to quantify since the abundance of
  this organelle is relatively low within the cell. In this study an
  organelle fractionation approach, targeting the Golgi apparatus, was
  combined with a label free quantitative mass spectrometry (MS),
  data-independent acquisition (DIA) method employing ion mobility
  separation known as LC-IMS-MSE (or HDMSE), to simultaneously localize
  proteins to the Golgi apparatus and assess their relative quantity. In
  total 102 Golgi localised proteins were quantified. These data provide
  new insight into Golgi apparatus organization and demonstrate that
  organelle fractionation in conjunction with label free quantitative MS
  is a powerful and relatively simple tool to access protein organelle
  localisation and their relative abundances. The findings presented
  open a unique view on the organization of the plant Golgi apparatus,
  leading towards novel hypotheses centered on the biochemical processes
  of this organelle.
  


  These data are a concatenation of 2 LOPIMS gradients, labelled
  gradient A and B, each with 10 fractions.

}

\source{
  Supplemental Data downloaded from 
  \url{http://www.plantphysiol.org/content/early/2014/08/13/pp.114.245589/suppl/DC1},
  also available in the package's \code{extdata} directory.
}

\references{
  Nikolovski N, Shliaha PV, Gatto L, Dupree P, Lilley KS. Label free
  protein quantification for plant Golgi protein localisation and
  abundance. Plant Physiol.  2014 Aug 13. pii: pp.114.245589. [Epub
  ahead of print] PubMed PMID: 25122472.
}

\examples{
data(nikolovski2014)
pData(nikolovski2014)
library("pRoloc")
plot2D(nikolovski2014)
addLegend(nikolovski2014, where = "topright", bty = "n", cex =.7)

A <- pData(nikolovski2014)$gradient == "A"
par(mfrow = c(1, 2))
plot2D(nikolovski2014[, A], main = "Gradient A")
plot2D(nikolovski2014[, !A], main = "Gradient B")
}

\keyword{datasets}
