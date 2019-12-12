\name{orre2019}
\alias{orre2019a431}
\alias{orre2019h322}
\alias{orre2019hcc827}
\alias{orre2019hcc827gef}
\alias{orre2019hcc827rep1}
\alias{orre2019hcc827rep2}
\alias{orre2019hcc827rep3}
\alias{orre2019mcf7}
\alias{orre2019u251}
\docType{data}

\title{
  SubCellBarCode: Proteome-wide Mapping
  of Protein Localization and Relocalization
}

\description{
  
  Data from 'SubCellBarCode: Proteome-wide Mapping of Protein Localization and Relocalization' 
  Molecular cell.
  
  Subcellular localization is a main determinant of protein function; 
  however, a global view of cellular proteome organization remains relatively unexplored. We
  have developed a robust mass spectrometry-based
  analysis pipeline to generate a proteome-wide view
  of subcellular localization for proteins mapping to
  12,418 individual genes across five cell lines. Based
  on more than 83,000 unique classifications and correlation profiling, we investigate the effect of alternative
  splicing and protein domains on localization, complex member co-localization, cell-type-specific localization,
  as well as protein relocalization after growth
  factor inhibition. Our analysis provides information
  about the cellular architecture and complexity of the
  spatial organization of the proteome; we show that
  the majority of proteins have a single main subcellular
  location, that alternative splicing rarely affects
  subcellular location, and that cell types are best
  distinguished by expression of proteins exposed to
  the surrounding environment. The resource is freely
  accessible via www.subcellbarcode.org.
  
}

\usage{
  data("orre2019a431")
  data("orre2019h322")
  data("orre2019hcc827")
  data("orre2019hcc827gef")
  data("orre2019hcc827rep1")
  data("orre2019hcc827rep2")
  data("orre2019hcc827rep3")
  data("orre2019mcf7")
  data("orre2019u251")
}

\format{
  The data is an instance of class \code{MSnSet} from package
  \code{MSnbase}.
}

\references{
  
  \emph{SubCellBarCode: Proteome-wide Mapping of Protein Localization and Relocalization} Lukas Minus Orre, 
  Mattias Vesterlund, Yanbo Pan, Taner Arslan, Yafeng Zhu,  Alejandro Fernandez Woodbridge, Oliver Frings,
  Erik Fredlund, and Janne Lehtio https://doi.org/10.1016/j.molcel.2018.11.035
  
}

\examples{
  data(orre2019a431)
  orre2019a431
  pData(orre2019a431)
  exprs(orre2019a431)[1:3,1:3]
  
  library("pRoloc")
  plot2D(orre2019a431,, main = "Orre 2019 A431")
}

\keyword{datasets}