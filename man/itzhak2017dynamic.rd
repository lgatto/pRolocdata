\name{itzhak2016dynamic}
\alias{helaCtrl}
\alias{helaEgf}
\docType{data}

\title{
  Global, quantitative and dynamic mapping of protein subcellular localization
}

\description{
  
  Data from 'Global, quantitative and dynamic mapping of protein subcellular localization.
  
  Subcellular localization critically influences protein function, and cells control
  protein localization to regulate biological processes. We have developed and applied
  Dynamic Organellar Maps, a proteomic method that allows global mapping of protein
  translocation events. We initially used maps statically to generate a database with
  localization and absolute copy number information for over 8700 proteins from HeLa cells,
  approaching comprehensive coverage. All major organelles were resolved, with exceptional
  prediction accuracy (estimated at >92\%). Combining spatial and abundance information
  yielded an unprecedented quantitative view of HeLa cell anatomy and organellar composition,
  at the protein level. We subsequently demonstrated the dynamic capabilities of the approach
  by capturing translocation events following EGF stimulation, which we integrated into a
  quantitative model. Dynamic Organellar Maps enable the proteome-wide analysis of
  physiological protein movements, without requiring any reagents specific to the
  investigated process, and will thus be widely applicable in cell biology.
  
}

\usage{
  data("itzhak2016helaCtrl")
  data("itzhak2016helaEgf")
}

\format{
  The data is an instance of class \code{MSnSet} from package
  \code{MSnbase}.
}

\references{
  
  \emph{Itzhak DN, Tyanova S, Cox J, Borner GH. Global, quantitative and dynamic mapping
  of protein subcellular localization. Elife. 2016 Jun 9;5:e16950.}
  
}

\examples{
  data("itzhak2016helaCtrl")
  helaCtrl <- itzhak2016helaCtrl
  pData(helaCtrl)
  exprs(helaCtrl)[1:3,1:3]
  
  library("pRoloc")
  plot2D(helaCtrl, main = "HeLa Ctrl", dims = c(1, 3))
}

\keyword{datasets}