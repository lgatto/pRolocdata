\name{ToxoLopit}
\alias{Barylyuk2020ToxoLopit}

\docType{data}

\title{
  Whole-cell spatial proteome of Toxoplasma: molecular anatomy of an apicomplexan cell
}

\description{
  
  Data from 'Whole-cell spatial proteome of Toxoplasma: molecular anatomy 
  of an apicomplexan cell' 
  
  Apicomplexan parasites are the causative agents of major human diseases
  and food insecurity and owe their considerable success to novel, highly
  specialized cell compartments and structures. These adaptations facilitate
  the recognition and non-destructive penetration of their host cells, and
  elaborate reengineering of these cells to promote growth, dissemination and
  active countering of host defense responses. The evolution of apicomplexan
  compartments is concomitant with great proteomic novelty that defines these
  new cell organizations and functions and, hence, approximately half of their
  proteins are unique and uncharacterized. Consequently, apicomplexan cells are
  relatively poorly understood. Here we employ the hyperLOPIT cell spatial
  proteomic method to the apicomplexan Toxoplasma gondii and define the
  steady-state subcellular location of thousands of proteins simultaneously
  giving comprehensive definition to these cells and their compartments. 
  These data, moreover, provide new insight into the spatial organizations
  of protein expression, adaptation to hosts, and the underlying evolutionary
  trajectories of these parasites. 
  
}

\usage{
  data("Barylyuk2020ToxoLopit")
}

\format{
  The data is an instance of class \code{MSnSet} from package
  \code{MSnbase}.
}

\examples{
  data(Barylyuk2020ToxoLopit)
  Barylyuk2020ToxoLopit
  pData(Barylyuk2020ToxoLopit)
  exprs(Barylyuk2020ToxoLopit)[1:3,1:3]
  
  library("pRoloc")
  plot2D(Barylyuk2020ToxoLopit, main = "Davies 2018 HeLa - wt")
}

\keyword{datasets}