\name{courtland_control}
\alias{courtland_control}
\alias{courtland_mutant}

\docType{data}

\title{
  Genetic Disruption of WASHC4 Drives Endo-lysosomal Dysfunction and Cognitive-Movement Impairments in Mice and Humans
}

\description{
  
  Data from ' Genetic Disruption of WASHC4 Drives Endo-lysosomal Dysfunction and Cognitive-Movement Impairments in Mice and Humans' 
  
  Mutation of the WASH complex subunit, SWIP, is
  implicated in human intellectual disability, but the cellular
  etiology of this association is unknown. We identify the neuronal
  WASH complex proteome, revealing a network of endosomal proteins.
  To uncover how dysfunction of endosomal SWIP leads to disease, we
  generate a mouse model of the human
  WASHC4 c.3056C>G mutation.  Quantitative spatial
  proteomics analysis of SWIP P1019R mouse brain
  reveals that this mutation destabilizes the WASH complex and
  uncovers significant  perturbations in both endosomal and lysosomal
  pathways.  Cellular and histological analyses confirm that
  SWIP P1019R results in  endo-lysosomal disruption
  and uncover indicators of neurodegeneration. We find that
  SWIP P1019R not only impacts cognition, but also
  causes significant progressive motor deficits in mice.  Remarkably,
  a retrospective analysis of SWIP P1019R patients
  confirms motor deficits in humans. Combined, these findings support
  the model that WASH complex destabilization, resulting from
  SWIP P1019R, drives cognitive and motor
  impairments via endo-lysosomal dysfunction in the brain.
  
}

\usage{
  data("courtland_control")
  data("courtland_mutant")
}

\format{
  The data is an instance of class \code{MSnSet} from package
  \code{MSnbase}.
}

\examples{
  data(courtland_control)
  courtland_control
  pData(courtland_control)
  exprs(courtland_control)[1:3,1:3]
  
  library("pRoloc")
  plot2D(courtland_control, main = "mouse brain control")
}

\keyword{datasets}