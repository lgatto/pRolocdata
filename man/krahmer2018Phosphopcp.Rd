\name{krahmer2018Phosphopcp}
\alias{krahmer2018Phosphopcp}
\docType{data}

\title{
  Subcellular Reorganization in Diet-Induced Hepatic Steatosis
}

\description{
  
  Data from 'Organellar Proteomics and Phospho-Proteomics Reveal Subcellular Reorganization in Diet-Induced Hepatic Steatosis' 
  Developmental cell.
  
  Lipid metabolism is highly compartmentalized between cellular organelles that
  dynamically adapt their compositions and interactions in response to metabolic 
  challenges. Here, we investigate how diet-induced hepatic lipid accumulation,
  observed in non-alcoholic fatty liver disease (NAFLD), affects protein localization,
  organelle organization, and protein phosphorylation in vivo. We develop a mass 
  spectrometric workflow for protein and phosphopeptide correlation profiling to 
  monitor levels and cellular distributions of 6,000 liver proteins and 16,000 
  phosphopeptides during development of steatosis. Several organelle contact 
  site proteins are targeted to lipid droplets (LDs) in steatotic liver,
  tethering organelles orchestrating lipid metabolism. Proteins of the
  secretory pathway dramatically redistribute, including the mis-localization
  of the COPI complex and sequestration of the Golgi apparatus at LDs. This 
  correlates with reduced hepatic protein secretion. Our systematic in vivo
  analysis of subcellular rearrangements and organelle-specific phosphorylation
  reveals how nutrient overload leads to organellar reorganization and cellular dysfunction.
  
}

\usage{
  data("krahmer2018Phosphopcp")
}

\format{
  The data is an instance of class \code{MSnSet} from package
  \code{MSnbase}.
}

\references{
  
  \emph{Organellar proteomics and phospho-proteomics reveal subcellular
    reorganization in diet-induced hepatic steatosis.} Krahmer N, Najafi B,
  Schueder F, Quagliarini F, Steger M, Seitz S, Kasper R, Salinas F, Cox J,
  Uhlenhaut NH, Walther TC.  Developmental Cell. 2018 Oct 22;47(2):205-21.
  https://doi.org/10.1016/j.devcel.2018.09.017
  
}

\examples{
  data(krahmer2018Phosphopcp)
  krahmer2018Phosphopcp
  pData(krahmer2018Phosphopcp)
  exprs(krahmer2018Phosphopcp)[1:3,1:3]
  
  library("pRoloc")
  plot2D(krahmer2018Phosphopcp, main = "Kramer 2018")
}

\keyword{datasets}