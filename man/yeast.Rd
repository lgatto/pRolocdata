\name{yeast2018}
\alias{yeast2018}
\docType{data}

\title{
  Saccharomyces cerevisiae spatial proteomics (2018)
}

\description{

 Data from 'The subcellular organisation of Saccharomyces cerevisiae' 
 (submitted).

 This dataset represents four biological replicate hyperLOPIT 
 experiments performed in Saccharomyces cerevisiae cultured to 
 early-mid exponential phase, in synthetic media with glucose as sole 
 carbon source (SD-His media (Breker et al 2013)). These were carried 
 out to produce a map of the spatial proteome of this organism under 
 no-perturbed conditions. The associated quantitation data from these 
 experiments were combined using the method described in reference. 
 This dataset contains quantitative information for 2,847 proteins 
 that were common across our four biological replicate experiments and 
 information regarding localisation for all of the proteins in the 
 combined experiment. Overall this dataset describes 936 proteins that 
 localise to one of 12 subcellular locations in S. cerevisiae under 
 our experimental conditions.

}

\usage{
data("yeast2018")
}


\examples{
data(yeast2018)

library("pRoloc")
par(mfrow = c(1, 2))
plot2D(yeast2018, main = "Markers")
addLegend(yeast2018, where = "bottomleft", cex = .7)
plot2D(yeast2018, fcol = "predicted.location", main = "Localisation")
}

\keyword{datasets}
