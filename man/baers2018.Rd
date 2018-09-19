\name{baers2018}
\alias{baers2018}
\alias{synechocystis}
\docType{data}

\title{
  Synechocystis spatial proteomics
}

\description{

 Data from 'Spatial mapping of a cyanobacterial proteome reveals
 distinct subcellular compartment organisation and dynamic metabolic
 pathways' (submitted).

 Cyanobacteria are complex prokaryotes, incorporating a Gram-negative
 cell wall and internal thylakoid membranes. However, localisation of
 proteins within cyanobacterial cells is poorly understood. Using
 subcellular fractionation and quantitative proteomics we report the
 most extensive subcellular map of the proteome of a cyanobacterial
 cell, identifying ~67\% of Synechocystis sp. PCC 6803 proteins, ~1000
 more than previous studies. 1711 proteins were assigned to six specific
 subcellular regions.

 This dataset is composed of two combined replicated 10-plex LOPIT
 experiments.

 Protein markers for the plasma membrane, thylakoid membrane, cytosol,
 and small and la rge ribosomal subunits were curated from a literature
 review. A Support Vector Machine (SVM) classifier was employed on the
 combined dataset, with a radial basis function kernel, using class
 specific weights for classification of unassigned proteins to one of
 the 5 defined sub-cellular niches, TM, PM, soluble, small ribosomal
 subunit, large ribosomal subunit. The weights used in classification
 were set to be inversely proportional to the subcellular class
 frequencies to account for class imbalance. Algorithmic performance of
 the SVM on the dataset was estimated (as described in Trotter et al
 2010). Scoring thresholds were calculated per subcellular niche and
 were set based on concordance with existing subcellular knowledge
 annotation to attain a 7.5\% false discovery rate (FDR). Unassigned
 proteins were then classified to 1 of the 5 compartments according to
 the SVM prediction if greater than the calcul ated class threshold.

}

\usage{
data("baers2018")
}


\examples{
data(baers2018)

library("pRoloc")
par(mfrow = c(1, 2))
plot2D(baers2018, main = "Markers")
addLegend(baers2018, where = "bottomright")
plot2D(baers2018, fcol = "Localisation", main = "Localisation")
}

\keyword{datasets}
