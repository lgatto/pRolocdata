\name{itzhak2016stcSILAC}
\alias{itzhak2016}
\alias{itzhak2016stcSILAC}

\docType{data}

\title{
  Data from Itzhak et al. 2016
}

\description{

  Data from Daniel N Itzhak, Stefka Tyanova, Jurgen Cox and Georg HH
  Borner. Global, quantitative and dynamic mapping of protein
  subcellular
  localization. DOI:\url{http://dx.doi.org/10.7554/eLife.16950}
  Published June 9, 2016 Cite as eLife 2016;10.7554/eLife.16950

  It currently contains

  \itemize{

    \item The second sheet contains the 6 replicates of the SILAC static
    data (\emph{Static} data were used to genrate six deep organellar
    maps) and is made available as \code{itzhak2016stcSILAC}.

  }
}

\usage{data("itzhak2016stcSILAC")}

\source{

  This data was generated from Supplementary file 9
  (\url{https://elifesciences.org/content/5/e16950/supp-material9}). See
  \code{inst/scripts/itzhak2016.R} for the exact code
  
}


\examples{
data(itzhak2016stcSILAC)
itzhak2016stcSILAC
dim(itzhak2016stcSILAC)
pData(itzhak2016stcSILAC)

## only 1st replicate
dim(itzhak2016stcSILAC[, itzhak2016stcSILAC$rep == 1])


## filter out features with missing values
itzhak2016stcSILAC <- filterNA(itzhak2016stcSILAC)

library("pRoloc")
## Cell map
plot2D(itzhak2016stcSILAC)
## as in the paper
plot2D(itzhak2016stcSILAC, dims = c(1, 3))
}

\keyword{datasets}
