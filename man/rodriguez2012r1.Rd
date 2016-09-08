\name{rodriguez2012r1}

\alias{rodriguez-pineiro2012}
\alias{rodriguez2012}
\alias{rodriguez2012r1}
\alias{rodriguez2012r2}
\alias{rodriguez2012r3}

\docType{data}

\title{

  Spatial proteomics of human inducible goblet-like LS174T cells.
  
}

\description{

  Data from Rodriguez-Pineiro AM, van der Post S, Johansson ME, Thomsson
  KA, Nesvizhskii AI, Hansson GC. Proteomic study of the mucin granulae
  in an intestinal goblet cell model. J Proteome Res. 2012 Mar
  2;11(3):1879-90. doi: 10.1021/pr2010988.  Epub 2012 Feb 2. PubMed
  PMID:22248381; PubMed Central PMCID:PMC3292267.

}

\usage{
data("rodriguez2012r1")
data("rodriguez2012r2")
data("rodriguez2012r3")
}



\details{

  As no marker were provided with the data, we transferred markers from
  the \code{hyperLOPIT2015} (mouse) data using gene names to match
  between experiments. To validate our marker annotations, we compared
  the relative distributions of our markers (see figure below) to the
  PCA plot provided by the authors (Figure 3). Both show a similar
  separation of mitochondion/ER vs the rest along PC1 and
  ribosomes/lysosome vs rest along PC2. The data do not match exactly as
  the different marker protein are used.

}

\source{

  The supplementary file is \code{pr2010988_si_003.xls}. See
  \code{scripts/rodriguez-pineiro2012.R} for data preparation.
  
}


\examples{
data(rodriguez2012r1)
data(rodriguez2012r2)
data(rodriguez2012r3)

library("pRoloc")
par(mfrow = c(2, 2))
plot2D(rodriguez2012r1)
addLegend(rodriguez2012r1, cex = .7, where = "topleft")
plot2D(rodriguez2012r2)
plot2D(rodriguez2012r3)


## compare to figure 3
dev.new()
plot2D(markerMSnSet(rodriguez2012r1),
       mirrorX = TRUE, mirrorY = TRUE,
       main = "Our markers")
addLegend(markerMSnSet(rodriguez2012r1), where = "bottomright")
}

\keyword{datasets}
