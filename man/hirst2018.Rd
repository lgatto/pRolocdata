\name{hirst2018}
\alias{hirst2018}


\docType{data}

\title{
  Data from Hirst et al. 2018
}

\description{

  From the supplementary file notes:

  These are the SILAC ratio data from 2046 proteins with complete
  profiles across all nine organellar maps.

  Each profile consists of five ratios, corresponding to five
  subcellular fractions obtained by differential centrifugation (3000 x
  g pellet, 6000 x g pellet, 12000 x g pellet, 24000 x pellet, 80000 x g
  pellet). The centrifugation speeds are available \code{centriguation}
  in the \code{MSnSet} object.

  Each ratio shows the abundance of the total membrane SILAC heavy
  spike-in relative to the abundance in a given subfraction.

  Maps were made from three cell lines (control HeLa, and two
  independent AP5Z1 KO HeLa cell lines, called AP5KNC2 and AP5KOC6),
  each in triplicate (replicates R1, R2, and R3). The sample are code as
  \code{"CTRL"} (HeLa control), \code{"C2"} (AP5KNC2 AP5Z1 KO cells) and
  \code{"C6"} (AP5KNC6 AP5Z1 KO cells).

  Marker proteins used to define organellar clusters in Supplemental
  Figure 1 in the manuscript are annotated as feature variable
  \code{markers}.

  Finally, the ratios in the \code{hirst2018} data where normalised by
  their sum (using \code{normalise(, method = "sum")}).

  The feature data also contains information about the comparison of
  organellar maps made from control or AP5 ablated cells, revealing
  putative proteins that undergo subcellular localisation shifts. Each
  protein receives an M score (magnitude of movement), and an R score
  (reproducibility of movement, i.e. correlation between replicates). In
  addition, the reproducibility of movement between the two AP5 KO cell
  lines is scored (Correlation C2 vs C6). Note however that the authors
  themselves claim that:

  \sQuote{The cutoffs chosen in Fig 1C (M > 1.5, R > 0.5) correspond to
    an estimated FDR of 23\%. Please note that the actual FDR is
    probably lower than this estimated FDR, because the mock data lack
    the additional cell line and the clonal correlation filter.}

  The re-localisation candidates are those that have an M score > 1.5
  and a R score > 0.5, and are marked with a hit feature variable set to
  \code{TRUE}.

}

\usage{
data(hirst2018)
}

\format{
  The data is an instance of class \code{MSnSet} from package
  \code{MSnbase}.
}


\references{
  Hirst J, Itzhak DN, Antrobus R, Borner GHH, Robinson MS. Role of the
  AP-5 adaptor protein complex in late endosome-to-Golgi retrieval. PLoS
  Biol. 2018 Jan 30;16(1):e2004411. doi:
  10.1371/journal.pbio.2004411. eCollection 2018 Jan. PubMed PMID:
  29381698; PubMed Central PMCID: PMC5806898.
}

\examples{
## load the two 24 hours datasets
data(hirst2018)
hirst2018

## experimental design
table(pData(hirst2018)[, -2])

## the expression data
exprs(hirst2018)[1:5, 1:3]

## abstract
abstract(hirst2018)

## split data by samples
x <- split(hirst2018, "sample")

## These are the relocalisation hits
hits <- which(fData(hirst2018)$Hits)
reloc <- FeaturesOfInterest(description = "Relocation hits",
			    featureNames(hirst2018)[hits])
reloc

## plotting
library("pRoloc")
par(mfrow = c(1, 3))
plot2D(x[[1]], main = "AP5KNC2")
highlightOnPlot(x[[1]], reloc)
plot2D(x[[2]], main = "AP5KNC6")
highlightOnPlot(x[[1]], reloc)
plot2D(x[[3]], main = "HeLa control")
highlightOnPlot(x[[1]], reloc)
addLegend(x[[3]], where = "topleft")
}

\keyword{datasets}
