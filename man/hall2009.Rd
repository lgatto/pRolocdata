\name{hall2009}
\alias{hall2009}
\docType{data}

\title{
  LOPIT data from Hall et al. 2009
}

\description{
  This is the data from Hall et al. \emph{The Organelle Proteome of the
  DT40 Lymphocyte Cell Line} Mol Cell Proteomics. 2009
  Jun;8(6):1295-305. (PMID: PMC2690488).
} 


\usage{data(hall2009)}


\format{
  An instance of class \code{MSnSet} from package \code{MSbase}.
}

\details{
  See reference.
}

\source{
  http://www.mcponline.org/content/8/6/1295.abstract
}

\references{
Hall SL, Hester S, Griffin JL, Lilley KS, Jackson AP. \emph{The
  organelle proteome of the DT40 lymphocyte cell line} Mol Cell
Proteomics. 2009 Jun;8(6):1295-305. doi:
10.1074/mcp.M800394-MCP200. Epub 2009 Jan 30. PubMed PMID: 19181659;
PubMed Central PMCID: PMC2690488.
}

\examples{
data(hall2009)
pData(hall2009)
library("pRoloc")
plot2D(hall2009)
}

\keyword{datasets}
