\name{nikolovski2012}
\alias{nikolovski2012}
\alias{nikolovski2012imp}

\docType{data}

\title{
  Data from Nikolovski et al. 2012
}


\description{
  This is the data from Dunkley et al., \emph{Mapping the Arabidopsis
    organelle proteome}, PNAS 2006 (PMID 16618929). See below for more
  details. 
}

\usage{
data(nikolovski2012)
data(nikolovski2012imp)
}

\format{
  The data is an instance of class \code{MSnSet} from package \code{MSnbase}.
}

\details{
  These data are a concatenation of 4 LOPIT experiments. Experiments 1
  and 2 are from Dunkley et al. 2006 (see also \code{dunkley2006}).
  Exepriments 3 and 4 are new.

  In the LOPIT experiments by Dunkley et al. (2006), peripheral membrane
  proteins were removed by carbonate washing of the isolated membranes,
  while for experiments 3 and 4, no carbonate wash was performed and
  are, as such, enriched in peripheral and luminal proteins.
  See figure 1 in Nikolovski 2012 for a description of the design.

  In \code{nikolovksi2012imp} missing values have been imputed using
  partial least-squares regression.

  The training set used for the Naive Bayesian classifier is available
  as the \code{markers} feature meta-data.
}

\source{
  Supporting Information on  
  \url{http://www.plantphysiol.org/content/160/2/1037.long},
  also available in the package's \code{extdata} directory.
}

\references{
Nikolovski N, Rubtsov D, Segura MP, Miles GP, Stevens TJ, Dunkley TP,
Munro S, Lilley KS, Dupree P. \emph{Putative glycosyltransferases and
  other plant Golgi apparatus proteins are revealed by LOPIT
  proteomics.} 
Plant Physiol. 2012 Oct;160(2):1037-51. doi: 10.1104/pp.112.204263. 
Epub 2012 Aug 24. PMID: 22923678; PMCID: PMC3461528.

Dunkley TP, Hester S, Shadforth IP, Runions J, Weimar T, Hanton SL, Griffin
JL, Bessant C, Brandizzi F, Hawes C, Watson RB, Dupree P, Lilley
KS. \emph{Mapping the Arabidopsis organelle proteome.} Proc Natl Acad
Sci U S A. 2006 Apr 25;103(17):6518-23. Epub 2006 Apr 17. PubMed PMID:
16618929; PubMed Central PMCID: PMC1458916.
}

\examples{
data(nikolovksi2012)
data(nikolovksi2012imp)
table(is.na(nikolovski2012))
table(is.na(nikolovski2012imp))
phenoData(nikolovski2012)
table(fData(nikolovski2012)$markers)
all.equal(sort(featureNames(nikolovski2012)),
          sort(featureNames(nikolovski2012imp)))
}

\keyword{datasets}
