\name{groen2014}
\alias{groen2014}
\alias{groen2014r1}
\alias{groen2014r2}
\alias{groen2014r3}
\alias{groen2014cmb}

\docType{data}

\title{
  Data from LOPIT experiments on Arabidopsis thaliana callus roots,
  taken from Groen et al (2014)
}

\description{ This is the data from Groen et al. \emph{Identification of
  Trans Golgi Network proteins in Arabidopsis thaliana root tissue}
  J. Proteome Res., Just Accepted. Publication Date (Web): December 17,
  2013 (Article) DOI: 10.1021/pr4008464 }


\usage{data(groen2014r1)
data(groen2014r2)
data(groen2014r3)
data(groen2014cmb)
}


\format{
  An instance of class \code{MSnSet} from package \code{MSbase}.
}

\details{ This is a LOPIT experiment. Normalised intensities for
  proteins for four iTRAQ 4-plex labelled fractions are available for 3
  replicates (\code{r1}, \code{r2} and \code{r3} respectively). The 3
  replicates have also been combined as described in Groen et al. and
  Trotter et al. (2010) to generate a fourth dataset (\code{cmb}), also
  shown in the example code below.  }

\source{
  http://pubs.acs.org/doi/abs/10.1021/pr4008464
}

\references{
  Groen AJ, Sancho-Andres G, Breckels LM, Gatto L, Aniento F,
  and Lilley KS. \emph{Identification of Trans Golgi Network proteins in
  Arabidopsis thaliana root tissue}. J. Proteome Res., Just
  Accepted. Publication Date (Web): December 17, 2013 (Article). DOI:
  10.1021/pr4008464

  Trotter MWB, Sadowski PG, Dunkley TPJ, Groen AJ and Lilley KS.
  \emph{Improved sub-cellular resolution via simultaneous analysis of
  organelle proteomics data across varied experimental conditions.}
  Proteomics 2010 10(23):4213-4219. PMID 21058340.
  
  Sadowski PG, Groen AJ, Dupree P and Lilley KS. \emph{Sub-cellular
  localization of membrane proteins.}  Proteomics 2008
  8(19):3991-4011. PMID 18780351.

  Dunkley TP, Hester S, Shadforth IP, Runions J, Weimar T, Hanton SL,
  Griffin JL, Bessant C, Brandizzi F, Hawes C, Watson RB, Dupree P,
  Lilley KS. \emph{Mapping the Arabidopsis organelle proteome.} Proc
  Natl Acad Sci U S A. 2006 Apr 25;103(17):6518-23. Epub 2006 Apr
  17. PubMed PMID: 16618929; PubMed Central PMCID: PMC1458916.

}


\examples{
data(groen2014r1)
data(groen2014r2)
data(groen2014r3)
data(groen2014cmb)

## The combine dataset can generated manually using 
cmb <- combine(groen2014r1, updateFvarLabels(groen2014r2))
cmb <- filterNA(cmb)
cmb <- combine(cmb, updateFvarLabels(groen2014r3))
cmb <- filterNA(cmb)
fData(cmb) <- fData(cmb)[, c(1,2,5)]
cmb

## or can simply be loaded directly
data(groen2014cmb)

## check datsets are the same
all.equal(cmb, groen2014cmb, check.attributes=FALSE)
}

\keyword{datasets}

