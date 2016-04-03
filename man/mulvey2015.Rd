\name{mulvey2015}
\alias{mulvey2015}
\alias{mulvey2015norm}

\docType{data}

\title{
  Data from Mulvey et al. 2015
}

\description{
  This is the data from Mulvey et al., \emph{Dynamic proteomic profiling
  of extra-embryonic endoderm differentiation in mouse embryonic stem
  cells. }, Stem Cell. (PMID 26059426). See below for more details.
}

\usage{
data(mulvey2015)
data(mulvey2015norm)
}

\format{
  The data is an instance of class \code{MSnSet} from package \code{MSnbase}.
}

\details{

  While not a spatial proteomics data, it was analysed with the
  \code{pRoloc} package.

  During mammalian preimplantation development, the cells of the
  blastocyst's inner cell mass differentiate into the epiblast and
  primitive endoderm lineages, which give rise to the fetus and
  extra-embryonic tissues, respectively. Extra-embryonic endoderm (XEN)
  differentiation can be modeled in vitro by induced expression of GATA
  transcription factors in mouse embryonic stem cells. Here, we use this
  GATA-inducible system to quantitatively monitor the dynamics of global
  proteomic changes during the early stages of this differentiation
  event and also investigate the fully differentiated phenotype, as
  represented by embryo-derived XEN cells. Using mass spectrometry-based
  quantitative proteomic profiling with multivariate data analysis
  tools, we reproducibly quantified 2,336 proteins across three
  biological replicates and have identified clusters of proteins
  characterized by distinct, dynamic temporal abundance profiles. We
  first used this approach to highlight novel marker candidates of the
  pluripotent state and XEN differentiation. Through functional
  annotation enrichment analysis, we have shown that the downregulation
  of chromatin-modifying enzymes, the reorganization of membrane
  trafficking machinery, and the breakdown of cell-cell adhesion are
  successive steps of the extra-embryonic differentiation process. Thus,
  applying a range of sophisticated clustering approaches to a
  time-resolved proteomic dataset has allowed the elucidation of complex
  biological processes which characterize stem cell differentiation and
  could establish a general paradigm for the investigation of these
  processes.

}

\source{
  Supporting Information on

}

\references{

  Mulvey CM, Schr\"oter C, Gatto L, Dikicioglu D, Fidaner IB,
  Christoforou A, Deery MJ, Cho LT, Niakan KK, Martinez-Arias A, Lilley
  KS. Dynamic Proteomic Profiling of Extra-Embryonic Endoderm
  Differentiation in Mouse Embryonic Stem Cells. Stem Cells. 2015
  Sep;33(9):2712-25. doi: 10.1002/stem.2067. Epub 2015 Jun 23. PubMed
  PMID: 26059426.  }

\examples{
data(mulvey2015)
plot2D(mulvey2015)

data(mulvey2015norm)
heatmap(exprs(mulvey2015))
}

\keyword{datasets}
