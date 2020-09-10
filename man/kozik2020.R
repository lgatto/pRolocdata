\name{Kozik_con}
\alias{Kozik_con}

\docType{data}

\title{
  Small molecule enhancers of endosome-to-cytosol import augment anti-tumour immunity
}

\description{
  
  Data from 'Small molecule enhancers of endosome-to-cytosol import augment anti-tumour immunity' 
  
  Efficient cross-presentation of antigens by dendritic cells (DCs) is critical for initiation
  of anti-tumour immune responses. Yet, several steps of antigen intracellular traffic during cross-presentation
  are incompletely understood: in particular, the molecular mechanisms and the relative importance of antigen import
  from endocytic compartments into the cytosol. Here, we asked whether antigen import into the cytosol is rate-limiting
  for cross-presentation and anti-tumour immunity. By screening 700 FDA-approved drugs, we identified 37 import
  enhancers. We focused on prazosin and tamoxifen, and generated proteomic organellar maps of drug-treated DCs,
  covering the subcellular localisations of over 2000 proteins. By combining organellar mapping, quantitative
  proteomics, microscopy, and bioinformatics, we conclude that import enhancers undergo lysosomal
  trapping leading to membrane permeation and antigen release into the cytosol. Enhancing antigen
  import facilitates cross-presentation of both soluble and cell-associated antigens. Systemic
  administration of prazosin also led to reduced growth of MC38 tumours and to a synergistic
  effect with checkpoint immunotherapy in a melanoma model. Thus, inefficient antigen import
  into the cytosol limits antigen cross-presentation, restraining the potency of anti-tumour
  immune responses and efficacy of checkpoint blockers.
  
}

\usage{
  data("Kozik_con")
  data("Kozik_pra")
  data("Kozik_tam")
  data("Kozik2020")
}

\format{
  The data is an instance of class \code{MSnSet} from package
  \code{MSnbase}.
}

\examples{
  data(Kozik_con)
  Kozik_con
  pData(Kozik_con)
  exprs(Kozik_con)[1:3,1:3]
  
  library("pRoloc")
  plot2D(Kozik_con, main = "denderitic cells control")
}

\keyword{datasets}