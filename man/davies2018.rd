\name{davies2018}
\alias{davies2018ap4b1}
\alias{davies2018ap4e1}
\alias{davies2018wt}

\docType{data}

\title{
  AP-4 vesicles contribute to spatial control of autophagy via RUSC-dependent peripheral delivery of ATG9A
}

\description{

  Data from 'AP-4 vesicles contribute to spatial control of autophagy
  via RUSC-dependent peripheral delivery of ATG9A' Nature
  Communications.

  Adaptor protein 4 (AP-4) is an ancient membrane trafficking complex,
  whose function has largely remained elusive. In humans, AP-4
  deficiency causes a severe neurological disorder of unknown aetiology.
  We apply unbiased proteomic methods, including 'Dynamic Organellar
  Maps', to find proteins whose subcellular localisation depends on
  AP-4.  We identify three transmembrane cargo proteins, ATG9A, SERINC1
  and SERINC3, and two AP-4 accessory proteins, RUSC1 and RUSC2. We
  demonstrate that AP-4 deficiency causes missorting of ATG9A in diverse
  cell types, including patient-derived cells, as well as dysregulation
  of autophagy. RUSC2 facilitates the transport of AP-4-derived,
  ATG9A-positive vesicles from the trans-Golgi network to the cell
  periphery.  These vesicles cluster in close association with
  autophagosomes, suggesting they are the 'ATG9A reservoir' required for
  autophagosome biogenesis.  Our study uncovers ATG9A trafficking as a
  ubiquitous function of the AP-4 pathway.  Furthermore, it provides a
  potential molecular pathomechanism of AP-4 deficiency, through
  dysregulated spatial control of autophagy.

}

\usage{
  data("davies2018ap4b1")
  data("davies2018ap4e1")
  data("davies2018wt")
}

\format{
  The data is an instance of class \code{MSnSet} from package
  \code{MSnbase}.
}

\references{

  \emph{AP-4 vesicles contribute to spatial control of autophagy via
  RUSC-dependent peripheral delivery of ATG9A} Alexandra K. Davies,
  Daniel N. Itzhak, James R. Edgar, Tara L. Archuleta, Jennifer Hirst,
  Lauren P. Jackson, Margaret S. Robinson & Georg H. H. Borner
  https://doi.org/10.1038/s41467-018-06172-7

}

\examples{
  data(davies2018wt)
  davies2018wt
  pData(davies2018wt)
  exprs(davies2018wt)[1:3,1:3]

  library("pRoloc")
  plot2D(davies2018wt,, main = "Davies 2018 HeLa - wt")
}

\keyword{datasets}