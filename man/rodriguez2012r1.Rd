\name{rodriguez2012r1}
\alias{rodriguez2012r1}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("rodriguez2012r1")}
\format{
  The format is:
Formal class 'MSnSet' [package "MSnbase"] with 9 slots
  ..@ experimentData   :Formal class 'MIAPE' [package "MSnbase"] with 30 slots
  .. .. ..@ title                   : chr ""
  .. .. ..@ url                     : chr ""
  .. .. ..@ abstract                : chr ""
  .. .. ..@ pubMedIds               : chr ""
  .. .. ..@ samples                 : list()
  .. .. ..@ preprocessing           : list()
  .. .. ..@ other                   : list()
  .. .. ..@ dateStamp               : chr(0) 
  .. .. ..@ name                    : chr ""
  .. .. ..@ lab                     : chr ""
  .. .. ..@ contact                 : chr ""
  .. .. ..@ email                   : chr ""
  .. .. ..@ instrumentModel         : chr(0) 
  .. .. ..@ instrumentManufacturer  : chr(0) 
  .. .. ..@ instrumentCustomisations: chr(0) 
  .. .. ..@ softwareName            : chr(0) 
  .. .. ..@ softwareVersion         : chr(0) 
  .. .. ..@ switchingCriteria       : chr(0) 
  .. .. ..@ isolationWidth          : num(0) 
  .. .. ..@ parameterFile           : chr(0) 
  .. .. ..@ ionSource               : chr(0) 
  .. .. ..@ ionSourceDetails        : chr(0) 
  .. .. ..@ analyser                : chr(0) 
  .. .. ..@ analyserDetails         : chr(0) 
  .. .. ..@ collisionGas            : chr(0) 
  .. .. ..@ collisionPressure       : num(0) 
  .. .. ..@ collisionEnergy         : chr(0) 
  .. .. ..@ detectorType            : chr(0) 
  .. .. ..@ detectorSensitivity     : chr(0) 
  .. .. ..@ .__classVersion__       :Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. ..@ .Data:List of 2
  .. .. .. .. .. ..$ : int [1:3] 1 0 0
  .. .. .. .. .. ..$ : int [1:3] 0 2 2
  ..@ processingData   :Formal class 'MSnProcess' [package "MSnbase"] with 10 slots
  .. .. ..@ files            : chr(0) 
  .. .. ..@ processing       : chr(0) 
  .. .. ..@ merged           : logi(0) 
  .. .. ..@ cleaned          : logi(0) 
  .. .. ..@ removedPeaks     : chr(0) 
  .. .. ..@ smoothed         : logi(0) 
  .. .. ..@ trimmed          : num(0) 
  .. .. ..@ normalised       : logi(0) 
  .. .. ..@ MSnbaseVersion   : chr "1.99.1"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 0 1 3
  ..@ qual             :'data.frame':	0 obs. of  0 variables
Formal class 'data.frame' [package "methods"] with 4 slots
  .. .. ..@ .Data    : list()
  .. .. ..@ names    : chr(0) 
  .. .. ..@ row.names: int(0) 
  .. .. ..@ .S3Class : chr "data.frame"
  ..@ assayData        :<environment: 0x10188f88> 
  ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	2 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr [1:2] NA NA
  .. .. ..@ data             :'data.frame':	11 obs. of  2 variables:
  .. .. .. ..$ fraction: num [1:11] 2 3 4 5 6 7 8 9 10 11 ...
  .. .. .. ..$ rep     : num [1:11] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	7 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr [1:7] NA NA NA NA ...
  .. .. ..@ data             :'data.frame':	2215 obs. of  7 variables:
  .. .. .. ..$ Protein Name              : chr [1:2215] "Double-stranded RNA-binding protein Staufen homolog 1" "GTPase NRas" "Splicing factor, arginine/serine-rich 4" "Low-density lipoprotein receptor" ...
  .. .. .. ..$ Genes                     : chr [1:2215] "STAU1;STAU;RP3-470L14.2-007;RP3-470L14.2-006" "NRAS;HRAS1;HRAS" "SFRS4;SRP75" "LDLR" ...
  .. .. .. ..$ Peptides identified       : num [1:2215] 28 5 8 2 60 4 11 3 18 6 ...
  .. .. .. ..$ Unique peptides identified: num [1:2215] 28 2 5 2 60 4 11 3 18 6 ...
  .. .. .. ..$ Sequence Coverage (%)     : num [1:2215] 37.4 36.5 15.6 3 65.8 9.2 20.4 22.7 31 24.6 ...
  .. .. .. ..$ Molecular Weight (kDa)    : num [1:2215] 63.2 21.2 56.7 95.4 99.3 ...
  .. .. .. ..$ markers                   : chr [1:2215] "unknown" "unknown" "unknown" "unknown" ...
  .. .. ..@ dimLabels        : chr [1:2] "featureNames" "featureColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ annotation       : chr(0) 
  ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr(0) 
  .. .. ..@ data             :'data.frame':	11 obs. of  0 variables
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. ..@ .Data:List of 5
  .. .. .. ..$ : int [1:3] 3 3 1
  .. .. .. ..$ : int [1:3] 2 33 2
  .. .. .. ..$ : int [1:3] 1 3 0
  .. .. .. ..$ : int [1:3] 0 1 1
  .. .. .. ..$ : int [1:3] 0 4 0
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(rodriguez2012r1)
## maybe str(rodriguez2012r1) ; plot(rodriguez2012r1) ...
}
\keyword{datasets}
