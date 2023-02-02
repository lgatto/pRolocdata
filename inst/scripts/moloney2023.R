library("MSnbase")
library("pRoloc")

## function to add experimental data
addExperimentInfo <- function(date = "Summer 2020",
                              instrument = "Orbitrap Fusion Lumos Tribrid",
                              species = "Trypanosoma brucei") {
  experiment <- new("MIAPE",
                    lab = "Cambridge Centre for Proteomics (CCP)",
                    name = "Nicola Moloney",
                    contact = "Kathryn S. Lilley",
                    email = "k.s.lilley@bioc.cam.ac.uk",
                    samples = list(
                      species = species,
                      operator = "Nicola Moloney"
                    ),
                    title = "Mapping diversity in African trypanosomes using high resolution spatial proteomics",
                    abstract = "",
                    pubMedIds = "",
                    url = "",
                    instrumentModel = instrument,
                    instrumentManufacturer = "ThermoScientific",
                    ionSource = "",
                    analyser = "Orbitrap",
                    detectorType = "Orbitrap",
                    softwareName = "Mascot Search Engine",
                    collisionEnergy = "",
                    dateStamp = "Summer 2020",
  )
}

## populate the msnset slots
makeMSnSet <- function(filename_data,
                       filename_pdata,
                       ecols = c(2:34),
                       species = "Trypanosoma brucei") {
  
  # csv <- read.csv(filename_data)
  data <- readMSnSet2(file = filename_data, 
                      ecol = ecols,
                      skip = 0, 
                      fnames = 1, 
                      sep = "\t")
  
  ## remove redundant columns in fdata
  if (any(fvarLabels(data) == "X")) {
    ind <- which(fvarLabels(data) == "X")
    fData(data) <- fData(data)[, -ind]
  }
  
  ## remove "X" prefix to sampleNames
  sampleNames(data) <- sapply(strsplit(sampleNames(data), 
                                       split = "X"), "[[", 2)
  
  
  if (!missing(filename_pdata)) {
    pdat <- read.csv(filename_pdata)
    if(all(sampleNames(data) != pdat$X)) stop("data does not match")
    rownames(pdat) <- pdat$X
    pdat <- pdat[, -1]
  }
  
  ## define slots of the MSnSet
  .experiment <- addExperimentInfo(species = species)
  .process <- new("MSnProcess",
                  processing=c(
                    paste("Loaded on ",date(),".",sep=""),
                    paste("Normalised to sum of intensities. Tue Jun 08 23:36:04 2021"),
                    paste("Combined MSnSets. Tue Jun 08 23:36:04 2021"), 
                    paste("Subset MSnSets. Tue Jun 08 23:36:04 2021"),
                    paste("Removed features with more than 0 NAs. Tue Jun 08 23:36:04 2021 "),
                    paste("Dropped featureData's levels. Tue Jun 08 23:36:04 2021."),
                    paste("Added markers from 20210610_Tb_markers.csv. Tue Jun 22 12:15:59 2021"),
                    paste("Performed TAGM-MCMC prediction. Jul 12 11:03:24 2021")))
  .exprs <- exprs(data)
  .pData <- new("AnnotatedDataFrame", pdat)
  .fData <- new("AnnotatedDataFrame", fData(data))
  
  ## create msnset
  obj <- new("MSnSet",
             exprs = .exprs,
             phenoData = .pData,
             experimentData = .experiment,
             featureData = .fData)
  obj@processingData <- .process

  if (validObject(obj))
    return(obj)
}

moloneyTbBSF <- makeMSnSet(
  filename_data = "../extdata/20220304_TbBSF_n3_output.csv", 
  filename_pdata = "../extdata/20220304_TbBSF_pdata.csv",
  species = "African Trypanosoma brucei, bloodstream")

moloneyTbPCF <- makeMSnSet(
  filename_data = "../extdata/20220304_TbPCF_n3_output.csv", 
  filename_pdata = "../extdata/20220304_TbPCF_pdata.csv",
  species = "African Trypanosoma brucei, procyclic form")

moloneyTcBSF <- makeMSnSet(
  filename_data = "../extdata/20220304_TcBSF_n3_output.csv", 
  filename_pdata = "../extdata/20220304_TcBSF_pdata.csv",
  species = "African Trypanosoma congolense, bloodstream form")

moloneyTcPCF <- makeMSnSet(
  filename_data = "../extdata/20220304_TcPCF_n3_output.csv", 
  filename_pdata = "../extdata/20220304_TcPCF_pdata.csv",
  species = "African Trypanosoma congolense, procyclic form")


stopifnot(validObject(moloneyTbBSF))
save(moloneyTbBSF, file = "../../data/moloneyTbBSF", 
     compress = "xz", compression_level = 9)

stopifnot(validObject(moloneyTcBSF))
save(moloneyTbBSF, file = "../../data/moloneyTcBSF",
     compress = "xz", compression_level = 9)

stopifnot(validObject(moloneyTbPCF))
save(moloneyTbBSF, file = "../../data/moloneyTbPCF", 
     compress = "xz", compression_level = 9)

stopifnot(validObject(moloneyTcPCF))
save(moloneyTbBSF, file = "../../data/moloneyTcPCF",
     compress = "xz", compression_level = 9)
