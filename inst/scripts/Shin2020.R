library("MSnbase")
library("pRoloc")


Shin2019MitoControl <- readRDS("../../inst/extdata/Shin2019MitoCon.rds")
Shin2019MitoGCC88 <- readRDS("../../inst/extdata/Shin2019MitoGCC88.rds")
Shin2019MitoGol97 <- readRDS("../../inst/extdata/Shin2019MitoGol97.rds")

Shin2019MitoControl <- Shin2019MitoControl$ExpA
Shin2019MitoGCC88 <- Shin2019MitoGCC88$ExpC
Shin2019MitoGol97 <- Shin2019MitoGol97$ExpB

## Experimental data to add
experiment <- new("MIAPE",
                  lab = "Medical Research Council, Laboratory for Molecular Biology",
                  name = "John Shin",
                  contact = "Sean Munro",
                  email = "sean@mrc-lmb.cam.ac.uk",
                  samples = list(
                    species = "Human HEK 293T",
                    operator = "John Shin"
                  ),
                  title = "Spatial proteomics defines the content of trafficking vesicles captured by golgin tethers",
                  abstract = "Intracellular traffic between compartments of the secretory and endocytic pathways is mediated by
                  vesicle-based carriers. The precise and complete proteomes of carriers destined for many
                  organelles are ill-defined because the vesicular intermediates are transient, low-abundance and
                  difficult to purify. Here, we combine vesicle relocalisation with organelle proteomics and Bayesian
                  analysis to define the content of different endosome-derived vesicles destined for the trans-Golgi
                  network (TGN). The golgin coiled-coil proteins golgin-97, golgin-245 and GCC88, shown previously
                  to capture endosome-derived vesicles at the TGN, were individually relocalised to mitochondria
                  and the content of subsequently re-routed vesicles was determined by organelle proteomics. Our
                  findings revealed 45 integral and 51 peripheral membrane proteins re-routed by golgin-97,
                  evidence for a distinct class of vesicles shared by golgin-97 and GCC88, and various cargoes
                  specific to individual golgins. These results illustrate a general strategy for analysing intracellular
                  sub-proteomes by combining acute cellular re-wiring with high-resolution spatial proteomics.",
                  pubMedIds = "",
                  url = "",
                  instrumentModel = "Q Exactive HF-X",
                  instrumentManufacturer = "ThermoScientific",
                  ionSource = "",
                  analyser = "Orbitrap",
                  detectorType = "Orbitrap",
                  softwareName = "MaxQuant ",
                  collisionEnergy = "",
                  dateStamp = "16 March 2020"
)


## Expression data
e1 <- exprs(Shin2019MitoControl)
e2 <- exprs(Shin2019MitoGCC88)
e3 <- exprs(Shin2019MitoGol97)


## Experiment info
toName <- paste0(colnames(e1)[1:33])
colnames(e1) <- toName
pd1 <- data.frame(toName,
                  row.names=colnames(e1))  
pd1 <- new("AnnotatedDataFrame", pd1)

## Experiment info
toName <- paste0(colnames(e2)[1:33])
colnames(e2) <- toName
pd2 <- data.frame(toName,
                  row.names=colnames(e2))  
pd2 <- new("AnnotatedDataFrame", pd2)


## Experiment info
toName <- paste0(colnames(e3)[1:33])
colnames(e3) <- toName
pd3 <- data.frame(toName,
                  row.names=colnames(e3))  
pd3 <- new("AnnotatedDataFrame", pd3)


## feature data
fd1 <- rownames(e1)
fd1 <- as.data.frame(fd1)
markerdata <- as.data.frame(fData(Shin2019MitoControl)$markers)
rownames(markerdata) <- rownames(fData(Shin2019MitoControl)$markers)
fd1$markers <- "unknown"
rownames(fd1) <- rownames(fData(Shin2019MitoControl))
fd1$markers[markerdata !=""] <- unlist(lapply(markerdata[markerdata !="",1], as.character))
fd1 <- new("AnnotatedDataFrame", fd1)

## feature data
fd2 <- rownames(e2)
fd2 <- as.data.frame(fd2)
markerdata <- as.data.frame(fData(Shin2019MitoGCC88)$markers)
rownames(markerdata) <- rownames(fData(Shin2019MitoGCC88)$markers)
fd2$markers <- "unknown"
rownames(fd2) <- rownames(fData(Shin2019MitoGCC88))
fd2$markers[markerdata !=""] <- unlist(lapply(markerdata[markerdata !="",1], as.character))
fd2 <- new("AnnotatedDataFrame", fd2)

## feature data
fd3 <- rownames(e3)
fd3 <- as.data.frame(fd3)
markerdata <- as.data.frame(fData(Shin2019MitoGol97)$markers)
rownames(markerdata) <- rownames(fData(Shin2019MitoGol97)$markers)
fd3$markers <- "unknown"
rownames(fd3) <- rownames(fData(Shin2019MitoGol97))
fd3$markers[markerdata !=""] <- unlist(lapply(markerdata[markerdata !="",1], as.character))
fd3 <- new("AnnotatedDataFrame", fd3)

process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("median Normalisation")),
               normalised=TRUE)

Shin2019MitoControl <- new("MSnSet",
                       exprs = e1,
                       phenoData = pd1,
                       experimentData = experiment,
                       featureData = fd1)

Shin2019MitoGCC88 <- new("MSnSet",
                       exprs = e2,
                       phenoData = pd2,
                       experimentData = experiment,
                       featureData = fd2)

Shin2019MitoGol97 <- new("MSnSet",
                    exprs = e3,
                    phenoData = pd3,
                    experimentData = experiment,
                    featureData = fd3)

##
plot2D(Shin2019MitoControl)


## Phenodata
pData(Shin2019MitoControl)$fraction <- c("Nuc", "1000g", "3000g", "5000g", "9000g", "12000g", "15000g", "30000g", "79000g", "120000g", "SN")
pData(Shin2019MitoControl)$replicate <- rep(c(1,2,3), each = 11)

pData(Shin2019MitoGCC88)$fraction <- c("Nuc", "1000g", "3000g", "5000g", "9000g", "12000g", "15000g", "30000g", "79000g", "120000g", "SN")
pData(Shin2019MitoGCC88)$replicate <- rep(c(1,2,3), each = 11)

pData(Shin2019MitoGol97)$fraction <- c("Nuc", "1000g", "3000g", "5000g", "9000g", "12000g", "15000g", "30000g", "79000g", "120000g", "SN")
pData(Shin2019MitoGol97)$replicate <- rep(c(1,2,3), each = 11)

## checks
stopifnot(length(pData(Shin2019MitoControl)$replicate) == ncol(e1)) # check columns and experiments match
stopifnot(length(pData(Shin2019MitoGCC88)$replicate) == ncol(e2)) # check columns and experiments match
stopifnot(length(pData(Shin2019MitoGol97)$replicate) == ncol(e1)) # check columns and experiments match

Shin2019MitoControl@processingData <- process
Shin2019MitoGCC88@processingData <- process
Shin2019MitoGol97@processingData <- process
stopifnot(validObject(Shin2019MitoControl))
stopifnot(validObject(Shin2019MitoGCC88))
stopifnot(validObject(Shin2019MitoGol97))

save(Shin2019MitoControl, file="../../data/Shin2019MitoControl.rda",
     compress = "xz", compression_level = 9)

save(Shin2019MitoGCC88, file="../../data/Shin2019MitoGCC88.rda",
     compress = "xz", compression_level = 9)

save(Shin2019MitoGol97, file="../../data/Shin2019MitoGol97.rda",
     compress = "xz", compression_level = 9)

