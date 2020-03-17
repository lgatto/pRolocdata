library("MSnbase")
library("pRoloc")

csvfileA1 <- "../../inst/extdata/Shin2019ExpA1_RAWData.csv"
csvfileA2 <- "../../inst/extdata/Shin2019ExpA2_RAWData.csv"
csvfileA3 <- "../../inst/extdata/Shin2019ExpA3_RAWData.csv"
csvfileB1 <- "../../inst/extdata/Shin2019ExpB1_RAWData.csv"
csvfileB2 <- "../../inst/extdata/Shin2019ExpB2_RAWData.csv"
csvfileB3 <- "../../inst/extdata/Shin2019ExpB3_RAWData.csv"
csvfileC1 <- "../../inst/extdata/Shin2019ExpC1_RAWData.csv"
csvfileC2 <- "../../inst/extdata/Shin2019ExpC2_RAWData.csv"
csvfileC3 <- "../../inst/extdata/Shin2019ExpC3_RAWData.csv"

markercsv <- read.csv(file = "../../inst/extdata/Shin2019newMarkers.csv")

csvfile <- c(csvfileC3, csvfileC2, csvfileC1, csvfileB3, csvfileB2, csvfileB1,
             csvfileA1, csvfileA2, csvfileA3)

csv <- read.csv(csvfile)


Shin2019MitoControlrep1 <- readMSnSet2(file = csvfile[7], ecol = 22:32, skip = 0, fnames = 1)
Shin2019MitoControlrep2 <- readMSnSet2(file = csvfile[8], ecol = 22:32, skip = 0, fnames = 1)
Shin2019MitoControlrep3 <- readMSnSet2(file = csvfile[9], ecol = 22:32, skip = 0, fnames = 1)
Shin2019MitoGcc88rep1 <- readMSnSet2(file = csvfile[1], ecol = 22:32, skip = 0, fnames = 1)
Shin2019MitoGcc88rep2 <- readMSnSet2(file = csvfile[2], ecol = 22:32, skip = 0, fnames = 1)
Shin2019MitoGcc88rep3 <- readMSnSet2(file = csvfile[3], ecol = 22:32, skip = 0, fnames = 1)
Shin2019MitoGol97rep1 <- readMSnSet2(file = csvfile[4], ecol = 22:32, skip = 0, fnames = 1)
Shin2019MitoGol97rep2 <- readMSnSet2(file = csvfile[5], ecol = 22:32, skip = 0, fnames = 1)
Shin2019MitoGol97rep3 <- readMSnSet2(file = csvfile[6], ecol = 22:32, skip = 0, fnames = 1)

Shin2019MitoControlrep1 <- filterNA(Shin2019MitoControlrep1)
Shin2019MitoControlrep2 <- filterNA(Shin2019MitoControlrep2)
Shin2019MitoControlrep3 <- filterNA(Shin2019MitoControlrep3)
Shin2019MitoGcc88rep1 <- filterNA(Shin2019MitoGcc88rep1)
Shin2019MitoGcc88rep2 <- filterNA(Shin2019MitoGcc88rep2)
Shin2019MitoGcc88rep3 <- filterNA(Shin2019MitoGcc88rep3)
Shin2019MitoGol97rep1 <- filterNA(Shin2019MitoGol97rep1)
Shin2019MitoGol97rep2 <- filterNA(Shin2019MitoGol97rep2)
Shin2019MitoGol97rep3 <- filterNA(Shin2019MitoGol97rep3)

Shin2019MitoControlrep1 <- updateSampleNames(Shin2019MitoControlrep1, 1)
Shin2019MitoControlrep2 <- updateSampleNames(Shin2019MitoControlrep2, 2)
Shin2019MitoControlrep3 <- updateSampleNames(Shin2019MitoControlrep3, 3)
Shin2019MitoGcc88rep1 <- updateSampleNames(Shin2019MitoGcc88rep1, 1)
Shin2019MitoGcc88rep2 <- updateSampleNames(Shin2019MitoGcc88rep2, 2)
Shin2019MitoGcc88rep3 <- updateSampleNames(Shin2019MitoGcc88rep3, 3)
Shin2019MitoGol97rep1 <- updateSampleNames(Shin2019MitoGol97rep1, 1)
Shin2019MitoGol97rep2 <- updateSampleNames(Shin2019MitoGol97rep2, 2)
Shin2019MitoGol97rep3 <- updateSampleNames(Shin2019MitoGol97rep3, 3)

Shin2019 <- MSnSetList(list(Shin2019MitoControlrep1, Shin2019MitoControlrep2, Shin2019MitoControlrep3,
                       Shin2019MitoGcc88rep1, Shin2019MitoGcc88rep2, Shin2019MitoGcc88rep3,
                       Shin2019MitoGol97rep1, Shin2019MitoGol97rep2, Shin2019MitoGol97rep3))

Shin2019 <- lapply(Shin2019, function(x) normalise(x, "sum"))

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
e <- lapply(Shin2019, exprs)

## Experiment info
pd <- list(length = length(Shin2019))
for (j in seq_along(Shin2019)) {
  
  toName <- paste0(colnames(e[[j]])[1:11])
  colnames(e[[j]]) <- toName
  pd[[j]] <- data.frame(toName,
                    row.names=colnames(e[[j]]))  
  pd[[j]] <- new("AnnotatedDataFrame", pd[[j]])
  
}

## feature data
fd <- list(length = length(Shin2019))
for (j in seq_along(Shin2019)) {
  fd[[j]] <- rownames(e[[j]])
  fd[[j]] <- as.data.frame(fd[[j]])
  markerdata <- as.data.frame(markercsv)
  rownames(markerdata) <- markercsv[,1]
  fd[[j]]$markers <- "unknown"
  rownames(fd[[j]]) <- rownames(Shin2019[[j]])
  fd[[j]][rownames(fd[[j]])[rownames(fd[[j]]) %in% rownames(markerdata)], "markers"] <-
    markerdata[rownames(fd[[j]])[rownames(fd[[j]]) %in% rownames(markerdata)],2]
  fd[[j]] <- new("AnnotatedDataFrame", fd[[j]])
}

process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("median Normalisation")),
               normalised=TRUE)

Shin2019MitoControlrep1 <- new("MSnSet",
                               exprs = e[[1]],
                               phenoData = pd[[1]],
                               experimentData = experiment,
                               featureData = fd[[1]])
Shin2019MitoControlrep2 <- new("MSnSet",
                               exprs = e[[2]],
                               phenoData = pd[[2]],
                               experimentData = experiment,
                               featureData = fd[[2]])
Shin2019MitoControlrep3 <- new("MSnSet",
                               exprs = e[[3]],
                               phenoData = pd[[3]],
                               experimentData = experiment,
                               featureData = fd[[3]])
Shin2019MitoGcc88rep1 <- new("MSnSet",
                             exprs = e[[4]],
                             phenoData = pd[[4]],
                             experimentData = experiment,
                             featureData = fd[[4]])
Shin2019MitoGcc88rep2 <- new("MSnSet",
                             exprs = e[[5]],
                             phenoData = pd[[5]],
                             experimentData = experiment,
                             featureData = fd[[5]])
Shin2019MitoGcc88rep3 <- new("MSnSet",
                             exprs = e[[6]],
                             phenoData = pd[[6]],
                             experimentData = experiment,
                             featureData = fd[[6]])
Shin2019MitoGol97rep1 <- new("MSnSet",
                             exprs = e[[7]],
                             phenoData = pd[[7]],
                             experimentData = experiment,
                             featureData = fd[[7]])
Shin2019MitoGol97rep2 <- new("MSnSet",
                             exprs = e[[8]],
                             phenoData = pd[[8]],
                             experimentData = experiment,
                             featureData = fd[[8]])
Shin2019MitoGol97rep3 <- new("MSnSet",
                             exprs = e[[9]],
                             phenoData = pd[[9]],
                             experimentData = experiment,
                             featureData = fd[[9]])

##
plot2D(Shin2019MitoControlrep1)


## Phenodata
pData(Shin2019MitoControlrep1)$fraction <- c("Nuc", "1000g", "3000g", "5000g", "9000g", "12000g", "15000g", "30000g", "79000g", "120000g", "SN")
pData(Shin2019MitoControlrep1)$replicate <- rep(c(1), each = 11)

pData(Shin2019MitoControlrep2)$fraction <- c("Nuc", "1000g", "3000g", "5000g", "9000g", "12000g", "15000g", "30000g", "79000g", "120000g", "SN")
pData(Shin2019MitoControlrep2)$replicate <- rep(c(2), each = 11)

pData(Shin2019MitoControlrep3)$fraction <- c("Nuc", "1000g", "3000g", "5000g", "9000g", "12000g", "15000g", "30000g", "79000g", "120000g", "SN")
pData(Shin2019MitoControlrep3)$replicate <- rep(c(3), each = 11)

pData(Shin2019MitoGcc88rep1)$fraction <- c("Nuc", "1000g", "3000g", "5000g", "9000g", "12000g", "15000g", "30000g", "79000g", "120000g", "SN")
pData(Shin2019MitoGcc88rep1)$replicate <- rep(c(1), each = 11)

pData(Shin2019MitoGcc88rep2)$fraction <- c("Nuc", "1000g", "3000g", "5000g", "9000g", "12000g", "15000g", "30000g", "79000g", "120000g", "SN")
pData(Shin2019MitoGcc88rep2)$replicate <- rep(c(2), each = 11)

pData(Shin2019MitoGcc88rep3)$fraction <- c("Nuc", "1000g", "3000g", "5000g", "9000g", "12000g", "15000g", "30000g", "79000g", "120000g", "SN")
pData(Shin2019MitoGcc88rep3)$replicate <- rep(c(3), each = 11)

pData(Shin2019MitoGol97rep1)$fraction <- c("Nuc", "1000g", "3000g", "5000g", "9000g", "12000g", "15000g", "30000g", "79000g", "120000g", "SN")
pData(Shin2019MitoGol97rep1)$replicate <- rep(c(1), each = 11)

pData(Shin2019MitoGol97rep2)$fraction <- c("Nuc", "1000g", "3000g", "5000g", "9000g", "12000g", "15000g", "30000g", "79000g", "120000g", "SN")
pData(Shin2019MitoGol97rep2)$replicate <- rep(c(2), each = 11)

pData(Shin2019MitoGol97rep3)$fraction <- c("Nuc", "1000g", "3000g", "5000g", "9000g", "12000g", "15000g", "30000g", "79000g", "120000g", "SN")
pData(Shin2019MitoGol97rep3)$replicate <- rep(c(3), each = 11)

## checks
stopifnot(length(pData(Shin2019MitoControlrep1)$replicate) == ncol(e[[1]])) # check columns and experiments match
stopifnot(length(pData(Shin2019MitoGcc88rep1)$replicate) == ncol(e[[4]])) # check columns and experiments match
stopifnot(length(pData(Shin2019MitoGol97rep3)$replicate) == ncol(e[[9]])) # check columns and experiments match

Shin2019MitoControlrep1@processingData <- process
Shin2019MitoControlrep2@processingData <- process
Shin2019MitoControlrep3@processingData <- process
Shin2019MitoGcc88rep1@processingData <- process
Shin2019MitoGcc88rep2@processingData <- process
Shin2019MitoGcc88rep3@processingData <- process
Shin2019MitoGol97rep1@processingData <- process
Shin2019MitoGol97rep2@processingData <- process
Shin2019MitoGol97rep3@processingData <- process

stopifnot(validObject(Shin2019MitoControlrep1))
stopifnot(validObject(Shin2019MitoGcc88rep1))
stopifnot(validObject(Shin2019MitoGol97rep1))

save(Shin2019MitoControlrep1, file="../../data/Shin2019MitoControlrep1.rda",
     compress = "xz", compression_level = 9)

save(Shin2019MitoControlrep2, file="../../data/Shin2019MitoControlrep2.rda",
     compress = "xz", compression_level = 9)

save(Shin2019MitoControlrep3, file="../../data/Shin2019MitoControlrep3.rda",
     compress = "xz", compression_level = 9)

save(Shin2019MitoGcc88rep1, file="../../data/Shin2019MitoGcc88rep1.rda",
     compress = "xz", compression_level = 9)

save(Shin2019MitoGcc88rep2, file="../../data/Shin2019MitoGcc88rep2.rda",
     compress = "xz", compression_level = 9)

save(Shin2019MitoGcc88rep3, file="../../data/Shin2019MitoGcc88rep3.rda",
     compress = "xz", compression_level = 9)

save(Shin2019MitoGol97rep1, file="../../data/Shin2019MitoGol97rep1.rda",
     compress = "xz", compression_level = 9)

save(Shin2019MitoGol97rep2, file="../../data/Shin2019MitoGol97rep2.rda",
     compress = "xz", compression_level = 9)

save(Shin2019MitoGol97rep3, file="../../data/Shin2019MitoGol97rep3.rda",
     compress = "xz", compression_level = 9)

