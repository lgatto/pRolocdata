library("MSnbase")
library("pRoloc")

## Function to extract assession no (an)
getAN <- function(x) {
  n <- vector("character", nrow(x))
  l <- vector("numeric", nrow(x))
  for (i in 1:nrow(x)) {
    spt1 <- unlist(strsplit(x[i,1], split="[", fixed=TRUE))
    l[i] <- length(spt1)
    if (length(spt1) > 2) {
      spt1 <- spt1[2:length(spt1)]
      an <- sapply(spt1, function(z) strsplit(z, split="]", fixed=TRUE))
      an <- unlist(lapply(an, function(z) z[1]), use.names=FALSE)
      aa <- paste(an[1], an[2], sep="|")
      if (length(an) > 2) {
        for (j in 3:length(an)) {
          aa <- paste(aa, an[j], sep="|")
        }
      }
    } else {
      aa <- unlist(strsplit(spt1, split="]", fixed=TRUE))[2]
    }
    n[i] <- aa
  }
  return(n)
}

## Experimental data to add
addExperimentInfo <- function(date = "Autumn 2013",
                              instrument = "Orbitrap Fusion Tribrid") {
  .experiment <- new("MIAPE",
                     lab = "Cambridge Centre for Proteomics (CCP)",
                     name = "Kathryn S. Lilley",
                     contact = "Kathryn S. Lilley",
                     email = "k.s.lilley@bioc.cam.ac.uk",
                     samples = list(
                       species = "Mouse",
                       tissue = "E14TG2a embryonic stem cells",
                       operator = c("Andy Christoforou",
                                    "Claire M. Mulvey")
                     ),
                     title = "A draft map of the mouse pluripotent stem cell spatial proteome",
                     abstract = "Knowledge of the subcellular distribution of proteins is vital for understanding cellular mechanisms. Capturing the subcellular proteome in a single experiment has proven challenging, with studies focusing on specific compartments or assigning proteins to subcellular niches with low resolution and/or accuracy. Here we introduce hyperLOPIT, a method that couples extensive fractionation, quantitative high-resolution accurate mass spectrometry with multivariate data analysis. We apply hyperLOPIT to a pluripotent stem cell population whose subcellular proteome has not been extensively studied. We provide localization data on over 5,000 proteins with unprecedented spatial resolution to reveal the organization of organelles, sub-organellar compartments, protein complexes, functional networks and steady-state dynamics of proteins and unexpected subcellular locations. The method paves the way for characterizing the impact of post-transcriptional and post-translational modification on protein location and studies involving proteome-level locational changes on cellular perturbation. An interactive open-source resource is presented that enables exploration of these data.",
                     pubMedIds = "26754106",
                     url = "",
                     instrumentModel = instrument,
                     instrumentManufacturer = "ThermoScientific",
                     ionSource = "ESI",
                     analyser = "Orbitrap",
                     detectorType = "Orbitrap",
                     softwareName = "Mascot Search Engine",
                     collisionEnergy = "",
                     dateStamp = date
  )
  return(.experiment)
}

## Experiment info
addPhenoData <- function(exprs, reps = c(1, 2), method = "MS3") {
  .pData <- data.frame(Replicate = rep(reps, each = 10),
                       TMT.Reagent = colnames(exprs),
                       Acquisiton.Method = method,
                       row.names=colnames(exprs))
  return(.pData)
}

finfo <- read.csv("../extdata/hyperLOPIT-SIData-fraction-info.csv",
                  row.names=1, skip = 1, stringsAsFactors = FALSE,
                  header = TRUE)
rownames(finfo) <- paste0("X", rownames(finfo))
  exp1 <- data.frame(finfo[, 1], finfo[, 4], stringsAsFactors = FALSE)
  exp2 <- data.frame(finfo[, 2], finfo[, 5], stringsAsFactors = FALSE)
  exp3 <- data.frame(finfo[, 3], finfo[, 6], stringsAsFactors = FALSE)
  colnames(exp1) <- colnames(exp2) <-
    colnames(exp3) <- c("Gradient.Fraction",
                        "Iodixonal.Density")
gradData <- vector("list", 3)
gradData[[1]] <- exp1
gradData[[2]] <- exp2
gradData[[3]] <- exp3

makeFusion <- function(filename, repNo = 1, method = "MS3",
                       date = "Autumn 2013",
                       instrument = "Orbitrap Fusion Tribrid") {
  csv <- read.csv(filename, row.names=1, skip = 1,
                  header = TRUE, stringsAsFactors=FALSE)
  EntryName <- getAN(csv)
  l <- colnames(csv)
  if (method != "MS3") {
    l[1:3] <- c("ProteinDescription", "Peptides", "PSMs")
  } else {
    l[1:4] <- c("ProteinDescription", "Peptides", "PSMs", "ProteinCoverage")
  }
  colnames(csv) <- l
  ind <- grep("X1", colnames(csv))
  .exprs <- csv[,ind]
  .exprs <- as.matrix(.exprs)
  .fData <- cbind(EntryName, csv[,-ind])
  .fData <- new("AnnotatedDataFrame", .fData)
  featuresinfo <- c("UniProt identifier for quantified protein group reported by Proteome Discoverer.",
                    "UniProt description for protein accession.",
                    "Number of quantified peptides. Only peptides that were unique to the a single protein group were used for quantification.",
                    "Number of quantified peptide-spectrum matches.",
                    "Protein Coverage %")
  if (method != "MS3") {
    .fData@varMetadata[,1] <- featuresinfo[-5]
  } else {
    .fData@varMetadata[,1] <- featuresinfo
  }
  .pData <- addPhenoData(.exprs, reps = repNo, method = method)
  if (method == "MS3") {
    .pData <- cbind(.pData, gradData[[repNo]])
  }
  .pData <- new("AnnotatedDataFrame", .pData)
  .experiment <- addExperimentInfo(date = date, instrument = instrument)
  .process <- new("MSnProcess",
                  processing=c(
                    paste("Loaded on ",date(),".",sep=""),
                    paste("Normalised to sum of intensities.")),
                  normalised=TRUE)
  obj <- new("MSnSet",
             exprs = .exprs,
             phenoData = .pData,
             experimentData = .experiment,
             featureData = .fData)
  obj@processingData <- .process
  if (validObject(obj))
    return (obj)
}

## Manually make intersect MSnSet from .csv
csv <- read.csv("../extdata/hyperLOPIT-SIData-ms3-rep12-intersect.csv.gz",
                row.names=1, header = TRUE, skip = 1, stringsAsFactors=FALSE)
l <- colnames(csv)
l <- c("entry.name", "protein.description", "Peptides.rep1", "Peptides.rep2",
       "PSMs.rep1", "PSMs.rep2", paste0(l[7:16], ".rep1"),
       paste0(l[7:16], ".rep2"), l[27:44])
l[grep("SVM.marker.set", l)] <- "markers"
colnames(csv) <- l
tokeep <- grep("X1", l)
.exprs <- csv[, tokeep]
.exprs <- as.matrix(.exprs)
.fData <- csv[, -tokeep]
uns <- which(.fData$Final.Localization.Assignment == "unclassified")
.fData$Final.Localization.Assignment[uns] <- "unknown"
.fData <- new("AnnotatedDataFrame", .fData)
info <- read.csv("../extdata/hyperLOPIT-SIData-info.csv.gz",      ## Add meta data
                 stringsAsFactors = FALSE)
metadata <- info[c(1:6, 17:34), 3]
metadata[3:6] <- c(paste0("Replicate 1: ", metadata[3]),
                   paste0("Replicate 2: ", metadata[3]),
                   paste0("Replicate 1: ", metadata[5]),
                   paste0("Replicate 2: ", metadata[5]))
.fData@varMetadata[,1] <- metadata
.pData <- addPhenoData(.exprs, reps = c(1, 2))
.pData[, 2] <- as.character(.pData[, 2])
.pData[, 2] <- sapply(1:length(.pData[, 2]),
                      function(z) strsplit(.pData[, 2][z], ".rep")[[1]][1])
.pData <- cbind(.pData, rbind(gradData[[1]], gradData[[2]]))
.pData <- new("AnnotatedDataFrame", .pData)
.experiment <- addExperimentInfo()
.process <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".",sep=""),
                  paste("Normalised to sum of intensities.")),
                normalised=TRUE)
hyperLOPIT2015 <- new("MSnSet",
           exprs = .exprs,
           phenoData = .pData,
           experimentData = .experiment,
           featureData = .fData)
hyperLOPIT2015@processingData <- .process

## Get markers for other replicates
mrk <- fData(markerMSnSet(hyperLOPIT2015, "markers"))$markers
mrk <- as.character(mrk)
ids <- featureNames(markerMSnSet(hyperLOPIT2015, "markers"))
names(mrk) <- ids

## Make MSnSets for each replicate
f1 <- "../extdata/hyperLOPIT-SIData-ms3-rep1.csv.gz"
f2 <- "../extdata/hyperLOPIT-SIData-ms3-rep2.csv.gz"
f3 <- "../extdata/hyperLOPIT-SIData-ms3-rep3.csv.gz"
f4 <- "../extdata/hyperLOPIT-SIData-ms2-rep1.csv.gz"

## make data
hyperLOPIT2015ms3r1 <- makeFusion(f1, 1, "MS3")
hyperLOPIT2015ms3r2 <- makeFusion(f2, 2, "MS3")
hyperLOPIT2015ms3r3 <- makeFusion(f3, 3, "MS3", date = "Summer 2015")
hyperLOPIT2015ms2 <- makeFusion(f4, 1, "MS2", instrument = "Q Exactive")

## add markers
hyperLOPIT2015ms3r1 <- addMarkers(hyperLOPIT2015ms3r1, markers = mrk, verbose = FALSE)
hyperLOPIT2015ms3r2 <- addMarkers(hyperLOPIT2015ms3r2, markers = mrk, verbose = FALSE)
hyperLOPIT2015ms3r3 <- addMarkers(hyperLOPIT2015ms3r3, markers = mrk, verbose = FALSE)
hyperLOPIT2015ms2 <- addMarkers(hyperLOPIT2015ms2, markers = mrk, verbose = FALSE)

## Update fvarLabel description
fvarMetadata(hyperLOPIT2015ms3r1)$labelDescription[6]  <-
  fvarMetadata(hyperLOPIT2015ms3r2)$labelDescription[6] <-
  fvarMetadata(hyperLOPIT2015ms3r3)$labelDescription[6] <-
  "Marker set, curated by AC and CMM, covering protein subcellular localizations to 14 subcellular compartments."

## Change fvarLabels
fvarLabels(hyperLOPIT2015) <- tolower(fvarLabels(hyperLOPIT2015))
fvarLabels(hyperLOPIT2015)[c(13, 14, 15, 20)] <- c("svm.top.quartile", "final.assignment", "first.evidence", "signalling.cascades")

## Update phenoData
pData(hyperLOPIT2015)[, 4] <- as.character(pData(hyperLOPIT2015)[, 4])
#pData(hyperLOPIT2015)[, 5] <- as.numeric(pData(hyperLOPIT2015)[, 5])

## Add unknown instead of "" for data plotting, change unclassified to unknown for completeness
hyperLOPIT2015 <- fDataToUnknown(hyperLOPIT2015, fcol = "svm.top.quartile", from = "unclassified")
aa <- fvarLabels(hyperLOPIT2015)
for (i in 16:24) {
  hyperLOPIT2015 <- fDataToUnknown(hyperLOPIT2015, fcol = aa[i])
}

## Rename markers to markers2015 and add new markers
fData(hyperLOPIT2015)$markers2015 <- fData(hyperLOPIT2015)$markers

toadd.ecm <- c("P47877", "P27808", "Q8C407")
toadd.end <- "Q62351"

fData(hyperLOPIT2015)[toadd.ecm, "markers"] <- "Extracellular matrix"
fData(hyperLOPIT2015)[toadd.end, "markers"] <- "Endosome"

## Adding TAGM allocations. The data.frames were generated by Olly
## during the development and implementation of the method
load("../extdata/hlMAPresults20180521.rda")  ## hlMAPresults20180521
load("../extdata/hlMCMCresults20180521.rda") ## hlMCMCresults20180521
stopifnot(identical(rownames(hlMAPresults20180521), featureNames(hyperLOPIT2015)))
stopifnot(identical(rownames(hlMCMCresults20180521), featureNames(hyperLOPIT2015)))

## change type from character to numeric
hlMCMCresults20180521$tagm.mcmc.probability <- as.numeric(hlMCMCresults20180521$tagm.mcmc.probability)
## Change convention for outlier so that markers have outlier of 0.
hlMCMCresults20180521$tagm.mcmc.outlier <- 1 - hlMCMCresults20180521$tagm.mcmc.outlier

fData(hyperLOPIT2015)$TAGM <-
                        cbind(hlMAPresults20180521,
                              hlMCMCresults20180521)




stopifnot(validObject(hyperLOPIT2015))

save(hyperLOPIT2015, file = "../../data/hyperLOPIT2015.RData", compress = "xz", compression_level = 9)
save(hyperLOPIT2015ms3r1, file = "../../data/hyperLOPIT2015ms3r1.RData", compress = "xz", compression_level = 9)
save(hyperLOPIT2015ms3r2, file = "../../data/hyperLOPIT2015ms3r2.RData", compress = "xz", compression_level = 9)
save(hyperLOPIT2015ms3r3, file = "../../data/hyperLOPIT2015ms3r3.RData", compress = "xz", compression_level = 9)
save(hyperLOPIT2015ms2, file = "../../data/hyperLOPIT2015ms2.RData", compress = "xz", compression_level = 9)


library("SummarizedExperiment")
hyperLOPIT2015_se <- as(hyperLOPIT2015, "SummarizedExperiment")
save(hyperLOPIT2015_se, file = "../../data/hyperLOPIT2015_se.rda", compress = "xz", compression_level = 9)
