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
                     abstract = "",
                     pubMedIds = "",
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
                       TMT.Reagent = rep(colnames(exprs)[1:10], each = length(reps)),
                       AcquisitonMethod = method,
                       row.names=colnames(exprs))
  return(.pData)
}

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
l <- c("EntryName", "ProteinDescription", "Peptides.r1", "Peptides.r2", 
       "PSMs.r1", "PSMs.r2", l[7:16], paste0(l[7:16], ".r2"), l[27:44])
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

hyperLOPIT2015ms3r1 <- makeFusion(f1, 1, "MS3")
hyperLOPIT2015ms3r2 <- makeFusion(f2, 2, "MS3")
hyperLOPIT2015ms3r3 <- makeFusion(f3, 1, "MS3", date = "Summer 2015")
hyperLOPIT2015ms2 <- makeFusion(f4, 1, "MS2", instrument = "Q Exactive") 
  
hyperLOPIT2015ms3r1 <- addMarkers(hyperLOPIT2015ms3r1, markers = mrk, verbose = FALSE)
hyperLOPIT2015ms3r2 <- addMarkers(hyperLOPIT2015ms3r2, markers = mrk, verbose = FALSE)
hyperLOPIT2015ms3r3 <- addMarkers(hyperLOPIT2015ms3r3, markers = mrk, verbose = FALSE)
hyperLOPIT2015ms2 <- addMarkers(hyperLOPIT2015ms2, markers = mrk, verbose = FALSE)

## Update fvarLabel description
fvarMetadata(hyperLOPIT2015ms3r1)$labelDescription[6]  <- 
  fvarMetadata(hyperLOPIT2015ms3r2)$labelDescription[6] <- 
  fvarMetadata(hyperLOPIT2015ms3r3)$labelDescription[6] <- 
  "Marker set, curated by AC and CMM, covering protein subcellular localizations to 14 subcellular compartments."

save(hyperLOPIT2015,file="../../data/hyperLOPIT2015.RData", compress = "xz", compression_level = 9)
save(hyperLOPIT2015ms3r1,file="../../data/hyperLOPIT2015ms3r1.RData", compress = "xz", compression_level = 9)
save(hyperLOPIT2015ms3r2,file="../../data/hyperLOPIT2015ms3r2.RData", compress = "xz", compression_level = 9)
save(hyperLOPIT2015ms3r3,file="../../data/hyperLOPIT2015ms3r3.RData", compress = "xz", compression_level = 9)
save(hyperLOPIT2015ms2,file="../../data/hyperLOPIT2015ms2.RData", compress = "xz", compression_level = 9)
