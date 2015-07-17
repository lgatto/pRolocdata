## Script for making Yloc featureSet MSnSet
library("MSnbase")
library("pRoloc")
library("pRolocdata")

makeYloc <- function(object, fileName, lopitCols) {
  csv <- read.csv(fileName, sep = "\t")
  csv <- csv[, -lopitCols]
  .ma <- match(featureNames(object), csv[,1])
  eset <- data.matrix(csv[.ma, 2:(ncol(csv) - 1)])
  rownames(eset) <- csv[.ma, 1]
  if (all(featureNames(object) != rownames(eset))) 
    stop("rownames do not match")
  fd <- new("AnnotatedDataFrame", fData(object))
  pd <- new("AnnotatedDataFrame", pData(object))
  object <- new("MSnSet", expr = eset, featureData = fd)
  if (!validObject(object)) stop("Not valid MSnSet")
  return(object)
}

data(E14TG2aS1)
E14TG2aS1 <- markerMSnSet(E14TG2aS1, "markers.tl")
f1 <- "../extdata/yloc_mouse_seqId_selected_features.csv"
lopitCols <- c(33:38)
E14TG2aS1yLoc <- makeYloc(E14TG2aS1, f1, lopitCols)
fData(E14TG2aS1yLoc) <- fData(E14TG2aS1yLoc)[, -c(10:15)] 

## Update fvarMetaData slots
fvarMetadata(E14TG2aS1yLoc)[, 1] <- c("UniProtKB accession number",
                                      "UniProt gene name",
                                      "Full protein description", 
                                      "Peptides", 
                                      "Peptide spectrum match", 
                                      "Localisation inferred from GO: Andy' output from his quickGO program",
                                      "Initial markers defined by Christoforou",
                                      "Hand curated updated marker list defined by Christoforou, Mulvey and Breckels",
                                      "Markers used for transfer learning in Breckels et al 2015, this is a subset of markers from using minMarkers function to set a minimum of 13 proteins per cluster")

## Add knntl results from Breckels et al (2015)
load("../extdata/tl-res/yloc-tl.rda")
experimentData(E14TG2aS1yLoc)@other$knntl <- yloc.opt

save(E14TG2aS1yLoc, file = "E14TG2aS1yLoc.rda")


