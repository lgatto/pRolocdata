library("MSnbase")
library("pRoloc")


csvfile <- "../../inst/extdata/hl-geladaki-2018.csv.gz"
#getEcols(csvfile)

hl <- readMSnSet2(file = csvfile, ecol = 2:58, fnames = 1)


## Experimental data to add
experiment <- new("MIAPE",
                  lab = "Cambridge Centre for Proteomics (CCP)",
                  name = "Kathryn S. Lilley",
                  contact = "Kathryn S. Lilley",
                  email = "k.s.lilley@bioc.cam.ac.uk",
                  samples = list(
                    species = "Human",
                    tissue = "U-2 OS human osteosarcoma cells",
                    operator = "Katerina Geladaki"
                  ),
                  title = "LOPIT-DC: A simpler approach to high-resolution spatial proteomics",
                  abstract = "Hyperplexed Localisation of Organelle Proteins by Isotope Tagging (hyperLOPIT) is a well-established method for studying protein subcellular localisation in complex biological samples. As a simpler alternative we developed a second workflow named Localisation of Organelle Proteins by Isotope Tagging after Differential ultraCentrifugation (LOPIT-DC) which is faster and less resource-intensive. We present the most comprehensive high-resolution mass spectrometry-based human dataset to date and deliver a flexible set of subcellular proteomics protocols for sample preparation and data analysis. For the first time, we methodically compare these two different mass spectrometry-based spatial proteomics methods within the same study and also apply QSep, the first tool that objectively and robustly quantifies subcellular resolution in spatial proteomics data. Using both approaches we highlight suborganellar resolution and isoform-specific subcellular niches as well as the locations of large protein complexes and proteins involved in signalling pathways which play important roles in cancer and metabolism. Finally, we showcase an extensive analysis of the multilocalising proteome identified via both methods",
                  pubMedIds = "",
                  url = "",
                  instrumentModel = "Orbitrap Fusion Lumos Tribrid",
                  instrumentManufacturer = "ThermoScientific",
                  ionSource = "ESI",
                  analyser = "Orbitrap",
                  detectorType = "Orbitrap",
                  softwareName = "Mascot Search Engine",
                  collisionEnergy = "",
                  dateStamp = "Spring 2018"
)
  

## Expression data
e <- exprs(hl)

ll1 <- strsplit(colnames(e), "_")
ll2 <- strsplit(colnames(e), "set")
f1 <- function(x) x[[1]]
f2 <- function(x) x[[2]]

## Experiment info
pd <- data.frame(Replicate = c(rep(1, 18), rep(2, 19), rep(3, 20)),
                 Set = sapply(ll2, f2),
                 TMT.Reagent = sapply(ll1, f1),
                 Acquisiton.Method = "MS3",
                 row.names=colnames(e))  
pd <- new("AnnotatedDataFrame", pd)

## feature data
fd <- fData(hl)[3:7]
names(fd) <- c("description", "markers", "svm", 
               "svm.scores", "final.assignment")
  
fd <- new("AnnotatedDataFrame", fd)

process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("Normalised to sum of intensities.")),
               normalised=TRUE)

hyperLOPITU2OS2018 <- new("MSnSet",
                          exprs = e,
                          phenoData = pd,
                          experimentData = experiment,
                          featureData = fd)

hyperLOPITU2OS2018@processingData <- process
if (validObject(hyperLOPITU2OS2018)) 
  return (hyperLOPITU2OS2018)

## Now create LOPIT-DC dataset
csvfile <- "../../inst/extdata/dc-geladaki-2018.csv.gz"
getEcols(csvfile, split = ",")

dc <- readMSnSet2(file = csvfile, ecol = 2:31, fnames = 1)


## Experimental data to add
experiment <- new("MIAPE",
                  lab = "Cambridge Centre for Proteomics (CCP)",
                  name = "Kathryn S. Lilley",
                  contact = "Kathryn S. Lilley",
                  email = "k.s.lilley@bioc.cam.ac.uk",
                  samples = list(
                    species = "Human",
                    tissue = "U-2 OS human osteosarcoma cells",
                    operator = "Nina Kocevar Britovsek"
                  ),
                  title = "LOPIT-DC: A simpler approach to high-resolution spatial proteomics",
                  abstract = "Hyperplexed Localisation of Organelle Proteins by Isotope Tagging (hyperLOPIT) is a well-established method for studying protein subcellular localisation in complex biological samples. As a simpler alternative we developed a second workflow named Localisation of Organelle Proteins by Isotope Tagging after Differential ultraCentrifugation (LOPIT-DC) which is faster and less resource-intensive. We present the most comprehensive high-resolution mass spectrometry-based human dataset to date and deliver a flexible set of subcellular proteomics protocols for sample preparation and data analysis. For the first time, we methodically compare these two different mass spectrometry-based spatial proteomics methods within the same study and also apply QSep, the first tool that objectively and robustly quantifies subcellular resolution in spatial proteomics data. Using both approaches we highlight suborganellar resolution and isoform-specific subcellular niches as well as the locations of large protein complexes and proteins involved in signalling pathways which play important roles in cancer and metabolism. Finally, we showcase an extensive analysis of the multilocalising proteome identified via both methods",
                  pubMedIds = "",
                  url = "",
                  instrumentModel = "Orbitrap Fusion Lumos Tribrid",
                  instrumentManufacturer = "ThermoScientific",
                  ionSource = "ESI",
                  analyser = "Orbitrap",
                  detectorType = "Orbitrap",
                  softwareName = "Mascot Search Engine",
                  collisionEnergy = "",
                  dateStamp = "Spring 2018"
)


## Expression data
e <- exprs(dc)

ll <- strsplit(colnames(e), "rep")
f1 <- function(x) x[[1]]
f2 <- function(x) x[[2]]

## Experiment info
pd <- data.frame(Replicate = sapply(ll, f2),
                 Fraction = sapply(ll, f1),
                 TMT.Reagent = rep(c("126", "127N", "127C", "128N",
                                   "128C", "129N", "129C", "130N",
                                   "130C", "131N"), 3))
rownames(pd) <- colnames(e)

pd <- new("AnnotatedDataFrame", pd)

## feature data
fd <- fData(dc)[3:11]
fd$markers <- fData(dc)$markers_10_classes
fd$final.assignment <- fData(dc)$curated_svm_predictions_10_classes
fd <- new("AnnotatedDataFrame", fd)

process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("Normalised to sum of intensities.")),
               normalised=TRUE)

lopitdcU2OS2018 <- new("MSnSet",
                       exprs = e,
                       phenoData = pd,
                       experimentData = experiment,
                       featureData = fd)

lopitdcU2OS2018@processingData <- process
if (validObject(lopitdcU2OS2018)) 
  return (lopitdcU2OS2018)

save(hyperLOPITU2OS2018, file="../../data/hyperLOPITU2OS2018.RData",
     compress = "xz", compression_level = 9)
save(lopitdcU2OS2018, file="../../data/lopitdcU2OS2018.RData",
     compress = "xz", compression_level = 9)