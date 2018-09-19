library("MSnbase")
library("pRoloc")


csvfile <- "../../inst/extdata/Supplementary_Data_1_yeast2018.csv.gz"
#getEcols(csvfile, split = ",", n = 3)

hl <- readMSnSet2(file = csvfile, ecol = 10:49, skip = 2)
featureNames(hl) <- fData(hl)[, "Accession"]

## Experimental data to add
experiment <- new("MIAPE",
                  lab = "Cambridge Centre for Proteomics (CCP)",
                  name = "Daniel Nightingale",
                  contact = "Kathryn S. Lilley",
                  email = "k.s.lilley@bioc.cam.ac.uk",
                  samples = list(
                    species = "Saccharomyces cerevisiae",
                    operator = "Daniel Nightingale"
                  ),
                  title = "The subcellular organisation of Saccharomyces cerevisiae",
                  abstract = "Subcellular protein localisation is essential for the mechanisms that govern cellular homeostasis. The ability to understand processes leading to this phenomenon will therefore enhance our understanding of cellular function. Here we review recent developments in this field with regard to mass spectrometry, fluorescence microscopy and computational prediction methods. We highlight relative strengths and limitations of current methodologies focussing particularly on studies in the yeast Saccharomyces cerevisiae. We further present the first cell-wide spatial proteome map of S. cerevisiae, generated using hyperLOPIT, a mass spectrometry-based protein correlation profiling technique. We compare protein subcellular localisation assignments from this map, with two published fluorescence microscopy studies and show that confidence in localisation assignment is attained using multiple orthogonal methods that provide complementary data.",
                  pubMedIds = "",
                  url = "",
                  instrumentModel = "Q Exactive mass spectrometer",
                  instrumentManufacturer = "ThermoScientific",
                  ionSource = "",
                  analyser = "Orbitrap",
                  detectorType = "Orbitrap",
                  softwareName = "Mascot Search Engine",
                  collisionEnergy = "",
                  dateStamp = "Summer 2018"
)


## Expression data
e <- exprs(hl)


## Experiment info
pd <- data.frame(Replicate = rep(1:4, each = 10),
                 Tag = rep(c("X126", "X127C", "X127N", "X128C", "X128N",
                             "X129C", "X129N", "X130C", "X130N", "X131"), 4),
                 row.names=colnames(e))  
pd <- new("AnnotatedDataFrame", pd)

## feature data
fd <- fData(hl)
names(fd)[9] <- "predicted.location"
fd <- new("AnnotatedDataFrame", fd)

process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("Normalised to sum of intensities.")),
               normalised=TRUE)

yeast2018 <- new("MSnSet",
                 exprs = e,
                 phenoData = pd,
                 experimentData = experiment,
                 featureData = fd)
                          

yeast2018@processingData <- process
if (validObject(yeast2018)) 
  return (yeast2018)

save(yeast2018, file="../../data/yeast2018.RData",
     compress = "xz", compression_level = 9)
