library("MSnbase")
library("pRoloc")


csvfile <- "../../inst/extdata/krahmer2018PhosphoPcp.csv"
csvfileFeature <- "../../inst/extdata/krahmer2018PhosphoPcpFeature.csv"
Phosphopcpcsv <- read.csv(csvfile)
featureData <- read.csv(csvfileFeature)
rownames(featureData) <- make.unique(as.character(featureData[,7]))
rownames(Phosphopcpcsv) <- make.unique(as.character(Phosphopcpcsv[,7]))
#getEcols(csvfile, split = ",", n = 3)

## There are 3 experiments in this dataset Low fat diet (LFD), HFD3 and HFD12
Phosphopcp <- readMSnSet2(file = csvfile, ecol = 20:195, skip = 0, fnames = 8)
rownames(Phosphopcp) <- rownames(Phosphopcpcsv)## match rownames

## Experimental data to add
experiment <- new("MIAPE",
                  lab = "Proteomics and Signal Transduction, Max-Planck Institute of Biochemistry",
                  name = "Natalie Krahmer",
                  contact = "Matthias Mann",
                  email = "mmann@biochem.mpg.de",
                  samples = list(
                    species = "Mouse Liver",
                    operator = "Natalie Krahmer"
                  ),
                  title = "Organellar Proteomics and Phospho-Proteomics Reveal Subcellular Reorganization in Diet-Induced Hepatic Steatosis",
                  abstract = "Lipid metabolism is highly compartmentalized between cellular organelles that dynamically adapt their compositions and interactions in response to metabolic challenges. Here, we investigate how diet-induced hepatic lipid accumulation, observed in non-alcoholic fatty liver disease (NAFLD), affects protein localization, organelle organization, and protein phosphorylation in vivo. We develop a mass spectrometric workflow for protein and phosphopeptide correlation profiling to monitor levels and cellular distributions of 6,000 liver proteins and 16,000 phosphopeptides during development of steatosis. Several organelle contact site proteins are targeted to lipid droplets (LDs) in steatotic liver, tethering organelles orchestrating lipid metabolism. Proteins of the secretory pathway dramatically redistribute, including the mis-localization of the COPI complex and sequestration of the Golgi apparatus at LDs. This correlates with reduced hepatic protein secretion. Our systematic in vivo analysis of subcellular rearrangements and organelle-specific phosphorylation reveals how nutrient overload leads to organellar reorganization and cellular dysfunction.",
                  pubMedIds = "",
                  url = "",
                  instrumentModel = "Q Exactive HF-X",
                  instrumentManufacturer = "ThermoScientific",
                  ionSource = "",
                  analyser = "Orbitrap",
                  detectorType = "Orbitrap",
                  softwareName = "MaxQuant ",
                  collisionEnergy = "",
                  dateStamp = "22 October 2018"
)

## Expression data
e <- exprs(Phosphopcp)

## Experiment info
toName <- paste0(colnames(Phosphopcpcsv)[20:195])
colnames(e) <- toName
pd <- data.frame(toName,
                 row.names=colnames(e))
pd <- new("AnnotatedDataFrame", pd)

## feature data
fd <- featureData[rownames(e),1:22]
colnames(fd)[1] <- "Ascession" #remove non ascii name
fd <- new("AnnotatedDataFrame", fd)

process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("Normalisation performed manually")),
               normalised=FALSE)

krahmer2018Phosphopcp <- new("MSnSet",
                      exprs = e,
                      phenoData = pd,
                      experimentData = experiment,
                      featureData = fd)

## Normalise
krahmer2018Phosphopcp <- normalise(krahmer2018Phosphopcp, method = "sum")

## Marker column
krahmer2018Phosphopcp <- fDataToUnknown(krahmer2018Phosphopcp, fcol = "Organelle_")
fd <- fData(krahmer2018Phosphopcp)
fd$markers <- fd$Organelle
fData(krahmer2018Phosphopcp) <- fd
krahmer2018Phosphopcp <- normalise(krahmer2018Phosphopcp, method = "sum")
#plot2D(krahmer2018Phosphopcp[,1:88])
#plot2D(krahmer2018Phosphopcp[,89:176])

## Phenodata
pData(krahmer2018Phosphopcp)$fraction <- as.numeric(sub("^.+FR", "", sampleNames(krahmer2018Phosphopcp)))
pData(krahmer2018Phosphopcp)$replicate <- as.numeric(sub("^.+_", "", sub("_FR.+$", "", sampleNames(krahmer2018Phosphopcp))))
pData(krahmer2018Phosphopcp)$experiment <- c(rep("LFD", 88), rep("HFD12", 88))

stopifnot(length(pData(krahmer2018Phosphopcp)$experiment) == ncol(e)) # check columns and experiments match


krahmer2018Phosphopcp@processingData <- process
stopifnot(validObject(krahmer2018Phosphopcp))

save(krahmer2018Phosphopcp,
     file="../../data/krahmer2018phosphopcp.rda",
     compress = "xz", compression_level = 9)
