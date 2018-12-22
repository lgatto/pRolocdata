library("MSnbase")
library("pRoloc")


csvfile <- "../../inst/extdata/krahmer2018pcp.csv"
csvfileFeature <- "../../inst/extdata/krahmer2018pcpFeature.csv"
pcpcsv <- read.csv(csvfile)
featureData <- read.csv(csvfileFeature)
rownames(featureData) <- featureData[,1]
#getEcols(csvfile, split = ",", n = 3)

## There are 3 experiments in this dataset Low fat diet (LFD), HFD3 and HFD12
pcp <- readMSnSet2(file = csvfile, ecol = 16:213, skip = 2, fnames = 1)

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
e <- exprs(pcp)

## Experiment info
toName <- paste0(colnames(pcpcsv)[16:213])
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
                 paste("No Normalisation")),
               normalised=FALSE)

krahmer2018pcp <- new("MSnSet",
                 exprs = e,
                 phenoData = pd,
                 experimentData = experiment,
                 featureData = fd)


krahmer2018pcp@processingData <- process
stopifnot(validObject(krahmer2018pcp))

save(krahmer2018pcp, file="../../data/krahmer2018pcp.rda",
     compress = "xz", compression_level = 9)
