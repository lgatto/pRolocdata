library("MSnbase")
library("pRoloc")
require("org.Hs.eg.db")
require("UniProt.ws")

csvfile <- "../../inst/extdata/HCC827rep2.csv"
csvfileMarker <- "../../inst/extdata/MarkersOrre2019.csv"
csv <- read.csv(csvfile)
makerData <- read.csv(csvfileMarker)
#getEcols(csvfile, split = ",", n = 3)

## There are 2 replicates in this dataset
HCC827rep2 <- readMSnSet2(file = csvfile, ecol = 2:11, skip = 0, fnames = 1)

## Convert rownames to Uniprot
genenames <- Rkeys(org.Hs.egSYMBOL)
key <- UniProt.ws::select(org.Hs.eg.db, genenames,"UNIPROT", "SYMBOL")
rownames(key) <- make.unique(key[, 1]) # get unique rownames
proteinnames <- key[rownames(HCC827rep2)[rownames(HCC827rep2) %in% key[, 1]], 2] # careful in wrong order

## Experimental data to add
experiment <- new("MIAPE",
                  lab = "Cancer Proteomics Mass Spectrometry, Karokinska",
                  name = "Lukas Orre",
                  contact = "Janne Lehtio",
                  email = "janne.lehtio@ki.se",
                  samples = list(
                    species = "Human lung cancer",
                    operator = "Lukas Orre"
                  ),
                  title = "SubCellBarCode: Proteome-wide Mapping of Protein Localization and Relocalization",
                  abstract = "Subcellular localization is a main determinant of protein function; however, a global view of cellular proteome organization remains relatively unexplored. We have developed a robust mass spectrometry-based analysis pipeline to generate a proteome-wide view of subcellular localization for proteins mapping to 12,418 individual genes across five cell lines. Based on more than 83,000 unique classifications and correlation profiling, we investigate the effect of alternative splicing and protein domains on localization, complex member co-localization, cell-type-specific localization, as well as protein relocalization after growth factor inhibition. Our analysis provides information about the cellular architecture and complexity of the spatial organization of the proteome; we show that the majority of proteins have a single main subcellular location, that alternative splicing rarely affects subcellular location, and that cell types are best distinguished by expression of proteins exposed to the surrounding environment. The resource is freely accessible via www.subcellbarcode.org.",
                  pubMedIds = "",
                  url = "",
                  instrumentModel = "(isoelectric focusing) Q Exactive HF-X",
                  instrumentManufacturer = "ThermoScientific",
                  ionSource = "",
                  analyser = "Orbitrap",
                  detectorType = "Orbitrap",
                  softwareName = "MaxQuant ",
                  collisionEnergy = "",
                  dateStamp = "3 January 2018"
)

## Expression data
e <- exprs(HCC827rep2)

## Experiment info
toName <- paste0(colnames(HCC827rep2)[1:10])
colnames(e) <- toName
pd <- data.frame(toName,
                 row.names=colnames(e))  
pd <- new("AnnotatedDataFrame", pd)

## feature data
fd <- rownames(e)
fd <- as.data.frame(fd)
rownames(makerData) <- makerData[,1]
fd$markers <- "unknown"
fd$markers[rownames(e) %in% makerData[,1]] <- as.character(makerData[rownames(e)[rownames(e) %in% makerData[,1]],2])
rownames(fd) <- rownames(e)
fd$protein <- NA
fd$protein[rownames(HCC827rep2) %in% key[, 1]] <- proteinnames
fd <- new("AnnotatedDataFrame", fd)

process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("median Normalisation")),
               normalised=FALSE)

orre2019hcc827rep2 <- new("MSnSet",
                          exprs = e,
                          phenoData = pd,
                          experimentData = experiment,
                          featureData = fd)

## Normalise
orre2019hcc827rep2 <- normalise(orre2019hcc827rep2, method = "sum")

## Map cluster to organelles (manually extracted)
fData(orre2019hcc827rep2)[,2][fData(orre2019hcc827rep2)[,2] == "mclust1"] <- c("Golgi, Endo/Lyosome")
fData(orre2019hcc827rep2)[,2][fData(orre2019hcc827rep2)[,2] == "mclust2"] <- c("ER, peroxisome")
fData(orre2019hcc827rep2)[,2][fData(orre2019hcc827rep2)[,2] == "mclust3"] <- c("ER, Mitochondrion")
fData(orre2019hcc827rep2)[,2][fData(orre2019hcc827rep2)[,2] == "mclust4"] <- c("Plasma Membrane")
fData(orre2019hcc827rep2)[,2][fData(orre2019hcc827rep2)[,2] == "mclust5"] <- c("Nucleosol, Ribosome")
fData(orre2019hcc827rep2)[,2][fData(orre2019hcc827rep2)[,2] == "mclust6"] <- c("Nucleus (speckles)")
fData(orre2019hcc827rep2)[,2][fData(orre2019hcc827rep2)[,2] == "mclust7"] <- c("Nucleus (Nucleolus)")
fData(orre2019hcc827rep2)[,2][fData(orre2019hcc827rep2)[,2] == "mclust8"] <- c("Nucleosol")
fData(orre2019hcc827rep2)[,2][fData(orre2019hcc827rep2)[,2] == "mclust9"] <- c("Cytosol_1")
fData(orre2019hcc827rep2)[,2][fData(orre2019hcc827rep2)[,2] == "mclust10"] <- c("Cytosol_2")
fData(orre2019hcc827rep2)[,2][fData(orre2019hcc827rep2)[,2] == "mclust11"] <- c("Cytosol_3")
fData(orre2019hcc827rep2)[,2][fData(orre2019hcc827rep2)[,2] == "mclust12"] <- c("Cytosol_4")
fData(orre2019hcc827rep2)[,2][fData(orre2019hcc827rep2)[,2] == "mclust13"] <- c("Cytosol_5")
fData(orre2019hcc827rep2)[,2][fData(orre2019hcc827rep2)[,2] == "mclust14"] <- c("Mitochondria matrix")
fData(orre2019hcc827rep2)[,2][fData(orre2019hcc827rep2)[,2] == "mclust15"] <- c("Mitochondria membrane")

plot2D(orre2019hcc827rep2, main = "HCC827rep2")
addLegend(orre2019hcc827rep2, where = "bottomleft", ncol = 1, cex = 0.9)

## Phenodata
pData(orre2019hcc827rep2)$fraction <- matrix(unlist(strsplit(sampleNames(orre2019hcc827rep2), "[.]")), 10, 4, byrow = T)[,1]
pData(orre2019hcc827rep2)$replicate <- matrix(unlist(strsplit(sampleNames(orre2019hcc827rep2), "[.]")), 10, 4, byrow = T)[,2]

stopifnot(length(pData(orre2019hcc827rep2)$replicate) == ncol(e)) # check columns and experiments match


orre2019hcc827rep2@processingData <- process
stopifnot(validObject(orre2019hcc827rep2))

save(orre2019hcc827rep2, file="../../data/orre2019hcc827rep2.rda",
     compress = "xz", compression_level = 9)