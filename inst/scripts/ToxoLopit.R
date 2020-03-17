## ToxoLopit data

library("MSnbase")
library("pRoloc")

csvfile <- "../../inst/extdata/TL123_complete.csv"
pDatacsv <- "../../inst/extdata/ToxoLOPIT_pData.csv"
csv <- read.csv(csvfile, sep = "")
csvpdata <- read.csv(pDatacsv)

Barylyuk2020ToxoLopit <- readMSnSet2(csvfile, ecol = 2:31, sep ="\t", fnames = 1)


## Experimental data to add
experiment <- new("MIAPE",
                  lab = "Department of Biochemistry, University of Cambridge, Cambridge, UK; Cambridge Centre for Proteomics, University of Cambridge, Cambridge, UK",
                  name = "Konstantin Barylyuk",
                  contact = "Konstantin Barylyuk, Ross F. Waller",
                  email = "kb601@cam.ac.uk, rfw26@cam.ac.uk",
                  samples = list(
                    species = "Toxoplasma gondii, strain RH",
                    operator = "Konstantin Barylyuk"
                  ),
                  title = "Whole-cell spatial proteome of Toxoplasma: molecular anatomy of an apicomplexan cell",
                  abstract = "",
                  pubMedIds = "",
                  url = "",
                  instrumentModel = "Orbitrap Fusion Lumos",
                  instrumentManufacturer = "ThermoScientific",
                  ionSource = "",
                  analyser = "Orbitrap",
                  detectorType = "Orbitrap",
                  softwareName = "Mascot, Proteome Discoverer",
                  collisionEnergy = "",
                  dateStamp = "16 March 2020"
)

## Expression data
e <- exprs(Barylyuk2020ToxoLopit)


## Experiment info
toName <- paste0(colnames(e)[1:30])
colnames(e) <- toName
pd <- data.frame(toName,
                  row.names=colnames(e), csvpdata)  
pd <- new("AnnotatedDataFrame", pd)


process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("VSN Normalisation")),
               normalised=TRUE)

Barylyuk2020ToxoLopit <- new("MSnSet",
                              exprs = e,
                              phenoData = pd,
                              experimentData = experiment,
                              featureData = new("AnnotatedDataFrame", fData(Barylyuk2020ToxoLopit)))

plot2D(Barylyuk2020ToxoLopit)

Barylyuk2020ToxoLopit@processingData <- process
stopifnot(validObject(Barylyuk2020ToxoLopit))

save(Barylyuk2020ToxoLopit, file="../../data/Barylyuk2020ToxoLopit.rda",
     compress = "xz", compression_level = 9)
