library("MSnbase")
library("pRoloc")

csvfile <- "../../inst/extdata/courtland2020pepitdedata.csv"
csv <- read.csv(csvfile)

## There are 2 conditions in this dataset
courtland_control <- readMSnSet2(file = csvfile, ecol = grepEcols(csvfile, pattern  = "Control", split = ","), skip = 0)
courtland_mutant <- readMSnSet2(file = csvfile, ecol = grepEcols(csvfile, pattern  = "Mutant", split = ","), skip = 0)

## Experimental data to add
experiment <- new("MIAPE",
                  lab = "Soderling",
                  name = "Tyler",
                  contact = "tyler",
                  email = "twesleyb10@gmail.com",
                  samples = list(
                    species = "mouse",
                    operator = "Jamie courtland"
                  ),
                  title = "Swip Spatial Proteomics",
                  abstract = "Mutation of the WASH complex subunit, SWIP, is
	     implicated in human intellectual disability, but the cellular
	     etiology of this association is unknown. We identify the neuronal
	     WASH complex proteome, revealing a network of endosomal proteins.
	     To uncover how dysfunction of endosomal SWIP leads to disease, we
	     generate a mouse model of the human
	     WASHC4 c.3056C>G mutation.  Quantitative spatial
	     proteomics analysis of SWIP P1019R mouse brain
	     reveals that this mutation destabilizes the WASH complex and
	     uncovers significant  perturbations in both endosomal and lysosomal
	     pathways.  Cellular and histological analyses confirm that
	     SWIP P1019R results in  endo-lysosomal disruption
	     and uncover indicators of neurodegeneration. We find that
	     SWIP P1019R not only impacts cognition, but also
	     causes significant progressive motor deficits in mice.  Remarkably,
	     a retrospective analysis of SWIP P1019R patients
	     confirms motor deficits in humans. Combined, these findings support
	     the model that WASH complex destabilization, resulting from
	     SWIP P1019R, drives cognitive and motor
	     impairments via endo-lysosomal dysfunction in the brain.",
                  pubMedIds = "",
                  url = "",
                  instrumentModel = "Lumos",
                  instrumentManufacturer = "ThermoScientific",
                  ionSource = "",
                  analyser = "Orbitrap",
                  detectorType = "Orbitrap",
                  softwareName = "MaxQuant ",
                  collisionEnergy = "",
                  dateStamp = "3 January 2018"
)

## Expression data
e <- exprs(courtland_control)

## Experiment info
toName <- paste0(colnames(courtland_control)[1:21])
colnames(e) <- toName
pd <- data.frame(toName,
                 row.names=colnames(e))  
pd <- new("AnnotatedDataFrame", pd)

## feature data
fd <- rownames(e)
fd <- as.data.frame(csv[,1:8])
fd
rownames(fd) <- rownames(e)
fd <- new("AnnotatedDataFrame", fd)

process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("")),
               normalised=FALSE)

courtland_control <- new("MSnSet",
                 exprs = e,
                 phenoData = pd,
                 experimentData = experiment,
                 featureData = fd)

## summarise to protein level
courtland_control <- combineFeatures(courtland_control, method = "median", fcol = "Accession", na.rm = TRUE)
courtland_control <- MSnbase::normalise(courtland_control, method = "sum")
fData(courtland_control)$markers <- "unknown"

## Expression data
e <- exprs(courtland_mutant)

## Experiment info
toName <- paste0(colnames(courtland_mutant)[1:21])
colnames(e) <- toName
pd <- data.frame(toName,
                 row.names=colnames(e))  
pd <- new("AnnotatedDataFrame", pd)

## feature data
fd <- rownames(e)
fd <- as.data.frame(csv[,1:8])
fd
rownames(fd) <- rownames(e)
fd <- new("AnnotatedDataFrame", fd)

process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("")),
               normalised=FALSE)

courtland_mutant <- new("MSnSet",
                         exprs = e,
                         phenoData = pd,
                         experimentData = experiment,
                         featureData = fd)

## summarise to protein level
courtland_mutant <- combineFeatures(courtland_mutant, method = "median", fcol = "Accession", na.rm = TRUE)
courtland_mutant <- MSnbase::normalise(courtland_mutant, method = "sum")
fData(courtland_mutant)$markers <- "unknown"

save(courtland_control, file="../../data/courtland_control.rda",
     compress = "xz", compression_level = 9)

save(courtland_mutant, file="../../data/courtland_mutant.rda",
     compress = "xz", compression_level = 9)


