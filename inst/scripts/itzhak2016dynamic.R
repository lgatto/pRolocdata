library("MSnbase")
library("pRoloc")

csvfile1 <- "../../inst/extdata/silacDynamicItzhakCtrl.csv"
csvfile2 <- "../../inst/extdata/silacDynamicItzhakEGF.csv"
csv1 <- read.csv(csvfile1)
csv2 <- read.csv(csvfile2)
#getEcols(csvfile, split = ",", n = 3)

## There are 3 replicates in each dataset
helaCtrl <- readMSnSet2(file = csvfile1, ecol = 7:21, skip = 0, fnames = 2)
helaEgf <- readMSnSet2(file = csvfile2, ecol = 7:21, skip = 0, fnames = 2)

## Experimental data to add
experiment <- new("MIAPE",
                  lab = "Max Planck Institute of Biochemistry, Germany",
                  name = "Daniel Itzhak",
                  contact = "Georg Borner",
                  email = "borner@biochem.mpg.de",
                  samples = list(
                    species = "Human HeLa Cells",
                    operator = "Daniel Itzhak"
                  ),
                  title = "Global, quantitative and dynamic mapping of protein subcellular localization",
                  abstract = "Subcellular localization critically influences protein function, and cells 
                  control protein localization to regulate biological processes. We have developed and 
                  applied Dynamic Organellar Maps, a proteomic method that allows global mapping of protein 
                  translocation events. We initially used maps statically to generate a database with 
                  localization and absolute copy number information for over 8700 proteins from HeLa cells, 
                  approaching comprehensive coverage. All major organelles were resolved, with exceptional
                  prediction accuracy (estimated at >92%). Combining spatial and abundance information yielded
                  an unprecedented quantitative view of HeLa cell anatomy and organellar composition,
                  at the protein level. We subsequently demonstrated the dynamic capabilities of the
                  approach by capturing translocation events following EGF stimulation, which we
                  integrated into a quantitative model. Dynamic Organellar Maps enable the proteome-wide
                  analysis of physiological protein movements, without requiring any reagents specific
                  to the investigated process, and will thus be widely applicable in cell biology.",
                  pubMedIds = "",
                  url = "",
                  instrumentModel = "Q Exactive HF-X",
                  instrumentManufacturer = "ThermoScientific",
                  ionSource = "",
                  analyser = "Orbitrap",
                  detectorType = "Orbitrap",
                  softwareName = "MaxQuant ",
                  collisionEnergy = "",
                  dateStamp = "3 January 2018"
)

## Expression data
e1 <- exprs(helaCtrl)
e2 <- exprs(helaEgf)


## Experiment info
toName <- paste0(colnames(e1)[1:15])
colnames(e1) <- toName
pd <- data.frame(toName,
                 row.names=colnames(e1))  
pd <- new("AnnotatedDataFrame", pd)

toName <- paste0(colnames(e2)[1:15])
colnames(e2) <- toName
pd2 <- data.frame(toName,
                 row.names=colnames(e2))  
pd2 <- new("AnnotatedDataFrame", pd2)

## feature data
fd <- rownames(e1)
fd <- as.data.frame(fd)
fd$markers <- "unknown"
fd$markers <- as.character(csv1$Organellar.markers)
fd$markers[fd$markers == ""] <- "unknown"
rownames(fd) <- rownames(e1)
fd <- new("AnnotatedDataFrame", fd)

fd2 <- rownames(e2)
fd2 <- as.data.frame(fd2)
fd2$markers <- "unknown"
fd2$markers <- as.character(csv2$Organellar.markers)
fd2$markers[fd2$markers == ""] <- "unknown"
rownames(fd2) <- rownames(e2)
fd2 <- new("AnnotatedDataFrame", fd2)



process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("median Normalisation")),
               normalised=FALSE)

helaCtrl <- new("MSnSet",
                    exprs = e1,
                    phenoData = pd,
                    experimentData = experiment,
                    featureData = fd)

helaEgf <- new("MSnSet",
                exprs = e2,
                phenoData = pd2,
                experimentData = experiment,
                featureData = fd2)

## Phenodata
pData(helaCtrl)$fraction <- matrix(unlist(strsplit(sampleNames(helaCtrl), "[.]")), 15, 5, byrow = T)[,5]
pData(helaCtrl)$replicate <- rep(c(1,2,3), each = 5)

pData(helaEgf)$fraction <- matrix(unlist(strsplit(sampleNames(helaEgf), "[.]")), 15, 5, byrow = T)[,5]
pData(helaEgf)$replicate <- rep(c(1,2,3), each = 5)


stopifnot(length(pData(helaCtrl)$replicate) == ncol(e1)) # check columns and experiments match
stopifnot(length(pData(helaEgf)$replicate) == ncol(e2)) # check columns and experiments match


helaCtrl@processingData <- process
helaEgf@processingData <- process
stopifnot(validObject(helaCtrl))
stopifnot(validObject(helaEgf))

itzhak2016helaCtrl <- helaCtrl
itzhak2016helaEgf <- helaEgf

save(itzhak2016helaCtrl, file="../../data/itzhak2016helaCtrl.rda", compress = "xz", compression_level = 9)
save(itzhak2016helaEgf, file="../../data/itzhak2016helaEgf.rda", compress = "xz", compression_level = 9)
