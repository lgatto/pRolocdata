library("MSnbase")
library("pRoloc")


csvfile <- "../../inst/extdata/Daviesetal.csv"
csv <- read.csv(csvfile)
#getEcols(csvfile, split = ",", n = 3)

## There are 2 replicates in this dataset
davies2018ap4b1 <- readMSnSet2(file = csvfile, ecol = 6:15, skip = 0, fnames = 1)
davies2018ap4e1 <- readMSnSet2(file = csvfile, ecol = 16:25, skip = 0, fnames = 1)
davies2018wt <- readMSnSet2(file = csvfile, ecol = 26:35, skip = 0, fnames = 1)


## Experimental data to add
experiment <- new("MIAPE",
                  lab = "Cambridge Institute for Medical Research",
                  name = "Alexandra K. Davies",
                  contact = "Margaret S. Robinson",
                  email = "msr12@cam.ac.uk",
                  samples = list(
                    species = "Human HeLa",
                    operator = "Alexandra K. Davies"
                  ),
                  title = "AP-4 vesicles contribute to spatial control of autophagy via RUSC-dependent peripheral delivery of ATG9A",
                  abstract = "Adaptor protein 4 (AP-4) is an ancient membrane trafficking complex, whose function has largely remained elusive. In humans, AP-4 deficiency causes a severe neurological disorder of unknown aetiology. We apply unbiased proteomic methods, including 'Dynamic Organellar Maps', to find proteins whose subcellular localisation depends on AP-4. We identify three transmembrane cargo proteins, ATG9A, SERINC1 and SERINC3, and two AP-4 accessory proteins, RUSC1 and RUSC2. We demonstrate that AP-4 deficiency causes missorting of ATG9A in diverse cell types, including patient-derived cells, as well as dysregulation of autophagy. RUSC2 facilitates the transport of AP-4-derived, ATG9A-positive vesicles from the trans-Golgi network to the cell periphery. These vesicles cluster in close association with autophagosomes, suggesting they are the 'ATG9A reservoir' required for autophagosome biogenesis. Our study uncovers ATG9A trafficking as a ubiquitous function of the AP-4 pathway. Furthermore, it provides a potential molecular pathomechanism of AP-4 deficiency, through dysregulated spatial control of autophagy.",
                  pubMedIds = "",
                  url = "",
                  instrumentModel = "Q Exactive HF-X",
                  instrumentManufacturer = "ThermoScientific",
                  ionSource = "",
                  analyser = "Orbitrap",
                  detectorType = "Orbitrap",
                  softwareName = "MaxQuant ",
                  collisionEnergy = "",
                  dateStamp = "27 September 2018"
)

## Expression data
e1 <- exprs(davies2018ap4b1)
e2 <- exprs(davies2018ap4e1)
e3 <- exprs(davies2018wt)
  
  
## Experiment info
toName <- paste0(colnames(e1)[1:10])
colnames(e1) <- toName
pd1 <- data.frame(toName,
                 row.names=colnames(e1))  
pd1 <- new("AnnotatedDataFrame", pd1)

## Experiment info
toName <- paste0(colnames(e2)[1:10])
colnames(e2) <- toName
pd2 <- data.frame(toName,
                  row.names=colnames(e2))  
pd2 <- new("AnnotatedDataFrame", pd2)


## Experiment info
toName <- paste0(colnames(e3)[1:10])
colnames(e3) <- toName
pd3 <- data.frame(toName,
                  row.names=colnames(e3))  
pd3 <- new("AnnotatedDataFrame", pd3)


## feature data
fd1 <- rownames(e1)
fd1 <- as.data.frame(fd1)
markerdata <- as.data.frame(csv[,5])
rownames(markerdata) <- csv[,1]
fd1$markers <- "unknown"
rownames(fd1) <- csv[,1]
fd1$markers[markerdata !=""] <- unlist(lapply(markerdata[markerdata !="",1], as.character))
fd1 <- new("AnnotatedDataFrame", fd1)

## feature data
fd2 <- rownames(e2)
fd2 <- as.data.frame(fd2)
markerdata <- as.data.frame(csv[,5])
rownames(markerdata) <- csv[,1]
fd2$markers <- "unknown"
rownames(fd2) <- csv[,1]
fd2$markers[markerdata !=""] <- unlist(lapply(markerdata[markerdata !="",1], as.character))
fd2 <- new("AnnotatedDataFrame", fd2)

## feature data
fd3 <- rownames(e3)
fd3 <- as.data.frame(fd3)
markerdata <- as.data.frame(csv[,5])
rownames(markerdata) <- csv[,1]
fd3$markers <- "unknown"
rownames(fd3) <- csv[,1]
fd3$markers[markerdata !=""] <- unlist(lapply(markerdata[markerdata !="",1], as.character))
fd3 <- new("AnnotatedDataFrame", fd3)


process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("median Normalisation")),
               normalised=FALSE)

davies2018ap4b1 <- new("MSnSet",
                    exprs = e1,
                    phenoData = pd1,
                    experimentData = experiment,
                    featureData = fd1)

davies2018ap4e1 <- new("MSnSet",
                       exprs = e2,
                       phenoData = pd2,
                       experimentData = experiment,
                       featureData = fd2)

davies2018wt <- new("MSnSet",
                       exprs = e3,
                       phenoData = pd3,
                       experimentData = experiment,
                       featureData = fd3)

## Normalise
davies2018ap4b1 <- normalise(davies2018ap4b1, method = "sum")
davies2018ap4e1 <- normalise(davies2018ap4e1, method = "sum")
davies2018wt <- normalise(davies2018wt, method = "sum")

#plot2D(orre2019a431, main = "A431", method = "t-SNE")
#addLegend(orre2019a431, where = "bottomleft", ncol = 1, cex = 0.9)

## Phenodata
pData(davies2018ap4b1)$fraction <- matrix(unlist(strsplit(sampleNames(davies2018ap4b1), "[_]")), 10, 3, byrow = T)[,3]
pData(davies2018ap4b1)$replicate <- matrix(unlist(strsplit(sampleNames(davies2018ap4b1), "[_]")), 10, 3, byrow = T)[,2]

pData(davies2018ap4e1)$fraction <- matrix(unlist(strsplit(sampleNames(davies2018ap4e1), "[_]")), 10, 3, byrow = T)[,3]
pData(davies2018ap4e1)$replicate <- matrix(unlist(strsplit(sampleNames(davies2018ap4e1), "[_]")), 10, 3, byrow = T)[,2]

pData(davies2018wt)$fraction <- matrix(unlist(strsplit(sampleNames(davies2018wt), "[_]")), 10, 3, byrow = T)[,3]
pData(davies2018wt)$replicate <- matrix(unlist(strsplit(sampleNames(davies2018wt), "[_]")), 10, 3, byrow = T)[,2]

stopifnot(length(pData(davies2018ap4b1)$replicate) == ncol(e1)) # check columns and experiments match
stopifnot(length(pData(davies2018ap4e1)$replicate) == ncol(e1)) # check columns and experiments match
stopifnot(length(pData(davies2018wt)$replicate) == ncol(e1)) # check columns and experiments match

davies2018ap4b1@processingData <- process
davies2018ap4e1@processingData <- process
davies2018wt@processingData <- process
stopifnot(validObject(davies2018ap4b1))
stopifnot(validObject(davies2018ap4e1))
stopifnot(validObject(davies2018wt))

save(davies2018ap4b1, file="../../data/davies2018ap4b1.rda",
     compress = "xz", compression_level = 9)

save(davies2018ap4e1, file="../../data/davies2018ap4e1.rda",
     compress = "xz", compression_level = 9)

save(davies2018wt, file="../../data/davies2018wt.rda",
     compress = "xz", compression_level = 9)
