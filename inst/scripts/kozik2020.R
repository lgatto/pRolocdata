library("MSnbase")
library("pRoloc")

csvfile <- "../../inst/extdata/Kozik_Dom_Data.csv"
csv <- read.csv(csvfile)
#getEcols(csvfile, split = ",", n = 3)

## There are 2 replicates in this dataset
Kozik_con <- readMSnSet2(file = csvfile, ecol = 8:17, skip = 0, fnames = 1)
Kozik_pra <- readMSnSet2(file = csvfile, ecol = 18:27, skip = 0, fnames = 1)
Kozik_tam <- readMSnSet2(file = csvfile, ecol = 28:37, skip = 0, fnames = 1)

## Experimental data to add
experiment <- new("MIAPE",
                  lab = "MRC LMB",
                  name = "Patrycja Kozik",
                  contact = "Patrycja Kozik",
                  email = "pkozik@mrc-lmb.cam.ac.uk <pkozik@mrc-lmb.cam.ac.uk>;",
                  samples = list(
                    species = "dendritic cells",
                    operator = "Patrycja Kozik"
                  ),
                  title = "Small Molecule Enhancers of Endosome-to-CytosolImport Augment Anti-tumor Immunity",
                  abstract = "Cross-presentation of antigens by dendritic cells (DCs) is critical for initiation of anti-tumor immune re-sponses. Yet, key steps involved in trafficking of antigens taken up by DCs remain incompletely understood.Here, we screen 700 US Food and Drug Administration (FDA)-approved drugs and identify 37 enhancers ofantigen import from endolysosomes into the cytosol. To reveal their mechanism of action, we generate pro-teomic organellar maps of control and drug-treated DCs (focusing on two compounds, prazosin and tamox-ifen). By combining organellar mapping, quantitative proteomics, and microscopy, we conclude that importenhancers undergo lysosomal trapping leading to membrane permeation and antigen release. Enhancing an-tigen import facilitates cross-presentation of soluble and cell-associated antigens. Systemic administrationof prazosin leads to reduced growth of MC38 tumors and to a synergistic effect with checkpoint immuno-therapy in a melanoma model. Thus, inefficient antigen import into the cytosol limits antigen cross-presen-tation, restraining the potency of anti-tumor immune responses and efficacy of checkpoint blockers",
                  pubMedIds = "",
                  url = "",
                  instrumentModel = "Q Exactive",
                  instrumentManufacturer = "ThermoScientific",
                  ionSource = "",
                  analyser = "Orbitrap",
                  detectorType = "Orbitrap",
                  softwareName = "MaxQuant ",
                  collisionEnergy = "",
                  dateStamp = "3 January 2018"
)

## Expression data
e <- exprs(Kozik_con)

## Experiment info
toName <- paste0(colnames(Kozik_con)[1:10])
colnames(e) <- toName
pd <- data.frame(toName,
                 row.names=colnames(e))  
pd <- new("AnnotatedDataFrame", pd)

## feature data
fd <- rownames(e)
fd <- as.data.frame(fd)
fd$markers <- "unknown"
fd$markers[csv$Marker != ""] <- csv$Marker[csv$Marker != ""]
rownames(fd) <- rownames(e)
fd <- new("AnnotatedDataFrame", fd)

process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("")),
               normalised=FALSE)

Kozik_con <- new("MSnSet",
                    exprs = e,
                    phenoData = pd,
                    experimentData = experiment,
                    featureData = fd)


## Expression data
e <- exprs(Kozik_pra)

## Experiment info
toName <- paste0(colnames(Kozik_pra)[1:10])
colnames(e) <- toName
pd <- data.frame(toName,
                 row.names=colnames(e))  
pd <- new("AnnotatedDataFrame", pd)

## feature data
fd <- rownames(e)
fd <- as.data.frame(fd)
fd$markers <- "unknown"
fd$markers[csv$Marker != ""] <- csv$Marker[csv$Marker != ""]
rownames(fd) <- rownames(e)
fd <- new("AnnotatedDataFrame", fd)

process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("")),
               normalised=FALSE)

Kozik_pra <- new("MSnSet",
                 exprs = e,
                 phenoData = pd,
                 experimentData = experiment,
                 featureData = fd)

## Expression data
e <- exprs(Kozik_tam)

## Experiment info
toName <- paste0(colnames(Kozik_tam)[1:10])
colnames(e) <- toName
pd <- data.frame(toName,
                 row.names=colnames(e))  
pd <- new("AnnotatedDataFrame", pd)

## feature data
fd <- rownames(e)
fd <- as.data.frame(fd)
fd$markers <- "unknown"
fd$markers[csv$Marker != ""] <- csv$Marker[csv$Marker != ""]
rownames(fd) <- rownames(e)
fd <- new("AnnotatedDataFrame", fd)

process <- new("MSnProcess",
               processing=c(
                 paste("Loaded on ",date(),".",sep=""),
                 paste("")),
               normalised=FALSE)

Kozik_tam <- new("MSnSet",
                 exprs = e,
                 phenoData = pd,
                 experimentData = experiment,
                 featureData = fd)


save(Kozik_con, file="../../data/Kozik_con.rda",
     compress = "xz", compression_level = 9)

save(Kozik_tam, file="../../data/Kozik_tam.rda",
     compress = "xz", compression_level = 9)

save(Kozik_pra, file="../../data/Kozik_pra.rda",
     compress = "xz", compression_level = 9)