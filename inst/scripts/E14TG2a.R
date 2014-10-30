library("MSnbase")

makeAndy <- function(filename, date, stringency) {
  csv <- read.csv(filename, row.names=1, header = TRUE, stringsAsFactors=FALSE)
  Uniprot.ID <- substr(rownames(csv), start=4, stop=9)
  Accession.No <- substr(rownames(csv), start=11, stop=nchar(rownames(csv)))
  rownames(csv) <- Accession.No
  ind <- grep("n1", colnames(csv))
  .exprs <- csv[,ind]
  ind <- c(grep("GO", colnames(csv)), grep("Train", colnames(csv)))
  .fData <- cbind(Uniprot.ID, csv[,c(1:3,ind)])
  null <- which(.fData$GOannotation=="null")
  .fData$GOannotation[null] <- rep("unknown", length(null))
  na <- which(.fData$Training=="")
  .fData$Training[na] <- rep("unknown", length(na))
  ind <- which(colnames(.fData) == "Training")
  colnames(.fData)[ind] <- "markers"
  .fData <- new("AnnotatedDataFrame", .fData)
  .fData@varMetadata[,1] <- c("UniProtKB accession number", 
                              "Full protein name", 
                              "Peptides", 
                              "Peptide spectrum match", 
                              "Localisation inferred from GO: Andy' output from his quickGO program",
                              "Andy's own curated training set")
  .pData <- new("AnnotatedDataFrame",
                data.frame(Fraction.information=c(rep("High Density", 7), 
                                                     "Soluble/cytosolic"), 
                              row.names=colnames(.exprs)))

  .experiment <- new("MIAPE",
                    lab = "Cambridge Centre for Proteomics (CCP)",
                    name = "Andy Christoforou",
                    contact = "Kathryn S. Lilley",
                    email = "k.s.lilley@bioc.cam.ac.uk",
                    samples = list(
                      species = "Mouse",
                      tissue = "E14TG2a embryonic stem cells",
                      operator = "Andy Christoforou"
                     ),
                    title="",
                    abstract = "",
                    pubMedIds = "",
                    url = "",
                    instrumentModel = "LTQ Orbitrap Velos",
                    instrumentManufacturer = "ThermoScientific",
                    ionSource = "New Objective PicoView nano-electrospray",
                    analyser = "Orbitrap",
                    detectorType = "Orbitrap",
                    softwareName = "Mascot Search Engine",
                    isolationWidth = 0.5,
                    collisionEnergy = "45% CE HCD (two-stepped, 10% width)",
                    other = list(
                      experimentType = "Experiment Type: Standard LOPIT experimental design on E14TG2a embryonic stem cells. Sample: E14TG2a mouse pluripotent embryonic stem cells cultured under conditions favouring self-renewal (serum+LIF). Label: iTRAQ 8-plex. Instrument: LTQ Orbitrap Velos. Scan Mode for Identification: MS2-HCD. Scan Mode for Quantitation: MS2-HCD. Scan Setup: Nth order double play, Top 10 HCD. MS1 Scan: FTMS, resolution = 30000, scan range m/z 380 - 1600. MS2 Scan: FTMS, resolution = 7500, scan range m/z 100 - 2000. Precursor Ion Selection Window: 0.5 Da on 'stringent' and 1.2 Da on 'relaxed' setting (see stringencySetting slot). Collision Energy: 45% CE HCD (two-stepped, 10% width)",
                      searchParameters = "Search Engine: Mascot. Search Database: UniProt Mouse. Fixed Modifications: iTRAQ 8-plex (N-term), iTRAQ-8plex (K), Methylthio (C). Variable Modifications: iTRAQ 8-plex (Y), Oxidation (M). Enzyme: Trypsin. Max. Missed Cleavages: 2. Decoy Type: None (see Percolator parameters below). Peptide Charge: 5+. Peptide Tolerance: +/- 10 ppm. MS/MS Tolerance: +/- 0.2 Da. Instrument: ESI-ORBITRAP-HCD",
                      postProcessing = "Unquantifiable spectra (E-value for PSM > 0.05, non-unique sequences, very low ion counts, >2 zero value reporter ions ) removed
                                        Spectra filtered based on several criteria (precursor relative signal, position of switch relative to peak apex, reporter ion intensity) to pick a single 'peptidotypic' spectrum per peptide
                                        Peptides were merged into proteins by intensity weighted mean",
                      percolatorParameters = "Percolator Parameters: Percolator Version Used: None. Percolator was misbehaving when this data was processed so Mascot E-values were used to benchmark ID significance rather than PEPs. Any E-values less than 0.05 were accepted for quantitation. The E-value is less stringent than PEP but my general impression is this has made little difference to the overall quality of data. Although no decoy searches were run for Percolating they were performed to check that the FDR was tolerable at the default Mascot p-value of 0.05.                      
                                        Additional filtering was applied to the data after percolating to refine protein inference within iSPY. Swissprot accessions were given precedence over trEMBL accessions, and isoforms from the same UniProt accession were collapsed together. These assumptions substantially reduce the redundancy of the database, allowing more 'unique' peptides to be taken forward for quantitation.",
                      stringencySetting = stringency,
                      MS = "iTRAQ8",
                      spatexp = "LOPIT",
                      markers.fcol = "markers",
                      prediction.fcol = NA
                      ),
                    dateStamp = date
  )               

  .process <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".",sep=""),
                  paste("Normalised to sum of intensities.")),
                normalised=TRUE,
                files=filename)

  obj <- new("MSnSet",
                exprs = .exprs,
                phenoData = .pData,
                experimentData = .experiment,
                featureData = .fData)
  obj@processingData <- .process
  if (validObject(obj))
    return (obj)
}

f1s <- "../extdata/E14TG2a_2011/E14TG2a_2011_12_stringentsettings/E14TG2a_2011_12_proteins.csv"
f2s <- "../extdata/E14TG2a_2011/E14TG2a_2012_01_stringentsettings/E14TG2a_2012_01_proteins.csv"
f1r <- "../extdata/E14TG2a_2011/E14TG2a_2011_12_relaxedsettings/E14TG2a_2011_12_1p2_proteins.csv"

date1 <- "December 2011"
date2 <- "January 2012"
H <- "High (stringent)"
L <- "Low (relaxed)"

E14TG2aS1 <- makeAndy(f1s, date1, H)
E14TG2aS2 <- makeAndy(f2s, date2, H)
E14TG2aR <- makeAndy(f1r, date1, L)

## Make uniprot accession number featureNames
foo <- function(data) {
  fData(data)$UniprotName <- featureNames(data) 
  featureNames(data) <- fData(data)$Uniprot.ID
  fData(data) <- fData(data)[, c(1, 7, 2:6)]
  fData(data)$markers.orig <- fData(data)$markers
  fData(data)$markers <- NULL
  return(data)
}

E14TG2aS1 <- foo(E14TG2aS1)
E14TG2aS2 <- foo(E14TG2aS2)
E14TG2aR <- foo(E14TG2aR)

## Add updated marker list
load("../extdata/markersE14.rda")

## --E14TG2aS1
E14TG2aS1 <- addMarkers(E14TG2aS1, markers = mrk, verbose = FALSE)
E14TG2aS1 <- minMarkers(E14TG2aS1, 6)
fData(E14TG2aS1)$markers <- fData(E14TG2aS1)$markers6
fData(E14TG2aS1)$markers6 <- NULL
## Remove annotation mismatches
#torm <- c("O55143", "Q9D1B9", "Q9CPX7", "P52503", "Q8C2E4", "Q8R0G7", "Q9CXW2", 
#          "Q99N87", "P19096", "Q8VDF2", "Q62315", "Q9JKR6", "Q6P5E4", "Q9CY27")
torm <- c(103, 946, 366, 545, 788, 821, 939, 172, 212, 346, 689, 97, 306, 369)
fData(E14TG2aS1)$markers[torm] <- rep("unknown", length(torm))

## --E14TG2aS2
E14TG2aS2 <- addMarkers(E14TG2aS2, mrk, verbose = FALSE)
E14TG2aS2 <- minMarkers(E14TG2aS2, 6)
fData(E14TG2aS2)$markers <- fData(E14TG2aS2)$markers6
fData(E14TG2aS2)$markers6 <- NULL
# c("Q9CXW2", "Q99N87", "P62858", "P52503", "Q8R0G7")
ind <- c(724,  763,  873,  958, 1042)
fData(E14TG2aS2)$markers[ind] <- rep("unknown", length(ind))

## --E14TG2aR
E14TG2aR <- addMarkers(E14TG2aR, mrk, verbose = FALSE)
E14TG2aR <- minMarkers(E14TG2aR, 6)
fData(E14TG2aR)$markers <- fData(E14TG2aR)$markers6
fData(E14TG2aR)$markers6 <- NULL
# c("P06745", "P52503", "Q99N87", "Q3UMR5", "Q8R0G7", "Q9CXW2", "Q8VDF2", "Q62315")
ind<- c(203,  810, 1029, 1402, 1643, 1799, 248, 1448)
fData(E14TG2aR)$markers[ind] <- rep("unknown", length(ind))

stopifnot(pRolocdata:::valid.pRolocmetadata(pRolocmetadata(E14TG2aR)))
stopifnot(pRolocdata:::valid.pRolocmetadata(pRolocmetadata(E14TG2aS1)))
stopifnot(pRolocdata:::valid.pRolocmetadata(pRolocmetadata(E14TG2aS2)))

if (validObject(E14TG2aR))
  save(E14TG2aR,file="../../data/E14TG2aR.RData", 
       compress = "xz", compression_level = 9)
if (validObject(E14TG2aS1))
save(E14TG2aS1,file="../../data/E14TG2aS1.RData", 
     compress = "xz", compression_level = 9)
if (validObject(E14TG2aS2)) 
save(E14TG2aS2,file="../../data/E14TG2aS2.RData", 
     compress = "xz", compression_level = 9)