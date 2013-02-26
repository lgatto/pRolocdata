library("MSnbase")

###########################################################
## From Andy's protein file "ispy_results_proteins_stringent.csv"
## Peptides were merged into proteins by intensity weighted mean
andy <- read.csv("../extdata/andy2011_ispy_results_proteins_stringent.csv",
                 row.names=1, sep="")
.exprs <- as.matrix(andy[,5:12])
.fData <- andy[,c(1:4,14,13)]
names(.fData)[5] <- "pd.2013"
.fData <- new("AnnotatedDataFrame",data=data.frame(as.matrix(.fData)))
.fData@varMetadata[,1] <- c("UniProtKB accession number", 
                            "Full protein name", 
                            "Peptides", 
                            "Peptide spectrum match", 
                            "PhenoDisco output as described in Breckels et al (2013) Journal of Proteomics. Accepted February 2013", 
                            "Protein markers")
.pData <- new("AnnotatedDataFrame",
              data=data.frame(Fraction.information=c(rep("High Density", 7), 
                                                     "Soluble/cytosolic"), 
                              row.names=colnames(.exprs)))

.experiment <- new("MIAPE",
                   lab = "Cambridge Centre for Proteomics (CCP)",
                   name = "Andy Christoforou",
                   contact = "Kathryn S. Lilley",
                   email = "k.s.lilley@bioc.cam.ac.uk",
                   samples = list(
                     species = "Human",
                     tissue = "Embryonic Kidney Fibroblast Cells (HEK293T)",
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
                     experimentType = "Standard LOPIT experimental design on HEK293T. Label: iTRAQ 8-plex. Instrument: LTQ Orbitrap Velos. Scan Mode for Identification: MS2-HCD. Scan Mode for Quantitation: MS2-HCD. Scan Setup: Nth order double play, Top 10 HCD. MS1 Scan: FTMS, resolution = 30000, scan range m/z 380 - 1600. MS2 Scan: FTMS, resolution = 7500, scan range m/z 100 - 2000. Precursor Ion Selection Window: 0.5 Da (very stringent). Collision Energy: 45% CE HCD (two-stepped, 10% width)",
                     searchParameters = "Search Engine: Mascot. Search Database: UniProt Human. Fixed Modifications: iTRAQ 8-plex (N-term), iTRAQ-8plex (K), MMTS (C). Variable Modifications: iTRAQ 8-plex (Y), Oxidation (M). Enzyme: Trypsin. Max. Missed Cleavages: 2. Decoy Type: None (see Percolator parameters below). Peptide Charge: 5+. Peptide Tolerance: ± 25 ppm. MS/MS Tolerance: ± 0.2 Da. Instrument: ESI-ORBITRAP-HCD",
                     postProcessing = "Maximum Parent Proteins: 1. Unquantifiable spectra (E-value for PSM > 0.05, non-unique sequences, very low ion counts, >2 zero value reporter ions ) removed. Spectra filtered based on several criteria (precursor relative signal, position of switch relative to peak apex, reporter ion intensity) to pick a single peptidotypic spectrum per peptide. Peptides were merged into proteins by intensity weighted mean"
                     ),
                   dateStamp = "2011-07"
                   )               

.process <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".",sep=""),
                  paste("Normalised to sum of intensities.")),
                 normalised=TRUE,
                files="andy2011_ispy_results_proteins_stringent.csv")

andy2011 <- new("MSnSet",
                   exprs = .exprs,
                   phenoData = .pData,
                   experimentData = .experiment,
                   featureData = .fData)
andy2011@processingData <- .process

if (validObject(andy2011))
  save(andy2011,file="../../data/andy2011.RData",
       compress = "xz", compression_level = 9)

