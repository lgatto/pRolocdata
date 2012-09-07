library(MSnbase)
xx <- read.csv("../extdata/AT_CHLORO_table_120906.csv", row.names = 1)

eset <- as.matrix(xx[, 6:9])
colnames(eset) <- c("ENV", "STR", "THY", "TOT")

pd <- new("AnnotatedDataFrame",
          data = data.frame(Fraction =
            c("envelope", "stroma", "thylakoids", "chloroplast"),
            row.names = colnames(eset)))

fd <- xx[, -(6:9)]
colnames(fd) <-
  c("Accession",
    "Percent_ENV", "Percent_STR", "Percent_THY",
    "TotalSpectralCount",
    "LOC_Training", "LOC_Test",
    "cTP_ChloroP", "TargetP_PPDB",
    "Description_TAIR",
    "MapManBin_PPDB",
    "CuratedLocation_PPDB",
    "Location_TAIR",
    "SubcellularLoc_MF120725", "SubplastidialLoc_MF120725",
    "SubsubplastidialLoc_120725",
    "TMHMM_PPDB",
    "CalcMW_PPDB",             
    "CalcPI_PPDB",
    "Length",
    "ProteomicsPub_PPDB",
    "GFP.YFP_PPDB")


fd$markers <- as.character(fd$LOC_Training)
fd$markers[fd$markers == " -"] <- "unknown"
fd$markers <- factor(fd$markers)



metaData <- data.frame(labelDescription = c(
                         "Gene Locus Accession",
                         "Percentage of occurrence in each of the ENV sub-chloroplast compartment, as calculated in Ferro et al., 2010",
                         "Percentage of occurrence in each of the STR sub-chloroplast compartment, as calculated in Ferro et al., 2010",
                         "Percentage of occurrence in each of the THY sub-chloroplast compartment, as calculated in Ferro et al., 2010",
                         "Total Spectral count, for the whole set of data",
                         "Annotated sub-plastidial localization for the training set used to set up the logistic regression model; see Ferro et al., 2010",
                         "Annotated sub-plastidial localization for the benchmark set used to test the logistic regression model; see Ferro et al., 2011",
                         "Prediction of a plastidial transit peptide by the ChloroP software (http://www.cbs.dtu.dk/services/ChloroP/)",
                         "Subcellular localization predicted by TargetP (information retrieved from the Plant Proteome Database - http://ppdb.tc.cornell.edu/)",
                         "Protein description as found in TAIR (information retrieved from PPDB)",
                         "Functional classification as found in MapManBin (information retrieved from PPDB)",
                         "Curated localisation as found in PPDB",
                         "Localization as found TAIR (information retrieved from PPDB)",
                         "Subcellular localisation MF120725 (curated annotations): chloroplast (C), cytosol, mitochondria (Mito), plasma membrane (PM), cell wall (CW), endoplasmic reticulum (ER), Nucleus, Peroxisome, Tonoplast, Unknown, vacuole", 
                         "Subplastidial localisation MF120725 (curated annotations): envelope (ENV), stroma (STR), thylakoids (THY), Not applicable (N/A), plastoglobules (PG), nucleoid, unknown",
                         "Subsubplastidial localisation 120725 (curated annotations): inner membrane (IM), stromal side, lumen, membrane (mb), Not applicable (N/A), Nucleoid, outer membrane (OM), thylakoids (thy), unknown",
                         "TMHMM prediction (PPDB)",
                         "Calculated MW (PPDB)",
                         "Calculated pI (PPDB)",
                         "Length of the protein (number of amino acids)",
                         "Proteomics publications (PPDB)",
                         "GFP/YFP publications (PPDB)",
                         "Marker proteins"))

fd <- new("AnnotatedDataFrame", data = fd, varMetadata = metaData)

ed <- new("MIAPE",
          title = "AT_CHLORO, a comprehensive chloroplast proteome database with subplastidial localization and curated information on envelope proteins.",
          url = "http://www.grenoble.prabi.fr/at_chloro/",
          abstract = "Recent advances in the proteomics field have allowed a series of high throughput experiments to be conducted on chloroplast samples, and the data are available in several public databases. However, the accurate localization of many chloroplast proteins often remains hypothetical. This is especially true for envelope proteins. We went a step further into the knowledge of the chloroplast proteome by focusing, in the same set of experiments, on the localization of proteins in the stroma, the thylakoids, and envelope membranes. LC-MS/MS-based analyses first allowed building the AT_CHLORO database (http://www.grenoble.prabi.fr/protehome/grenoble-plant-proteomics/), a comprehensive repertoire of the 1323 proteins, identified by 10,654 unique peptide sequences, present in highly purified chloroplasts and their subfractions prepared from Arabidopsis thaliana leaves. This database also provides extensive proteomics information (peptide sequences and molecular weight, chromatographic retention times, MS/MS spectra, and spectral count) for a unique chloroplast protein accurate mass and time tag database gathering identified peptides with their respective and precise analytical coordinates, molecular weight, and retention time. We assessed the partitioning of each protein in the three chloroplast compartments by using a semiquantitative proteomics approach (spectral count). These data together with an in-depth investigation of the literature were compiled to provide accurate subplastidial localization of previously known and newly identified proteins. A unique knowledge base containing extensive information on the proteins identified in envelope fractions was thus obtained, allowing new insights into this membrane system to be revealed. Altogether, the data we obtained provide unexpected information about plastidial or subplastidial localization of some proteins that were not suspected to be associated to this membrane system. The spectral counting-based strategy was further validated as the compartmentation of well known pathways (for instance, photosynthesis and amino acid, fatty acid, or glycerolipid biosynthesis) within chloroplasts could be dissected. It also allowed revisiting the compartmentation of the chloroplast metabolism and functions.",
          pubMedIds = "20061580")


at_chloro <- new("MSnSet",
                 exprs = eset,
                 phenoData = pd,
                 featureData = fd,
                 experimentData = ed)

stopifnot(rowSums(exprs(at_chloro)) == fData(at_chloro)$TotalSpectralCount)

if (validObject(at_chloro))
  save(at_chloro,
       file = "../../data/at_chloro.rda",
       compression_level = 9, compress = "xz")
