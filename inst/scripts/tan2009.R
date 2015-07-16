library("pRolocdata")
library("pRoloc")

mrk <- read.csv("../extdata/pr800866n_si_007.csv.gz",
                row.names = 2,
                stringsAsFactors = FALSE)

makeTan <- function(csvfile, markers = mrk) {
  require("MSnbase")
  xx <- read.csv(csvfile, row.names = 2, stringsAsFactors = FALSE)
  eset <- as.matrix(xx[, grep("area", names(xx))])
  colnames(eset) <- paste0("X", 114:117)
  fd <- xx[, c(1:5,10)]
  names(fd)[3] <- "No.peptide.IDs"
  names(fd)[5] <- "No.peptide.quantified"
  names(fd)[6] <- "PLSDA"
  fd$PLSDA[fd$PLSDA == ""] <- "unknown"
  fd$PLSDA[grep("asma membrane", fd$PLSDA)] <- "PM" ## r2/3 and 'Plasma membrane'/'plasma membrane'
  fd$PLSDA <- factor(fd$PLSDA)
  fd$markers <- "unknown"
  getNA <- is.na(match(rownames(mrk), rownames(xx)))
  if (any(getNA)) {
    foo <- na.omit(match(rownames(mrk), rownames(xx)))
    fd$markers[foo] <- mrk$Organelle[-na.action(foo)]
  } else {
    fd$markers[match(rownames(mrk), rownames(xx))] <- mrk$Organelle
  }
  fd$markers[fd$markers == "mito"] <- "mitochondrion"
  fd$markers <- factor(fd$markers)
  if (any(names(xx)=="pd.2013")) {
    fd$pd.2013 <- xx[,14]
  }
  if (any(names(xx)=="pd.markers")) {
    fd$pd.markers <- xx[,15]
  }
  pd <- data.frame(Fractions = c("4/5", "12/13", "19", "21"),
                   row.names = paste0("X", 114:117))
  exp <- new("MIAPE",
             lab = "Cambridge Centre for Proteomics (CCP)",
             name = "Kathryn S. Lilley",
             email = "k.s.lilley@bioc.cam.ac.uk",
             samples = list(
               species = "Drosophila melanogaster",
               tissue = "Embryos"
               ),
             title = "Mapping Organelle Proteins and Protein Complexes in Drosophila melanogaster",
             abstract = "Many proteins within eukaryotic cells are organized spatially and functionally into membrane bound organelles and complexes. A protein's location thus provides information about its function. Here, we apply LOPIT, a mass-spectrometry based technique that simultaneously maps proteins to specific subcellular compartments, to Drosophila embryos. We determine the subcellular distribution of hundreds of proteins, and protein complexes. Our results reveal the potential of LOPIT to provide average snapshots of cells.",
             pubMedIds = "19317464",
             url = "http://www.bio.cam.ac.uk/proteomics/",
             instrumentModel = "QSTAR",
             instrumentManufacturer = "Applied Biosystems",
             ionSource = "ESI",
             analyser = "TOF",
             detectorType = "PMT",
             other = list(
                 MS = "iTRAQ4",
                 spatexp = "LOPIT",
                 markers.fcol = "markers",
                 prediction.fcol = "PLSDA"))
  fd.Ann <- new("AnnotatedDataFrame", data = fd)
  if (any(names(xx)=="pd.markers")) {
  	fd.Ann@varMetadata[,1] <- c("CG number",
                            "FlyBase symbol/name",
                            "Peptides",
                            "Mascot score",
                            "Number of peptides quantified",
                            "Protein localisation assigned by PLSDA experiment described in Tan et al 2009",
                            "Original protein markers used in PLSDA analysis in Tan et al 2009",
			     "PhenoDisco output as described in Breckels et al (2013) Journal of Proteomics. Accepted February 2013",
                            "Updated protein markers (original markers see column 'markers' plus new protein markers found by phenoDisco and verified by UniprotKB/literature as described in Breckels et al)") 
  } else {
  fd.Ann@varMetadata[,1] <- c("CG number",
                            "FlyBase symbol/name",
                            "Peptides",
                            "Mascot score",
                            "Number of peptides quantified",
                            "Protein localisation assigned by PLSDA experiment described in Tan et al 2009",
                            "Original protein markers used in PLSDA analysis in Tan et al 2009")
  }
  ans <- new("MSnSet",
             exprs = eset,
             experimentData = exp,
             featureData = fd.Ann,
             phenoData = new("AnnotatedDataFrame", data = pd))  
  if (validObject(ans))
    return(ans)
}

tan2009r1 <- makeTan("../extdata/pr800866n_si_004-rep1.csv.gz")
tan2009r2 <- makeTan("../extdata/pr800866n_si_004-rep2.csv.gz")
tan2009r3 <- makeTan("../extdata/pr800866n_si_004-rep3.csv.gz")

## Function for adding Uniprot IDs to the Tan datasets
source("addTanIds.R")
tan2009r1 <- addTanIds(tan2009r1)
tan2009r2 <- addTanIds(tan2009r2)
tan2009r3 <- addTanIds(tan2009r3)

fData(tan2009r1) <- fData(tan2009r1)[, c(1:2, 10:13, 3:9)]
fData(tan2009r2) <- fData(tan2009r2)[, c(1:2, 8:11, 3:7)]
fData(tan2009r3) <- fData(tan2009r3)[, c(1:2, 8:11, 3:7)]

fvarMetadata(tan2009r1)$labelDescription[c(3, 5)] <- rep("Uniprot Accession Number converted from FlyBase using CG numbers as input (see README for details)", 2) 
fvarMetadata(tan2009r2)$labelDescription[c(3, 5)] <- rep("Uniprot Accession Number converted from FlyBase using CG numbers as input (see README for details)", 2) 
fvarMetadata(tan2009r3)$labelDescription[c(3, 5)] <- rep("Uniprot Accession Number converted from FlyBase using CG numbers as input (see README for details)", 2) 
fvarMetadata(tan2009r1)$labelDescription[c(4, 6)] <- rep("Uniprot Entry Name converted from FlyBase using CG numbers as input (see README for details)", 2) 
fvarMetadata(tan2009r2)$labelDescription[c(4, 6)] <- rep("Uniprot Entry Name converted from FlyBase using CG numbers as input (see README for details)", 2) 
fvarMetadata(tan2009r3)$labelDescription[c(4, 6)] <- rep("Uniprot Entry Name converted from FlyBase using CG numbers as input (see README for details)", 2) 

tan2009r1@experimentData@other$markers.fcol <- "pd.markers"
tan2009r1@experimentData@other$prediction.fcol <- "pd.2013"

stopifnot(pRolocdata:::valid.pRolocmetadata(pRolocmetadata(tan2009r1)))
stopifnot(pRolocdata:::valid.pRolocmetadata(pRolocmetadata(tan2009r2)))
stopifnot(pRolocdata:::valid.pRolocmetadata(pRolocmetadata(tan2009r3)))

## Make Uniprot accession number the featureNames as it is the most stable ID
fData(tan2009r1)$FBgn <- featureNames(tan2009r1)
fData(tan2009r2)$FBgn <- featureNames(tan2009r2)
fData(tan2009r3)$FBgn <- featureNames(tan2009r3)
featureNames(tan2009r1) <- fData(tan2009r1)$AccessionNo
featureNames(tan2009r2) <- fData(tan2009r2)$AccessionNo
featureNames(tan2009r3) <- fData(tan2009r3)$AccessionNo
fData(tan2009r1) <- fData(tan2009r1)[c(length(fvarLabels(tan2009r1)), 1:(length(fvarLabels(tan2009r1))-1))]
fData(tan2009r2) <- fData(tan2009r2)[c(length(fvarLabels(tan2009r2)), 1:(length(fvarLabels(tan2009r2))-1))]
fData(tan2009r3) <- fData(tan2009r3)[c(length(fvarLabels(tan2009r3)), 1:(length(fvarLabels(tan2009r3))-1))]

## Add updated marker lists
load("../extdata/markersTan.rda")
fData(tan2009r1)$markers.orig <- fData(tan2009r1)$markers
fData(tan2009r1)$markers <- NULL
tan2009r1 <- addMarkers(tan2009r1, mrk, verbose = FALSE)
# remove any marker mismatches
# c("P38979", "O76927")
ind <- c(195, 339, 654)
fData(tan2009r1)$markers[ind] <- rep("unknown", length(ind))

fData(tan2009r2)$markers.orig <- fData(tan2009r2)$markers
fData(tan2009r2)$markers <- NULL
tan2009r2 <- addMarkers(tan2009r2, mrk, verbose = FALSE)
# remove the below as only 4 members in each that do not cluster
nuc <- which(fData(tan2009r2)$markers == "Nucleus")
cyt <- which(fData(tan2009r2)$markers == "Cytoskeleton")
per <- which(fData(tan2009r2)$markers == "Peroxisome")
fData(tan2009r2)$markers[c(cyt, nuc, per)] <- rep("unknown", length(c(cyt,nuc, per)))

fData(tan2009r3)$markers.orig <- fData(tan2009r3)$markers
fData(tan2009r3)$markers <- NULL
tan2009r3 <- addMarkers(tan2009r3, mrk, verbose = FALSE)
ga <- which(fData(tan2009r3)$markers == "Golgi")
cyt <- which(fData(tan2009r3)$markers == "Cytoskeleton")
per <- which(fData(tan2009r3)$markers == "Peroxisome")
fData(tan2009r3)$markers[c(cyt, ga, per)] <- rep("unknown", length(c(cyt, ga, per)))
ind <- c(132, 650)
fData(tan2009r3)$markers[ind] <- rep("unknown", length(ind))

## Add transfer learning results from Breckels et al 2015
load("../extdata/tl-res/tan-tl.rda")
tan2009r1 <- minMarkers(tan2009r1, 13)
ind <- which(fvarLabels(tan2009r1) == "markers13")
fvarLabels(tan2009r1)[ind] <- "markers.tl"
#experimentData(tan2009r1)@other$knntl <- tan.tl$knntl
#experimentData(tan2009r1)@other$svmtl <- tan.tl$svmtl
# Now added to GO CC equivalent set

## Update fvarMetaData slots
fvarMetadata(tan2009r1)[c(1, 14:16),1]  <- c("FBgn number from FlyBase",
                                             "Original markers from Tan (2009) et al",
                                             "Updated and curated protein markers, same as pd.markers",
                                             "Markers used in transfer learning paper by Breckels et al (2015)") 

stopifnot(pRolocdata:::valid.pRolocmetadata(pRolocmetadata(tan2009r1)))
stopifnot(pRolocdata:::valid.pRolocmetadata(pRolocmetadata(tan2009r2)))
stopifnot(pRolocdata:::valid.pRolocmetadata(pRolocmetadata(tan2009r3)))

save(tan2009r1, file="../../data/tan2009r1.RData",
     compress = "xz", compression_level = 9)
save(tan2009r2, file="../../data/tan2009r2.RData",
     compress = "xz", compression_level = 9)
save(tan2009r3, file="../../data/tan2009r3.RData",
     compress = "xz", compression_level = 9)

