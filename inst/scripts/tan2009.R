mrk <- read.csv("../extdata/pr800866n_si_007.csv",
                row.names = 1,
                stringsAsFactors = FALSE)

makeTan <- function(csvfile, markers = mrk) {
  require("MSnbase")
  xx <- read.csv(csvfile, row.names = 1, stringsAsFactors = FALSE)
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
  pd <- data.frame(Fractions = c("4/5", "12/13", "19", "21"),
                   row.names = paste0("X", 114:117))
  exp <- new("MIAPE",
             lab = "Cambridge Centre for Proteomics (CCP)",
             name = "Kathryn S. Lilley",
             email = "k.s.lilley@bioc.cam.ac.uk",
             samples = list(
               species = "Drosophila melanogaster",
               tissue = "Embryos (0-16 h old)"
               ),
             title = "Mapping Organelle Proteins and Protein Complexes in Drosophila melanogaster",
             abstract = "Many proteins within eukaryotic cells are organized spatially and functionally into membrane bound organelles and complexes. A protein’s location thus provides information about its function. Here, we apply LOPIT, a mass-spectrometry based technique that simultaneously maps proteins to specific subcellular compartments, to Drosophila embryos. We determine the subcellular distribution of hundreds of proteins, and protein complexes. Our results reveal the potential of LOPIT to provide average snapshots of cells.",
             pubMedIds = "19317464",
             url = "http://www.bio.cam.ac.uk/proteomics/",
             instrumentModel = "QSTAR",
             instrumentManufacturer = "Applied Biosystems",
             ionSource = "ESI",
             analyser = "TOF",
             detectorType = "PMT")                           
  ans <- new("MSnSet",
             exprs = eset,
             experimentData = exp,
             featureData = new("AnnotatedDataFrame", data = fd),
             phenoData = new("AnnotatedDataFrame", data = pd))  
  if (validObject(ans))
    return(ans)
}

tan2009r1 <- makeTan("../extdata/pr800866n_si_004-rep1.csv")
tan2009r2 <- makeTan("../extdata/pr800866n_si_004-rep2.csv")
tan2009r3 <- makeTan("../extdata/pr800866n_si_004-rep3.csv")

save(tan2009r1, file="../../data/tan2009r1.RData",
     compress = "xz", compression_level = 9)
save(tan2009r2, file="../../data/tan2009r2.RData",
     compress = "xz", compression_level = 9)
save(tan2009r3, file="../../data/tan2009r3.RData",
     compress = "xz", compression_level = 9)


