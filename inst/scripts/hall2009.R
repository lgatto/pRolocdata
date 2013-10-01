library("MSnbase")

plsda <- read.csv("../extdata/Hall2009_Table_S5.csv", skip = 1, row.names = 1, stringsAsFactor = FALSE)
markers <- read.csv("../extdata/Hall2009_Table_S4.csv", skip = 1, row.names = 1, stringsAsFactor = FALSE)

xx <- read.csv("../extdata/Hall2009_Table_S1.csv", row.names = 1)

e <- as.matrix(xx[, c(6:13, 18:25)])
fd <- xx[, c(1:5, 14:17)] 
fd$PLSDA <- fd$markers <- "unknown"

fd$markers[match(rownames(markers), rownames(fd))] <- markers[, "Known.localization"]
fd$PLSDA[match(rownames(plsda), rownames(fd))] <- plsda[, "PLS.DA.prediction"]

pd <- data.frame(
    fractions = rep(c('1', '4', '13+14', '21',
        '1', '9+10', '16', '18'), 2),
    replicate = rep(1:2, each = 8),
    row.names = colnames(e))

.experiment <- new("MIAPE",
                   lab = "Cambridge Centre for Proteomics (CCP)",
                   name = "Kathryn S. Lilley",
                   email = "k.s.lilley@bioc.cam.ac.uk",
                   samples = list(
                     species = "Gallus gallus",
                     tissue = "pre-B cell line DT40"),
                   title="The Organelle Proteome of the DT40 Lymphocyte Cell Line.",
                   abstract = "A major challenge in eukaryotic cell biology is to understand the roles of individual proteins and the subcellular compartments in which they reside. Here, we use the localization of organelle proteins by isotope tagging technique to complete the first proteomic analysis of the major organelles of the DT40 lymphocyte cell line. This cell line is emerging as an important research tool because of the ease with which gene knockouts can be generated. We identify 1090 proteins through the analysis of preparations enriched for integral membrane or soluble and peripherally associated proteins and localize 223 proteins to the endoplasmic reticulum, Golgi, lysosome, mitochondrion, or plasma membrane by matching their density gradient distributions to those of known organelle residents. A striking finding is that within the secretory and endocytic pathway a high proportion of proteins are not uniquely localized to a single organelle, emphasizing the dynamic steady-state nature of intracellular compartments in eukaryotic cells.",
                   pubMedIds = "19181659",
                   url = "http://www.bio.cam.ac.uk/proteomics/",
                   instrumentModel = "QSTAR XL",
                   instrumentManufacturer = "Applied Biosystems",
                   ionSource = "ESI",
                   analyser = "TOF",
                   detectorType = "PMT")

.process <- new("MSnProcess",
                processing=c(
                    paste("Loaded on ",date(),".",sep=""),
                    paste("Normalised to sum of intensities.")),
                normalised=TRUE,
                files = "../extdata/Hall2009_Table_S1.csv")



hall2009 <- new("MSnSet", exprs = e,
                featureData = new("AnnotatedDataFrame", data = fd),
                phenoData = new("AnnotatedDataFrame", data = pd),
                experimentData = .experiment)
hall2009@processingData <- .process


stopifnot(validObject(hall2009))

save(hall2009, file = "../../data/hall2009.rda",
     compress = "xz", compression_level = 9)
