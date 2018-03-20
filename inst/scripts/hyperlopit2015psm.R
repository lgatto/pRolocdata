library("MSnbase")
library("pRoloc")
library("pRolocdata")
library("dplyr")

.experiment <- new("MIAPE",
                   lab = "Cambridge Centre for Proteomics (CCP)",
                   name = "Kathryn S. Lilley",
                   contact = "Kathryn S. Lilley",
                   email = "k.s.lilley@bioc.cam.ac.uk",
                   samples = list(
                       species = "Mouse",
                       tissue = "E14TG2a embryonic stem cells",
                       operator = c("Andy Christoforou",
                                    "Claire M. Mulvey")
                   ),
                   title = "A draft map of the mouse pluripotent stem cell spatial proteome",
                   abstract = "Knowledge of the subcellular distribution of proteins is vital for understanding cellular mechanisms. Capturing the subcellular proteome in a single experiment has proven challenging, with studies focusing on specific compartments or assigning proteins to subcellular niches with low resolution and/or accuracy. Here we introduce hyperLOPIT, a method that couples extensive fractionation, quantitative high-resolution accurate mass spectrometry with multivariate data analysis. We apply hyperLOPIT to a pluripotent stem cell population whose subcellular proteome has not been extensively studied. We provide localization data on over 5,000 proteins with unprecedented spatial resolution to reveal the organization of organelles, sub-organellar compartments, protein complexes, functional networks and steady-state dynamics of proteins and unexpected subcellular locations. The method paves the way for characterizing the impact of post-transcriptional and post-translational modification on protein location and studies involving proteome-level locational changes on cellular perturbation. An interactive open-source resource is presented that enables exploration of these data.",
                   pubMedIds = "26754106",
                   url = "",
                   instrumentModel = "Orbitrap Fusion Tribrid",
                   instrumentManufacturer = "ThermoScientific",
                   ionSource = "ESI",
                   analyser = "Orbitrap",
                   detectorType = "Orbitrap",
                   softwareName = "Mascot Search Engine",
                   collisionEnergy = "")

e <- c(17, 19, 22, 25, 28, 31, 34, 37, 40, 43)

data(hyperLOPIT2015)
m <- markerMSnSet(hyperLOPIT2015)
mrk <- as.character(fData(m)$markers)
names(mrk) <- featureNames(m)

keepQuanUsage <- function(x) {
    sel <- fData(x)$Quan.Usage == "Used"
    x <- x[sel, ]
    x <- MSnbase:::logging(x, "Keep only PSMs that have Quan.Usage as Used.")
    x
}

fixSampleNames <- function(x) {
    sampleNames(x) <- sub("_", "", sampleNames(x))
    x
}


addPdata <- function(x, y) {
    i <- match(sampleNames(x), sampleNames(y))
    pData(x) <- pData(y)[i, ]
    x
}

keepTMT6plexMod <- function(x) {
    ## Note that this function should be redundant with keepQuanUsage,
    ## as PSMs without any TMT modification shouldn't be used for
    ## quantitation
    hastmt <- grepl("TMT6plex", fData(x)$Modification)
    x[hastmt, ]
}

addMarkersToPSMs <- function(x, mrk) {
    fData(x)$markers <- "unknown"
    i <- match(fData(x)$Protein.Group.Accession, names(mrk))
    fData(x)$markers[!is.na(i)] <- mrk[na.omit(i)]
    x
}


data(hyperLOPIT2015ms3r1)
hyperLOPIT2015ms3r1psm <-
    readMSnSet2("../extdata/Regrouped_MS3_rep1_full.csv.gz", e) %>%
    filterNA %>%
    keepQuanUsage %>%
    fixSampleNames %>%
    addPdata(hyperLOPIT2015ms3r1) %>%
    keepTMT6plexMod %>%
    droplevels %>%
    addMarkersToPSMs(mrk) %>%
    normalise(method = "sum")

hyperLOPIT2015ms3r1psm@experimentData <- .experiment

save(hyperLOPIT2015ms3r1psm,
     file = "../../data/hyperLOPIT2015ms3r1psm.rda",
     compress = "xz", compression_level = 9)


data(hyperLOPIT2015ms3r2)
hyperLOPIT2015ms3r2psm <-
    readMSnSet2("../extdata/Regrouped_MS3_rep2_inj2_full.csv.gz", e) %>%
    filterNA %>%
    keepQuanUsage %>%
    fixSampleNames %>%
    addPdata(hyperLOPIT2015ms3r2) %>%
    keepTMT6plexMod %>%
    droplevels %>%
    addMarkersToPSMs(mrk) %>%
    normalise(method = "sum")

hyperLOPIT2015ms3r2psm@experimentData <- .experiment

save(hyperLOPIT2015ms3r2psm,
     file = "../../data/hyperLOPIT2015ms3r2psm.rda",
     compress = "xz", compression_level = 9)

data(hyperLOPIT2015ms3r1)
hyperLOPIT2015ms2psm <-
    readMSnSet2("../extdata/Regrouped_MS2_full.csv.gz", e) %>%
    filterNA %>%
    keepQuanUsage %>%
    fixSampleNames %>%
    addPdata(hyperLOPIT2015ms3r1) %>%
    keepTMT6plexMod %>%
    droplevels %>%
    addMarkersToPSMs(mrk) %>%
    normalise(method = "sum")

hyperLOPIT2015ms2psm@experimentData <- .experiment
hyperLOPIT2015ms2psm@experimentData@instrumentModel <- "Q Exactive"

save(hyperLOPIT2015ms2psm,
     file = "../../data/hyperLOPIT2015ms2psm.rda",
     compress = "xz", compression_level = 9)
