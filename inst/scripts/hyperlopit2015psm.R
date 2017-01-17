library("MSnbase")
library("pRoloc")
library("pRolocdata")
library("dplyr")

e <- c(17, 19, 22, 25, 28, 31, 34, 37, 40, 43)

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
    
save(hyperLOPIT2015ms3r1psm,
     file = "../extdata/hyperLOPIT2015ms3r1psm.rda",
     compress = "xz", compression_level = 9)


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
    
save(hyperLOPIT2015ms3r2psm,
     file = "../extdata/hyperLOPIT2015ms3r2psm.rda",
     compress = "xz", compression_level = 9)

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
    
save(hyperLOPIT2015ms2psm,
     file = "../extdata/hyperLOPIT2015ms2psm.rda",
     compress = "xz", compression_level = 9)
