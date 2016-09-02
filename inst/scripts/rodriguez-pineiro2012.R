library("readxl")
library("pRoloc")

r1 <- read_excel("../extdata/pr2010988_si_003.xls", sheet = 2, skip = 1)
rodriguez2012r1 <- readMSnSet2(r1, grep("fraction", names(r1)))
rodriguez2012r1$fraction <- as.numeric(sub("fraction ", "", sampleNames(rodriguez2012r1)))
rodriguez2012r1$rep <- 1
featureNames(rodriguez2012r1) <- fData(rodriguez2012r1)[, 1]
fData(rodriguez2012r1) <- fData(rodriguez2012r1)[, -1]
stopifnot(validObject(rodriguez2012r1))

r2 <- read_excel("../extdata/pr2010988_si_003.xls", sheet = 3, skip = 1)
rodriguez2012r2 <- readMSnSet2(r2, grep("fraction", names(r2)))
rodriguez2012r2$fraction <- as.numeric(sub("fraction ", "", sampleNames(rodriguez2012r2)))
rodriguez2012r2$rep <- 2
featureNames(rodriguez2012r2) <- fData(rodriguez2012r2)[, 1]
fData(rodriguez2012r2) <- fData(rodriguez2012r2)[, -1]
stopifnot(validObject(rodriguez2012r2))

r3 <- read_excel("../extdata/pr2010988_si_003.xls", sheet = 4, skip = 1)
rodriguez2012r3 <- readMSnSet2(r3, grep("fraction", names(r3)))
rodriguez2012r3$fraction <- as.numeric(sub("fraction ", "", sampleNames(rodriguez2012r3)))
rodriguez2012r3$rep <- 3
featureNames(rodriguez2012r3) <- fData(rodriguez2012r3)[, 1]
fData(rodriguez2012r3) <- fData(rodriguez2012r3)[, -1]
stopifnot(validObject(rodriguez2012r3))

library("pRolocdata")
data(hyperLOPIT2015)

.addMarkers <- function(obj) {
    gn <- fData(obj)$Genes
    hlgn <- sub(" .+$", "", sub("^.+GN=", "",
                                fData(hyperLOPIT2015)$protein.description))
    ii <- match(tolower(gn), tolower(hlgn))
    fData(obj)$markers <- "unknown"
    fData(obj)[!is.na(ii), "markers"] <-
        fData(hyperLOPIT2015)[na.omit(ii), "markers"]
    if (validObject(obj)) obj
}

rodriguez2012r1 <- .addMarkers(rodriguez2012r1)
rodriguez2012r2 <- .addMarkers(rodriguez2012r2)
rodriguez2012r3 <- .addMarkers(rodriguez2012r3)

save(rodriguez2012r1, file = "../../data/rodriguez2012r1.rda", compress = "xz", compression_level = 9)
save(rodriguez2012r2, file = "../../data/rodriguez2012r2.rda", compress = "xz", compression_level = 9)
save(rodriguez2012r3, file = "../../data/rodriguez2012r3.rda", compress = "xz", compression_level = 9)
