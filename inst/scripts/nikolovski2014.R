library("pRolocdata")

f <- c("../extdata/245589ST2_protein_distributions.csv",
       "../extdata/245589ST3_MarkerList_250614.csv",
       "../extdata/245589ST4_SVMLocalisation_280514.csv")


## Expression data
## getEcols(f[1], ",")
e <- 9:28
nikolovski2014 <- readMSnSet2(f[1], e-1, row.names = 1)
fData(nikolovski2014)$markers <- 
    as.character(fData(nikolovski2014)[, "Subcellular.localisation.marker..known.localisation."])
fData(nikolovski2014)$markers[fData(nikolovski2014)$markers == " "] <- "unknown"
fData(nikolovski2014)[, "Subcellular.localisation.marker..known.localisation."] <- NULL
pData(nikolovski2014) <- 
    data.frame(gradient = substr(sampleNames(nikolovski2014), 1, 1), 
               fraction = as.numeric(sub("[AB]_", "", sampleNames(nikolovski2014))),
               row.names = sampleNames(nikolovski2014))


## Adding markers
markers <- read.csv(f[2], row.names = 1, stringsAsFactors = FALSE)
i <- match(rownames(markers), featureNames(nikolovski2014))
fData(nikolovski2014)$PMID <-
    fData(nikolovski2014)$Localisation.Method <- NA
fvarLabels(nikolovski2014)[8:9] <- 
    paste(fvarLabels(nikolovski2014)[8:9], "Markers", sep = ".")
fData(nikolovski2014)[i, c("Localisation.Method.Markers", "PMID.Markers")] <-
    markers[, c("Localisation.Method", "PMID")]

## Adding localisation results
res <- read.csv(f[3], row.names = 1)[, -(1:2)]
i <- match(rownames(res), featureNames(nikolovski2014))
fData(nikolovski2014)[, names(res)] <- NA
fData(nikolovski2014)[i, names(res)] <- res

stopifnot(validObject(nikolovski2014))

save(nikolovski2014, file = "../../data/nikolovski2014.rda",
     compress = "xz", compression_level = 9)
