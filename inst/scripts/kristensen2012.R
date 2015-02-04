library("MSnbase")
library("pRoloc")

i <- grepEcols("../extdata/22863883_Rep1.txt.gz", pattern = "Ratio", split = "\t")
kristensen2012r1 <- readMSnSet2("../extdata/22863883_Rep1.txt.gz", sep = "\t",
                                ecol = i, as.is = TRUE)

i <- fData(kristensen2012r1)$AC == ""
fData(kristensen2012r1)$AC[i] <- fData(kristensen2012r1)$Protein.IDs[i]
featureNames(kristensen2012r1) <- fData(kristensen2012r1)$AC 


i <- grepEcols("../extdata/22863883_Rep2.txt.gz", pattern = "Ratio", split = "\t")
kristensen2012r2 <- readMSnSet2("../extdata/22863883_Rep2.txt.gz", sep = "\t",
                                ecol = i, as.is = TRUE)
i <- fData(kristensen2012r2)$AC == ""
fData(kristensen2012r2)$AC[i] <- fData(kristensen2012r2)$Protein.IDs[i]
featureNames(kristensen2012r2) <- fData(kristensen2012r2)$AC 


i <- grepEcols("../extdata/22863883_Rep3.txt.gz", pattern = "Ratio", split = "\t")
kristensen2012r3 <- readMSnSet2("../extdata/22863883_Rep3.txt.gz", sep = "\t",
                                ecol = i, as.is = TRUE)
fvarLabels(kristensen2012r3)[2] <- "Protein.IDs"
i <- fData(kristensen2012r3)$AC == ""
fData(kristensen2012r3)$AC[i] <- fData(kristensen2012r3)$Protein.IDs[i]
featureNames(kristensen2012r3) <- fData(kristensen2012r3)$AC 

experimentData(kristensen2012r1) <-
    experimentData(kristensen2012r2) <-
        experimentData(kristensen2012r3) <-
            new("MIAPE",
                title = "A high-throughput approach for measuring temporal changes in the interactome.",
                abstract = "Interactomes are often measured using affinity purification-mass spectrometry (AP-MS) or yeast two-hybrid approaches, but these methods do not provide stoichiometric or temporal information. We combine quantitative proteomics and size-exclusion chromatography to map 291 coeluting complexes. This method allows mapping of an interactome to the same depth and accuracy as AP-MS with less work and without overexpression or tagging. The use of triplex labeling enables monitoring of interactome rearrangements.",
                name = "Kristensen AR, Gsponer J and Foster LJ",
                pubMedIds = "22863883",
                samples  = list('HeLa cells (Human)'),
                other = list(Quantification = 'SILAC',
                    Fractionation = 'Size Exclusion Chromatography'))

load("../extdata/protein-complex-markers.rda")
kristensen2012r1 <- addMarkers(kristensen2012r1, pcmrk)
kristensen2012r2 <- addMarkers(kristensen2012r2, pcmrk)
kristensen2012r3 <- addMarkers(kristensen2012r3, pcmrk)

stopifnot(validObject(kristensen2012r1))
stopifnot(validObject(kristensen2012r2))
stopifnot(validObject(kristensen2012r3))

save(kristensen2012r1, file = "../../data/kristensen2012r1.rda")
save(kristensen2012r2, file = "../../data/kristensen2012r2.rda")
save(kristensen2012r3, file = "../../data/kristensen2012r3.rda")
