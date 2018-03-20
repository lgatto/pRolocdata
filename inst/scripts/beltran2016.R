library("MSnbase")
library("pRoloc")
library("readxl")

fm <- "../extdata/1-s2.0-S2405471216302897-mmc4.xlsx"
m <- read_xlsx(fm, sheet = 1, skip = 1)
beltran2016markers <- m[[3]]
names(beltran2016markers) <- m[[1]]

csvfls <- dir("../extdata/", pattern = "OrganelleProfiles\\.TMT",
           full.names = TRUE)
nms <- sub("\\.h", "", sub("\\.csv\\.gz", "",
                           sub("^.+TMT\\.", "beltran2016", csvfls)))
rdafls <- file.path("../../data", paste0(nms, ".rda"))

exp <- new("MIAPE",
           title = "A Portrait of the Human Organelle Proteome In Space and Time during Cytomegalovirus Infection",
           abstract = "The organelles within a eukaryotic host are manipulated by viruses to support successful virus replication and spread of infection, yet the global impact of viral infection on host organelles is poorly understood. Integrating microscopy, subcellular fractionation, mass spectrometry, and functional analyses, we conducted a cell-wide study of organelles in primary fibroblasts throughout the time course of human cytomegalovirus (HCMV) infection. We used label-free and isobaric-labeling proteomics to characterize nearly 4,000 host and 100 viral proteins, then classified their specific subcellular locations over time using machine learning. We observed a global reorganization of proteins across the secretory pathway, plasma membrane, and mitochondria, including reorganization and processing of lysosomal proteins into distinct subpopulations and translocations of individual proteins between organelles at specific time points. We also demonstrate that MYO18A, an unconventional myosin that translocates from the plasma membrane to the viral assembly complex, is necessary for efficient HCMV replication. This study provides a comprehensive resource for understanding host and virus biology during HCMV pathogenesis.",
           pubMedIds = "27641956")

for (i in seq_along(csvfls)) {
    x <- readMSnSet2(csvfls[i], ecol = 2:7)
    x@experimentData <- exp
    featureNames(x) <- fData(x)[, 1]
    x <- addMarkers(x, beltran2016markers, verbose = FALSE)
    fData(x) <- fData(x)[, -1, drop = FALSE]
    stopifnot(validObject(x))
    assign(nms[i], x)
    save(list = nms[i], file = rdafls[i],
         compress = "xz", compression_level = 9)
}
