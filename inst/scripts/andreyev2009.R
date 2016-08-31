library("MSnbase")
library("pRoloc")
library("readr")

xx <- read_csv("../../inst/extdata/Andreyev2009_Supplemental_Table_S1-L-C.csv.gz",
               skip = 1)
xx <- xx[, !is.na(names(xx))] ## get rid of empty columns
names(xx)[2] <- "Accession"
    
andreyev2009 <- readMSnSet2(as.data.frame(xx), ecol = 7:78)
andreyev2009$fractions <- sub("-.+$", "", sampleNames(andreyev2009))
andreyev2009$condition <- ifelse(grepl("-C", sampleNames(andreyev2009)),
                                 "Resting", "Endotoxin-activated")
andreyev2009$Var1 <- ifelse(grepl("a$", sampleNames(andreyev2009)), "A", "B")
andreyev2009$Rep <- 0
andreyev2009$Rep[grepl("1", sampleNames(andreyev2009))] <- 1
andreyev2009$Rep[grepl("2", sampleNames(andreyev2009))] <- 2
andreyev2009$Rep[grepl("3", sampleNames(andreyev2009))] <- 3
featureNames(andreyev2009) <- fData(andreyev2009)[, 2]
andreyev2009 <- impute(andreyev2009, "zero")



mm <- read.csv("../../inst/extdata/Andreyev2009_Supplemental_Table_S2-L-C-EAD.csv.gz",
               skip = 2, stringsAsFactors = FALSE)
names(mm) <- c("index", "GeneSymbol", "Accession", "Coverage",
               "UniquePepNum", "TotalScore", "ProtName", "markers")
mrks <- mm$markers
names(mrks) <- make.unique(mm[, "Accession"])

andreyev2009 <- addMarkers(andreyev2009, mrks)


sel <- andreyev2009$condition == "Resting"
andreyev2009rest <- normalise(andreyev2009[, sel], "sum")
andreyev2009activ <- normalise(andreyev2009[, !sel], "sum")

save(andreyev2009, file = "../../data/andrevey2009.rda",
     compression_level = 9, compress = "xz")
save(andreyev2009rest, file = "../../data/andrevey2009rest.rda",
     compression_level = 9, compress = "xz")
save(andreyev2009activ, file = "../../data/andrevey2009activ.rda",
     compression_level = 9, compress = "xz")
