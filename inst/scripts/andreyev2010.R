library("MSnbase")
library("pRoloc")
library("readr")

xx <- read_csv("../../inst/extdata/Andreyev2009_Supplemental_Table_S1-L-C.csv.gz",
               skip = 1)
xx <- xx[, !is.na(names(xx))] ## get rid of empty columns
names(xx)[2] <- "Accession"
    
andreyev2010 <- readMSnSet2(as.data.frame(xx), ecol = 7:78)
andreyev2010$fractions <- sub("-.+$", "", sampleNames(andreyev2010))
andreyev2010$condition <- ifelse(grepl("-C", sampleNames(andreyev2010)),
                                 "Resting", "Endotoxin-activated")
andreyev2010$Var1 <- ifelse(grepl("a$", sampleNames(andreyev2010)), "A", "B")
andreyev2010$Rep <- 0
andreyev2010$Rep[grepl("1", sampleNames(andreyev2010))] <- 1
andreyev2010$Rep[grepl("2", sampleNames(andreyev2010))] <- 2
andreyev2010$Rep[grepl("3", sampleNames(andreyev2010))] <- 3
featureNames(andreyev2010) <- fData(andreyev2010)[, 2]
andreyev2010 <- impute(andreyev2010, "zero")



mm <- read.csv("../../inst/extdata/Andreyev2009_Supplemental_Table_S2-L-C-EAD.csv.gz",
               skip = 2, stringsAsFactors = FALSE)
names(mm) <- c("index", "GeneSymbol", "Accession", "Coverage",
               "UniquePepNum", "TotalScore", "ProtName", "markers")
mrks <- mm$markers
names(mrks) <- make.unique(mm[, "Accession"])

andreyev2010 <- addMarkers(andreyev2010, mrks)


sel <- andreyev2010$condition == "Resting"
andreyev2010rest <- normalise(andreyev2010[, sel], "sum")
andreyev2010activ <- normalise(andreyev2010[, !sel], "sum")

stopifnot(validObject(andreyev2010activ))
stopifnot(validObject(andreyev2010rest))
stopifnot(validObject(andreyev2010))

save(andreyev2010, file = "../../data/andreyev2010.rda",
     compression_level = 9, compress = "xz")
save(andreyev2010rest, file = "../../data/andreyev2010rest.rda",
     compression_level = 9, compress = "xz")
save(andreyev2010activ, file = "../../data/andreyev2010activ.rda",
     compression_level = 9, compress = "xz")
