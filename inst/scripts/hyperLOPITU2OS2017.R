library("MSnbase")

x <- read.csv("..//extdata/hyperLOPIT_U2OS_201702.csv")
names(x) <- gsub("(\\.)+", "\\.", names(x))

hyperLOPITU2OS2017 <- readMSnSet2(x, 2:41, fnames = "Accession")

sn <- sampleNames(hyperLOPITU2OS2017)
pd <- pData(hyperLOPITU2OS2017)
pd$replicate <- as.integer(sub("\\.TMT.+$", "", sub("Replicate", "", sn)))
pd$set <- as.integer(sub("\\..+$", "", sub("Replicate.+Set", "", sn)))
pd$tmt <- sub("^.+Set[1-2]\\.", "", sn)
pData(hyperLOPITU2OS2017) <- pd

sn <- sub("\\.Set", "S", sub("TMT\\.", "", sub("Replicate", "R", sn)))
sampleNames(hyperLOPITU2OS2017) <- sn

stopifnot(validObject(hyperLOPITU2OS2017))
save(hyperLOPITU2OS2017, file = "../../data/hyperLOPITU2OS2017.rda",
     compress = "xz", compression_level = 9)


i <- order(colSums(exprs(hyperLOPITU2OS2017)))[1:3]
hyperLOPITU2OS2017b <- hyperLOPITU2OS2017[, -i]

stopifnot(validObject(hyperLOPITU2OS2017b))
save(hyperLOPITU2OS2017b, file = "../../data/hyperLOPITU2OS2017b.rda",
     compress = "xz", compression_level = 9)
