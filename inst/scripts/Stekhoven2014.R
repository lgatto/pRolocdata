library("MSnbase")
library("pRoloc")
library("magrittr")
stekhoven2014 <- readMSnSet2("../extdata/mmc3.csv", ecol = 16:37, sep = ",")
fData(stekhoven2014)$markers <-
                       fData(stekhoven2014)[, "selected.marker.protein"]

stekhoven2014 <- stekhoven2014 %>%
    normalise(method = "sum") %>%
    fDataToUnknown(from = "")

stopifnot(validObject(stekhoven2014))

save(stekhoven2014, file = "../../data/stekhoven2014.rda",
     compress = "xz", compression_level = 9)
