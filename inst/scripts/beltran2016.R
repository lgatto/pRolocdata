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


for (i in seq_along(csvfls)) {
    x <- readMSnSet2(csvfls[i], ecol = 2:7)
    featureNames(x) <- fData(x)[, 1]
    x <- addMarkers(x, beltran2016markers, verbose = FALSE)
    fData(x) <- fData(x)[, -1, drop = FALSE]
    assign(nms[i], x)
    save(list = nms[i], file = rdafls[i],
         compress = "xz", compression_level = 9)
}
