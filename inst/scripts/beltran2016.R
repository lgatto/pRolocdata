library("MSnbase")
library("readxl")

fm <- "../extdata/1-s2.0-S2405471216302897-mmc4.xlsx"
m <- read_xlsx(fm, sheet = 1, skip = 1)
beltran2016markers <- m[[3]]
names(beltran2016markers) <- m[[1]]

## These data are corrupted: some 31 values are formatted as number
## and dates (see [1] for details). An email has been sent to the
## authors. In the meantime, these phony values are set to NA.
## [1] https://twitter.com/lgatt0/status/974756030793043969)

f1 <- "../extdata/1-s2.0-S2405471216302897-mmc3.xlsx"

tabs <- c(paste("HCMV", paste0(seq(24, 120, 24), "hpi")),
          paste("Mock", paste0(seq(24, 120, 24), "hpi")))

coltypes <- c("text", "text", "text", "numeric", "numeric", "numeric",
              "numeric", "numeric", "numeric")

mlst <- lapply(tabs, function(ti) {
    x <- read_xlsx(f1, sheet = ti, skip = 2, col_types = coltypes)
    ms <- readMSnSet2(x, ecol = 4:9)
    fvarLabels(ms) <- gsub(" ", "_", fvarLabels(ms))
    featureNames(ms) <- fData(ms)$Uniprot_Accession
    addMarkers(ms, beltran2016markers, verbose = FALSE)
})

nms <- paste0("beltran2016", sub("hpi", "", sub(" ", "", tabs)))
fls <- file.path("../../data", paste0(nms, ".rda"))

for (i in 1:10) {
    assign(nms[i], mlst[[i]])
    save(list = nms[i], file = fls[i],
         compress = "xz", compression_level = 9)
}
