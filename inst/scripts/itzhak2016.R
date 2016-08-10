## download it if not available in the current working directory
url <- "https://elife-publishing-cdn.s3.amazonaws.com/16950/elife-16950-supp9-v3-download.xlsx"
f <- basename(url)
f <- file.path("../extdata", f)
if (!file.exists(f))
    download.file(url, f)


library("readxl")
xx <- read_excel(f, sheet = 2)
suppressPackageStartupMessages(library("pRoloc"))
itzhak2016stcSILAC <- readMSnSet2(xx,
                                  ecol = grep("log", colnames(xx)),
                                  fnames = 2)
itzhak2016stcSILAC$rep <- as.numeric(sub("_.*$", "",
                                         sub("log.+MAP", "",
                                             sampleNames(itzhak2016stcSILAC))))
fData(itzhak2016stcSILAC)$markers <-
                            fData(itzhak2016stcSILAC)[, "Organellar markers"]
stopifnot(packageVersion("pRoloc") >= "1.13.8") ## support for 'from = NA'

itzhak2016stcSILAC <- fDataToUnknown(itzhak2016stcSILAC,
                                     from = NA, to = "unknown")

stopifnot(validObject(itzhak2016stcSILAC))
save(itzhak2016stcSILAC,
     file = "../extdata/itzhak2016stcSILAC.rda",
     compress = "xz", compression_level = 9)
