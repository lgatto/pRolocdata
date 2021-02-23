library("MSnbase")
library("pRoloc")


## ============================================================
## ============================================================
## ==============LPS TIMECOURSE PSM LEVEL DATA=================
## ============================================================
## ============================================================

f <- "../../inst/extdata/thp1/LPS_timecourse_3reps_PSMs_301018.csv.gz"
allpsms <- readMSnSet2(f, 42:47, stringsAsFactors = FALSE)
# head(exprs(allpsms))
# dim(allpsms)
sf <- fData(allpsms)$Spectrum.File
psms1 <- allpsms[grep("THP1", sf),  ]
psms2 <- allpsms[grep("replicate1", sf),  ]
psms3 <- allpsms[grep("replicate2", sf),  ]
tags <- c("X126", "X127", "X128", "X129", "X130", "X131")
makepdata <- function(msnset, t, tr, repNo) { 
  pData(msnset)$Tag <- tags 
  pData(msnset)$Time <- t
  # pData(msnset)$Treatment <- tr
  pData(msnset)$Replicate <- rep(repNo, length(t))
  return(msnset)
}
psms1 <- makepdata(psms1, t = c(0, 2, 4, 6, 12, 24), 
                   r = 1)
psms2 <- makepdata(psms2, t = c(0, 2, 4, 6, 12, 24), 
                   r = 2)
psms3 <- makepdata(psms3, t = c(0, 2, 4, 6, 12, 24), 
                   r = 3)
## Make sampleNames treatment not tag
sn <- function(z) {
  sampleNames(z) <- paste(pData(z)[, 2], "hr", 
                          "rep", pData(z)[, 3])
  return(z)
}
psms1 <- sn(psms1)
psms2 <- sn(psms2)
psms3 <- sn(psms3)
getPSMs <- function(object) {
  fData(object)$Annotated.Sequence <- toupper(fData(object)$Annotated.Sequence)
  seqs <- fData(object)$Annotated.Sequence
  counts <- sapply(seqs, function(z) length(which(z == seqs)))
  fData(object)$PSM.count <- counts
  return(object)
} 
psms_lpsTimecourse_rep1_mulvey2021 <- getPSMs(psms1)
psms_lpsTimecourse_rep2_mulvey2021 <- getPSMs(psms2)
psms_lpsTimecourse_rep3_mulvey2021 <- getPSMs(psms3)

save(psms_lpsTimecourse_rep1_mulvey2021, file = "../../data/psms_lpsTimecourse_rep1_mulvey2021.RData", 
     compress = "xz", compression_level = 9)
stopifnot(validObject(psms_lpsTimecourse_rep1_mulvey2021))
save(psms_lpsTimecourse_rep2_mulvey2021, file = "../../data/psms_lpsTimecourse_rep2_mulvey2021.RData", 
     compress = "xz", compression_level = 9)
stopifnot(validObject(psms_lpsTimecourse_rep2_mulvey2021))
save(psms_lpsTimecourse_rep3_mulvey2021, file = "../../data/psms_lpsTimecourse_rep3_mulvey2021.RData", 
     compress = "xz", compression_level = 9)
stopifnot(validObject(psms_lpsTimecourse_rep3_mulvey2021))


## ============================================================
## ============================================================
## ============LPS TIMECOURSE PROTEIN LEVEL DATA===============
## ============================================================
## ============================================================

f <- "../../inst/extdata/thp1/LPS_timecourse_r1.csv.gz"
lpsTimecourse_rep1_mulvey2021 <- readMSnSet2(f, 2:7, stringsAsFactors = FALSE)
lpsTimecourse_rep1_mulvey2021 <- makepdata(lpsTimecourse_rep1_mulvey2021, 
                                           t = c(0, 2, 4, 6, 12, 24),
                                           r = 1)
featureNames(lpsTimecourse_rep1_mulvey2021) <- fData(lpsTimecourse_rep1_mulvey2021)$X
fData(lpsTimecourse_rep1_mulvey2021)$X <- NULL
save(lpsTimecourse_rep1_mulvey2021, file = "../../data/lpsTimecourse_rep1_mulvey2021.RData", 
     compress = "xz", compression_level = 9)
stopifnot(validObject(lpsTimecourse_rep1_mulvey2021))

f <- "../../inst/extdata/thp1/LPS_timecourse_r2.csv.gz"
lpsTimecourse_rep2_mulvey2021 <- readMSnSet2(f, 2:7, stringsAsFactors = FALSE)
lpsTimecourse_rep2_mulvey2021 <- makepdata(lpsTimecourse_rep2_mulvey2021, 
                                           t = c(0, 2, 4, 6, 12, 24),
                                           r = 2)
featureNames(lpsTimecourse_rep2_mulvey2021) <- fData(lpsTimecourse_rep2_mulvey2021)$X
fData(lpsTimecourse_rep2_mulvey2021)$X <- NULL
save(lpsTimecourse_rep2_mulvey2021, file = "../../data/lpsTimecourse_rep2_mulvey2021.RData", 
     compress = "xz", compression_level = 9)
stopifnot(validObject(lpsTimecourse_rep2_mulvey2021))

f <- "../../inst/extdata/thp1/LPS_timecourse_r3.csv.gz"
lpsTimecourse_rep3_mulvey2021 <- readMSnSet2(f, 2:7, stringsAsFactors = FALSE)
lpsTimecourse_rep3_mulvey2021 <- makepdata(lpsTimecourse_rep3_mulvey2021, 
                                           t = c(0, 2, 4, 6, 12, 24),
                                           r = 3)
featureNames(lpsTimecourse_rep3_mulvey2021) <- fData(lpsTimecourse_rep3_mulvey2021)$X
fData(lpsTimecourse_rep3_mulvey2021)$X <- NULL
save(lpsTimecourse_rep3_mulvey2021, file = "../../data/lpsTimecourse_rep3_mulvey2021.RData", 
     compress = "xz", compression_level = 9)
stopifnot(validObject(lpsTimecourse_rep3_mulvey2021))


## DATA from manuscript with limma results

f <- "../../inst/extdata/thp1/LPS_timecourse_allreps_with_stats.csv.gz"
lpsTimecourse_mulvey2021 <- readMSnSet2(f, 4:21, stringsAsFactors = FALSE)
featureNames(lpsTimecourse_mulvey2021) <- fData(lpsTimecourse_mulvey2021)$Accession
pData(lpsTimecourse_mulvey2021) <- rbind(pData(lpsTimecourse_rep1_mulvey2021), 
                                         pData(lpsTimecourse_rep2_mulvey2021), 
                                         pData(lpsTimecourse_rep3_mulvey2021))
save(lpsTimecourse_mulvey2021, file = "../../data/lpsTimecourse_mulvey2021.RData", 
     compress = "xz", compression_level = 9)
stopifnot(validObject(lpsTimecourse_mulvey2021))


