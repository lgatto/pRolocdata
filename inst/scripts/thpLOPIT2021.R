library("MSnbase")
library("pRoloc")

## =================== function to add experimental data =================== 
addExperimentInfo <- function(date = "Summer 2016/2017",
                              instrument = "Orbitrap Fusion Lumos Tribrid") {
  experiment <- new("MIAPE",
                    lab = "Cambridge Centre for Proteomics (CCP)",
                    name = "Claire M. Mulvey",
                    contact = "Kathryn S. Lilley",
                    email = "k.s.lilley@bioc.cam.ac.uk",
                    samples = list(
                      species = "thp1 cells",
                      operator = "Claire M. Mulvey"
                    ),
                    title = "Spatiotemporal proteomic profiling of the dynamic pro-inflammatory response to lipopolysaccharide in the THP-1 human leukaemia cell line",
                    abstract = "Protein localisation and translocation between intracellular compartments underlie key biological mechanisms. Recent advances in the spatial proteomics field allow comprehensive insights into dynamic (patho)physiological processes at the subcellular level. The hyperLOPIT proteomics platform combines mass-spectrometry with state-of-the-art machine learning for simultaneous “mapping” of the steady-state subcellular location of thousands of proteins. Here, we use a synergistic approach and combine global proteome analysis with hyperLOPIT in a fully Bayesian framework to elucidate the spatio-temporal changes that occur during the pro-inflammatory response to lipopolysaccharide in the human monocytic leukaemia cell-line THP-1.  We report cell-wide protein relocalisations upon LPS stimulation, including many proteins that have known roles in cell migration and the endo-lysosomal autophagy system. By harnessing and quantifying proteome-wide uncertainty through Bayesian modelling we are able to distinguish vital translocation events, revealing a necessary role for protein relocalisation through the most extensive insight into the LPS-driven innate immune response, to date.",
                    pubMedIds = "",
                    url = "",
                    instrumentModel = instrument,
                    instrumentManufacturer = "ThermoScientific",
                    ionSource = "",
                    analyser = "Orbitrap",
                    detectorType = "Orbitrap",
                    softwareName = "Mascot Search Engine",
                    collisionEnergy = "",
                    dateStamp = date
  )
}

## ======================= function to add Experiment info =======================
# makepdata <- function(data, ) {
#   sn <- sampleNames(data)
#   r1 <- grep("rep1", sn)
#   r2 <- grep("rep2", sn)
#   r3 <- grep("rep3", sn)
#   Replicate = rep(reps, each = 10)
#   
#   .pData <- data.frame(Replicate = rep(reps, each = 10),
#                        TMT.Reagent = sampleNames(msnset),
#                        Acquisiton.Method = method,
#                        row.names= sampleNames(msnset))
#   return(.pData)
# }

## populate the msnset slots
makeTHP <- function(filename_data,
                    filename_pdata,
                    ecols = c(2:60),
                    date = "Summer 2016/2017",
                    instrument = "Orbitrap Fusion Lumos Tribrid") {

  csv <- read.csv(filename_data)
  data <- readMSnSet2(file = filename_data, ecol = ecols, skip = 0, fnames = 1)
  
  ## remove redundant columns in fdata
  if (any(fvarLabels(data) == "markers.old")) {
    ind <- which(fvarLabels(data) == "markers.old")
    fData(data) <- fData(data)[, -ind]
  }
  if (any(fvarLabels(data) == "X")) {
    ind <- which(fvarLabels(data) == "X")
    fData(data) <- fData(data)[, -ind]
  }
  
  if (!missing(filename_pdata)) {
    pdat <- read.csv(filename_pdata)
    if(all(sampleNames(data) != pdat$X)) stop("data does not match")
    rownames(pdat) <- pdat$X
    pdat <- pdat[, -1]
    pData(data) <- pdat
  } else {
    data <- makepdata(data)
  }
  
  .experiment <- addExperimentInfo(date = date, instrument = instrument)
  experimentData(data) <- .experiment
  .process <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".",sep=""),
                  paste("Normalised to sum of intensities.")),
                normalised=TRUE)
  data@processingData <- .process
  # obj <- new("MSnSet",
  #            exprs = .exprs,
  #            phenoData = .pData,
  #            experimentData = .experiment,
  #            featureData = .fData)

  if (validObject(data))
  return (data)
}
makepdata <- function(msnset) { 
  nb <- strsplit(sampleNames(msnset), "_")
  pData(msnset)$Tag <- sapply(nb, function(z) z[4])
  pData(msnset)$Treatment <- sapply(nb, function(z) z[1])
  pData(msnset)$Replicate <- sapply(nb, function(z) z[2])
  pData(msnset)$Set <- sapply(nb, function(z) z[3])
  pData(msnset)$Fraction <- sapply(nb, function(z) z[5])
  return(msnset)
}


## ================================================================
## ======================MAKE DATA=================================
## ================================================================

## --------------full datasets------------------
## 3882 for the unstimulated concatenated
## 4067 for the THP stimulated set
file_full_unst <- "../../inst/extdata/thp1/thp1_LOPIT_unst_protein_full.csv.gz"
pdat_full_unst <- "../../inst/extdata/thp1/thp1_pdata_unst_full.csv.gz"

thpLOPIT_unstimulated_mulvey2021 <- makeTHP(filename_data = file_full_unst, 
                                        filename_pdata = pdat_full_unst, 
                                        date = "August 2016 - August 2017",
                                        instrument = "Orbitrap Fusion Lumos Tribrid")
stopifnot(validObject(thpLOPIT_unstimulated_mulvey2021))
save(thpLOPIT_unstimulated_mulvey2021, file = "../../data/thpLOPIT_unstimulated_mulvey2021.RData", 
     compress = "xz", compression_level = 9)

file_full_lps <- "../../inst/extdata/thp1/thp1_LOPIT_lps_protein_full.csv.gz"
pdat_full_lps <- "../../inst/extdata/thp1/thp1_pdata_lps_full.csv.gz"

thpLOPIT_lps_mulvey2021 <- makeTHP(filename_data = file_full_lps, 
                                        filename_pdata = pdat_full_lps, 
                                        date = "December 2016 - August 2017",
                                        instrument = "Orbitrap Fusion Lumos Tribrid")
stopifnot(validObject(thpLOPIT_lps_mulvey2021))
save(thpLOPIT_lps_mulvey2021, file = "../../data/thpLOPIT_lps_mulvey2021.RData", 
     compress = "xz", compression_level = 9)


## ============================================================
## ============================================================
## =========== PROTEIN LEVEL DATA - 3 REPLICATES===============
## ============================================================
## ============================================================

## 5107 proteins unst r1
fn1 <- "../../inst/extdata/thp1/thp1_LOPIT_unst_protein_r1.csv.gz"
thpLOPIT_unstimulated_rep1_mulvey2021 <- makeTHP(filename_data = fn1, 
                                             date = "August 2016",
                                             ecols = 2:21,
                                             instrument = "Orbitrap Fusion Lumos Tribrid")
stopifnot(validObject(thpLOPIT_unstimulated_rep1_mulvey2021))
save(thpLOPIT_unstimulated_rep1_mulvey2021, file = "../../data/thpLOPIT_unstimulated_rep1_mulvey2021.RData", 
     compress = "xz", compression_level = 9)

## 4838 proteins in unstim r2
fn2 <- "../../inst/extdata/thp1/thp1_LOPIT_unst_protein_r2.csv.gz"
thpLOPIT_unstimulated_rep2_mulvey2021 <- makeTHP(filename_data = fn2, 
                                             date = "October 2016",
                                             ecols = 2:20,
                                             instrument = "Orbitrap Fusion Lumos Tribrid")
stopifnot(validObject(thpLOPIT_unstimulated_rep2_mulvey2021))
save(thpLOPIT_unstimulated_rep2_mulvey2021, file = "../../data/thpLOPIT_unstimulated_rep2_mulvey2021.RData", 
     compress = "xz", compression_level = 9)

## 5733 proteins in unstim r3
fn3 <- "../../inst/extdata/thp1/thp1_LOPIT_unst_protein_r3.csv.gz"
thpLOPIT_unstimulated_rep3_mulvey2021 <- makeTHP(filename_data = fn3, 
                                             date = "August 2017",
                                             ecols = 2:21,
                                             instrument = "Orbitrap Fusion Lumos Tribrid")
stopifnot(validObject(thpLOPIT_unstimulated_rep3_mulvey2021))
save(thpLOPIT_unstimulated_rep3_mulvey2021, file = "../../data/thpLOPIT_unstimulated_rep3_mulvey2021.RData", 
     compress = "xz", compression_level = 9)

## 4879 proteins LPS r1
fn1 <- "../../inst/extdata/thp1/thp1_LOPIT_lps_protein_r1.csv.gz"
thpLOPIT_lps_rep1_mulvey2021 <- makeTHP(filename_data = fn1, 
                                             date = "December 2016",
                                             ecols = 2:21,
                                             instrument = "Orbitrap Fusion Lumos Tribrid")
stopifnot(validObject(thpLOPIT_lps_rep1_mulvey2021))
save(thpLOPIT_lps_rep1_mulvey2021, file = "../../data/thpLOPIT_lps_rep1_mulvey2021.RData", 
     compress = "xz", compression_level = 9)

## 4866 proteins in lpsim r2
fn2 <- "../../inst/extdata/thp1/thp1_LOPIT_lps_protein_r2.csv.gz"
thpLOPIT_lps_rep2_mulvey2021 <- makeTHP(filename_data = fn2, 
                                             date = "April 2017",
                                             ecols = 2:21,
                                             instrument = "Orbitrap Fusion Lumos Tribrid")
stopifnot(validObject(thpLOPIT_lps_rep2_mulvey2021))
save(thpLOPIT_lps_rep2_mulvey2021, file = "../../data/thpLOPIT_lps_rep2_mulvey2021.RData", 
     compress = "xz", compression_level = 9)

## 5848 proteins in lpsim r3
fn3 <- "../../inst/extdata/thp1/thp1_LOPIT_lps_protein_r3.csv.gz"
thpLOPIT_lps_rep3_mulvey2021 <- makeTHP(filename_data = fn3, 
                                             date = "August 2017",
                                             ecols = 2:21,
                                             instrument = "Orbitrap Fusion Lumos Tribrid")
stopifnot(validObject(thpLOPIT_lps_rep3_mulvey2021))
save(thpLOPIT_lps_rep3_mulvey2021, file = "../../data/thpLOPIT_lps_rep3_mulvey2021.RData", 
     compress = "xz", compression_level = 9)


## ============================================================
## ====================PSM LEVEL DATA =========================
## ============ 2 CONDITIONS X 3 REPS X 2 TMT SETS=============
## ================  A TOTAL OF 12 DATASETS ===================
## ============================================================

f <- paste0("../../inst/extdata/thp1/", list.files("../../inst/extdata/thp1/"))
f <- f[grep("1218.csv.gz", f)]

## create msnsets
sel <- 41:50
psms <- sapply(1:length(f), function(z) readMSnSet2(f[z], 41:50, 
                                                    stringsAsFactors = FALSE))

## Note: Claire does not want to use channel X126 in unstimulated rep2 set2 as something went wrong with this
# sampleNames(psms[[10]])
psms[[10]] <- psms[[10]][, -1]
# sampleNames(psms[[10]])
fvarLabels(psms[[10]])
fData(psms[[10]])$Average.Reporter.S.N <- fData(psms[[10]])$new.ave.rep.SN
fData(psms[[10]])$new.ave.rep.SN <- NULL

## update sampleNames with fraction info as per previous datasets sent november 2018
sn <- read.table("../../inst/extdata/thp1/sampleNames_timecourse.csv.gz")
sn <- sn[,1]
sampleNames(psms[[1]]) <- sn[1:10]
sampleNames(psms[[2]]) <- sn[11:20]
sampleNames(psms[[3]]) <- sn[21:30]
sampleNames(psms[[4]]) <- sn[31:40]
sampleNames(psms[[5]]) <- sn[41:50]
sampleNames(psms[[6]]) <- sn[51:60]
sampleNames(psms[[7]]) <- sn[61:70]
sampleNames(psms[[8]]) <- sn[71:80]
sampleNames(psms[[9]]) <- sn[81:90]
sampleNames(psms[[10]]) <- sn[91:99]
sampleNames(psms[[11]]) <- sn[100:109]
sampleNames(psms[[12]]) <- sn[110:119]

## add pData
makeTC_pdata <- function(msnset, tr, repNo, setNo) {
  nb <- strsplit(sampleNames(msnset), "_")
  pData(msnset)$Tag <- sapply(nb, function(z) z[4])
  pData(msnset)$Treatment <- tr
  pData(msnset)$Replicate <- repNo
  pData(msnset)$Set <- setNo
  pData(msnset)$Fraction <- sapply(nb, function(z) z[5])
  return(msnset)
}

## Count number of PSMs with the same sequence
# getPSMsCount <- function(object) {
#   fData(object)$Annotated.Sequence <- toupper(fData(object)$Annotated.Sequence)
#   seqs <- fData(object)$Annotated.Sequence
#   counts <- sapply(seqs, function(z) length(which(z == seqs)))
#   fData(object)$PSM.count <- counts
#   return(object)
# } 

## record the number of NA's across the fractions per PSM
recordPSMsNA <- function(object) {
  r <- apply(exprs(object), 1, function(z) sum(is.na(z)))
  fData(object)$PSM.count.na <- r
  return(object)
}

## record the maximum number of NA's in a PSM for a given peptide
## e.g. we have 7 PSMs, 1 of which had 2 missing values, 3 of which 
## had 1 missing value, we record the maximum number of NA's here as 3
maxPSMsNA <- function(object) {
  fData(object)$Annotated.Sequence <- toupper(fData(object)$Annotated.Sequence)
  seqs <- fData(object)$Annotated.Sequence
  ind <- lapply(seqs, function(z) 
    which(z == fData(object)$Annotated.Sequence))
  mymax <- sapply(ind, function(z) max(fData(object)[z, "PSM.count.na"]))
  fData(object)$PSM.count.na.max <- mymax
  return(object)
} 

t <- rep(c("LPS", "Unstimulated"), each = 6)
r <- rep(1:3, 2, each = 2)
s <- rep(1:2, 6)
for (i in 1:length(t)) {
  psms[[i]] <- makeTC_pdata(psms[[i]], t[i], r[i], s[i])
}

## Claire has imputed 0.01 for every missing value. We need to return these to NA's.
lapply(psms, function(z) exprs(z)[exprs(z) == 0.01] <- NA)
for (i in 1:length(psms)) {
  exprs(psms[[i]])[exprs(psms[[i]]) == 0.01] <- NA
}
psms <- lapply(psms, function(z) recordPSMsNA(z))

## record max number of missing values i.e. if 3 PSMs and 2 have 1 
## missing values, and 1 has no missing values, we record 1 as the 
## max (so missing value information is not lost when combining to
## peptide)
# psms <- lapply(psms, function(z) maxPSMsNA(z))

psms_thpLOPIT_lps_rep1_set1 <- psms[[1]]
stopifnot(validObject(psms_thpLOPIT_lps_rep1_set1))
save(psms_thpLOPIT_lps_rep1_set1, file = "../../data/psms_thpLOPIT_lps_rep1_set1.RData", 
     compress = "xz", compression_level = 9)
psms_thpLOPIT_lps_rep1_set2 <- psms[[2]]
stopifnot(validObject(psms_thpLOPIT_lps_rep1_set2))
save(psms_thpLOPIT_lps_rep1_set2, file = "../../data/psms_thpLOPIT_lps_rep1_set2.RData", 
     compress = "xz", compression_level = 9)
psms_thpLOPIT_lps_rep2_set1 <- psms[[3]]
stopifnot(validObject(psms_thpLOPIT_lps_rep2_set1))
save(psms_thpLOPIT_lps_rep2_set1, file = "../../data/psms_thpLOPIT_lps_rep2_set1.RData", 
     compress = "xz", compression_level = 9)
psms_thpLOPIT_lps_rep2_set2 <- psms[[4]]
stopifnot(validObject(psms_thpLOPIT_lps_rep2_set2))
save(psms_thpLOPIT_lps_rep2_set2, file = "../../data/psms_thpLOPIT_lps_rep2_set2.RData", 
     compress = "xz", compression_level = 9)
psms_thpLOPIT_lps_rep3_set1 <- psms[[5]]
stopifnot(validObject(psms_thpLOPIT_lps_rep3_set1))
save(psms_thpLOPIT_lps_rep3_set1, file = "../../data/psms_thpLOPIT_lps_rep3_set1.RData", 
     compress = "xz", compression_level = 9)
psms_thpLOPIT_lps_rep3_set2 <- psms[[6]]
stopifnot(validObject(psms_thpLOPIT_lps_rep3_set2))
save(psms_thpLOPIT_lps_rep3_set2, file = "../../data/psms_thpLOPIT_lps_rep3_set2.RData", 
     compress = "xz", compression_level = 9)

psms_thpLOPIT_unstim_rep1_set1 <- psms[[7]]
stopifnot(validObject(psms_thpLOPIT_unstim_rep1_set1))
save(psms_thpLOPIT_unstim_rep1_set1, file = "../../data/psms_thpLOPIT_unstim_rep1_set1.RData", 
     compress = "xz", compression_level = 9)
psms_thpLOPIT_unstim_rep1_set2 <- psms[[8]]
stopifnot(validObject(psms_thpLOPIT_unstim_rep1_set2))
save(psms_thpLOPIT_unstim_rep1_set2, file = "../../data/psms_thpLOPIT_unstim_rep1_set2.RData", 
     compress = "xz", compression_level = 9)
psms_thpLOPIT_unstim_rep2_set1 <- psms[[9]]
stopifnot(validObject(psms_thpLOPIT_unstim_rep2_set1))
save(psms_thpLOPIT_unstim_rep2_set1, file = "../../data/psms_thpLOPIT_unstim_rep2_set1.RData", 
     compress = "xz", compression_level = 9)
psms_thpLOPIT_unstim_rep2_set2 <- psms[[10]]
stopifnot(validObject(psms_thpLOPIT_unstim_rep2_set2))
save(psms_thpLOPIT_unstim_rep2_set2, file = "../../data/psms_thpLOPIT_unstim_rep2_set2.RData", 
     compress = "xz", compression_level = 9)
psms_thpLOPIT_unstim_rep3_set1 <- psms[[11]]
stopifnot(validObject(psms_thpLOPIT_unstim_rep3_set1))
save(psms_thpLOPIT_unstim_rep3_set1, file = "../../data/psms_thpLOPIT_unstim_rep3_set1.RData", 
     compress = "xz", compression_level = 9)
psms_thpLOPIT_unstim_rep3_set2 <- psms[[12]]
stopifnot(validObject(psms_thpLOPIT_unstim_rep3_set2))
save(psms_thpLOPIT_unstim_rep3_set2, file = "../../data/psms_thpLOPIT_unstim_rep3_set2.RData", 
     compress = "xz", compression_level = 9)


## =====================Have a look at missing values=================
# xx <- lapply(psms, function(z) which(apply(exprs(z), 1, anyNA)))
# for (i in seq(xx)) {
#   message("Replicate ", i, " contains ", length(xx[[i]]), 
#           " psms with missing values")
#   naplot(psms[[i]][xx[[i]], ], col = "black", main = paste0("replicate ", i), las = 2)
# }
# 
# If we examine the protein group of each PSM with a missing value, we find that t
# here are other PSMs available that can be used for quantitation for most proteins. 
# There are still several hundred cases where we only have 1 PSM for a given protein 
# group, so we would lose several hundred proteins per replicate if we were to just 
# remove missing values.  
# 
# numpeps <- sapply(seq(psms), function(y)
# {sapply(fData(psms[[y]])$Master.Protein.Accessions[xx[[y]]],
#         function(z) length(which(fData(psms[[y]])$Master.Protein.Accessions 
#                                  == z)))})
# sapply(numpeps, function(z) length(which(z == 1)))
# 
# So from the above we see that we have 564 singeton PSMs with a missing value in 
# replicate 1 set 1 (LPS), 337 in replicate 1 set 2 (LPS), 438 in replicate 2 set 1 (LPS) etc. 
# Let's have a look at all missing values (not just the ones for the above singletons) 
# to see if there is a trend in where they are missing.
# 
# par(mar=c(10,4,2,2), mfrow = c(4,3))
# for (i in 1:length(psms)) 
#   barplot(apply(exprs(psms[[i]]), 2, function(z) length(which(is.na(z)))), las = 2, 
#           main = strsplit(sampleNames(psms[[i]])[1], split = "_1")[[1]][1])
# 
# The barplots suggest missing values mostly accumulate in the first few fractions - these are 
# not missing at random according to Claire and are reflected from the gradient distributions.
# 
# 
## 
# psms.imputed <- lapply(psms, function(z) impute(z, method = "MinDet"))
# 
# Now check do the imputed values look okay? Yes
# 
# 1. Normalise
# 2. Combine features to protein
# psms.imputed.norm <- lapply(psms.imputed, function(z) normalise(z,  "sum"))
# 
# prots <- lapply(psms.imputed.norm, function(z) 
#   combineFeatures(z, 
#                   fun = "median",
#                   groupBy = fData(z)[, "Master.Protein.Accessions"],
#                   verbose = FALSE))
# 
# 