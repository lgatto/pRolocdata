library("MSnbase")
library("pRoloc")

i <- grepEcols("../extdata/22939629.txt.gz",
               pattern = "Fraction", split = "\t")

havugimana2012 <-
    readMSnSet2("../extdata/22939629.txt.gz", sep = "\t",
                ecol = i, fnames = "AC")

fData(havugimana2012)$X <- NULL ## empty column in tsv file

experimentData(havugimana2012) <- 
    new("MIAPE",
        title = "A census of human soluble protein complexes",
        name = "Havugimana PC et al.",
        pubMedIds = "22939629",
        abstract = "Cellular processes often depend on stable physical associations between proteins. Despite recent progress, knowledge of the composition of human protein complexes remains limited. To close this gap, we applied an integrative global proteomic profiling approach, based on chromatographic separation of cultured human cell extracts into more than one thousand biochemical fractions that were subsequently analyzed by quantitative tandem mass spectrometry, to systematically identify a network of 13,993 high-confidence physical interactions among 3,006 stably associated soluble human proteins. Most of the 622 putative protein complexes we report are linked to core biological processes and encompass both candidate disease genes and unannotated proteins to inform on mechanism. Strikingly, whereas larger multiprotein assemblies tend to be more extensively annotated and evolutionarily conserved, human protein complexes with five or fewer subunits are far more likely to be functionally unannotated or restricted to vertebrates, suggesting more recent functional innovations.",
        samples = list("cytosolic and nuclear frations from HeLa and HEK cells (Human)"),
        other = list(
            Quantification = "Label-Free (Spectral Counting)",
            Fractionation = c(
                "HCW (Dual Heparin-IEX)",
                "SAF (Single phase-Heparin)", 
                "TCS (Triple phase-IEX)",
                "SGF (Sucrose gradient)",
                "IEF_pH 5 to 8 (Isoelectric focusing, pH 5 to 8)",
                "IEF_pH 3 to 10 (Isoelectric focusing, pH 3 to 10)",
                "WAX (Single phase-IEX)"),
            NOTE = c("extensive fractionation",
                "very high number of missing values")))


load("../extdata/protein-complex-markers.rda")

havugimana2012 <- addMarkers(havugimana2012, pcmrk)


## Nomenclature paper Havuginama et al:
## Example: LTQ_HeLaCE_WAX_Fraction01
## - LTQ or Orbitrap are for the mass-spectrometer type used to analyzed the sample
## - HeLaCE, HeLaNE and 293NE are the cell type from which the sample has
##   been produced and CE and NE stand for Cytosolic Extract and Nuclear
##   Extract respectively (a number can be added after CE or NE
##   corresponding to the number of biological replicate).
## - HCW (Dual Heparin-IEX), SAF (Single phase-Heparin), TCS (Triple
##   phase-IEX), SGF (Sucrose gradient), IEF_pH 5 to 8 (Isoelectric
##   focusing, pH 5 to 8), IEF_pH 3 to 10 (Isoelectric focusing, pH 3 to
##   10) and WAX (Single phase-IEX) are the separation methods used to
##   fractionate the samples.
## - Fraction correspond to the number of the fraction of a defined
##   separation method.
## - The order of the terms in the nomenclature is the following: MS
##   used_cell type and subcellular compartment_separation method
##   used_number of the fraction.
## 
## Example:
##   LTQ_HeLaCE12_TCS_Fraction238 correspond to Liner Trap Quadrupole_HeLa
##   Cytosolic Extract replicate 12_Triple phase-IEX_Fraction number 238.
##
## Summary table:
## |--------------------------------------+---------------------------------------+-----------------------+----------------+--------------------+---------------|
## | Cell and extracts                    | Fractionation approach                | # fractions collected | MS instrument  | # of MS replicates | Total MS runs |
## |--------------------------------------+---------------------------------------+-----------------------+----------------+--------------------+---------------|
## | HEK 293 nuclear extract (293NE)      | Dual Heparin-IEX (HCW)                |                   120 | LTQ instrument |                  2 |           240 |
## | HeLa S3  nuclear extract (HeLaNE)    | Single phase-Heparin (SAF)            |                    48 | LTQ instrument |                  1 |            48 |
## |                                      | Dual Heparin-IEX  (HCW)               |                   120 | LTQ instrument |                  2 |           240 |
## |                                      | Dual Heparin-IEX  (HCW)               |                   120 | LTQ-Orbitrap   |                  1 |           120 |
## |                                      | Triple phase-IEX (TCS)                |                   375 | LTQ instrument |                  2 |           750 |
## |                                      | Sucrose gradient (SGF)                |                    14 | LTQ-Orbitrap   |                  1 |            14 |
## |                                      | Isoelectric focusing (IEF_pH 5 to 8)  |                    10 | LTQ-Orbitrap   |                  1 |            10 |
## |                                      | Isoelectric focusing (IEF_pH 3 to 10) |                    10 | LTQ-Orbitrap   |                  1 |            10 |
## | La S3  cytoplasmic  extract (HeLaCE) | Isoelectric focusing (IEFpH 5 to 8)   |                    10 | LTQ-Orbitrap   |                  1 |            10 |
## |                                      | Isoelectric focusing (IEFpH 3 to 10)  |                    10 | LTQ-Orbitrap   |                  2 |            20 |
## |                                      | Sucrose gradient (SGF)                |                    14 | LTQ-Orbitrap   |                  1 |            14 |
## |                                      | Single phase-IEX (WAX)                |                    43 | LTQ instrument |                  1 |            43 |
## |                                      | Triple phase-IEX (TCS)                |                   269 | LTQ instrument |                  2 |           538 |
## |--------------------------------------+---------------------------------------+-----------------------+----------------+--------------------+---------------|

x <- sampleNames(havugimana2012)
x <- strsplit(x, "_Fraction")

pData(havugimana2012)$fraction <- as.integer(sapply(x, "[", 2))
x <- sapply(x, "[", 1)

pData(havugimana2012)$instrument <- sub("_.+$", "", x)
x <- sub("^[A-Za-z]+_", "", x)

cells <- sapply(strsplit(x, "_"), "[", 1)

celltype <- sapply(strsplit(cells, "E"), "[", 1)
celltype <- paste0(celltype, "E")

pData(havugimana2012)$celltype <- 
    ifelse(celltype == "293NE", "HEK293", "HeLa")

extract <- sapply(celltype,
                  function(x) {
                      x <- strsplit(x, "")[[1]]
                      lx <- length(x)
                      paste(x[(lx-1):lx], collapse = "")
                  })
pData(havugimana2012)$extract <-
    ifelse(extract == "NE", "Nuclear", "Cytosolic")

pData(havugimana2012)$rep <-
    as.numeric(sapply(strsplit(cells, "E"), "[", 2))

fractionation <- sub("^[A-Za-z0-9]+_", "", x)

fractionation[fractionation == "ph3_to_10"] <- "IEF_pH3_to_10"
fractionation[fractionation == "pH3_to10"] <- "IEF_pH3_to_10"

pData(havugimana2012)$fractionation <- fractionation

stopifnot(validObject(havugimana2012))

save(havugimana2012, file = "../../data/havugimana2012.rda")
