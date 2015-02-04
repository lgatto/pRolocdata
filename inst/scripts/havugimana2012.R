library("MSnbase")
library("pRoloc")

i <- grepEcols("../extdata/22939629.txt.gz", pattern = "Fraction", split = "\t")

havugimana2012 <- readMSnSet2("../extdata/22939629.txt.gz", sep = "\t",
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
            Fractionation = "Ion Exchange Chromatography, Isoelectric Focusing, Sucrose Density Gradient Centrifugation (extensive fractionation)",
            NOTE = "very high number of missing values"))


load("../extdata/protein-complex-markers.rda")

havugimana2012 <- addMarkers(havugimana2012, pcmrk)

stopifnot(validObject(havugimana2012))
save(havugimana2012, file = "../../data/havugimana2012.rda")
