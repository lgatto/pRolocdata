library("MSnbase")
library("pRoloc")

f <- "../extdata/journal.pbio.2004411.s003.tab3.csv"
i <- grepEcols(f, split = ",", pattern = "Ratio")

## Prepare and annotate data

hirst2018 <- readMSnSet2(f, i)
fvarLabels(hirst2018)[4] <- "markers"
hirst2018 <- fDataToUnknown(hirst2018)
hirst2018 <- normalise(hirst2018, method = "sum")

featureNames(hirst2018) <- fData(hirst2018)[, 1]

sampleNames(hirst2018) <- sub("Ratio\\.H\\.L\\.", "",
                              sampleNames(hirst2018))

hirst2018$rep <- as.numeric(sub("_.+$", "",
                                sub("R", "", sampleNames(hirst2018))))
hirst2018$centrifugation <- factor(paste0(sub("^.+_", "",
                                              sampleNames(hirst2018)), "K"))
sample <- rep("CTRL", ncol(hirst2018))
sample[grep("C2", sampleNames(hirst2018))] <- "C2"
sample[grep("C6", sampleNames(hirst2018))] <- "C6"
hirst2018$sample <- factor(sample)

## Add publication information
miape <- new("MIAPE",
             abstract = "The AP-5 adaptor protein complex is presumed to function in membrane traffic, but so far nothing is known about its pathway or its cargo. We have used CRISPR-Cas9 to knock out the AP-5 zeta subunit gene, AP5Z1, in HeLa cells, and then analysed the phenotype by subcellular fractionation profiling and quantitative mass spectrometry. The retromer complex had an altered steady-state distribution in the knockout cells, and several Golgi proteins, including GOLIM4 and GOLM1, were depleted from vesicle-enriched fractions. Immunolocalisation showed that loss of AP-5 led to impaired retrieval of the cation-independent mannose 6-phosphate receptor (CIMPR), GOLIM4, and GOLM1 from endosomes back to the Golgi region. Knocking down the retromer complex exacerbated this phenotype. Both the CIMPR and sortilin interacted with the AP-5-associated protein SPG15 in pull-down assays, and we propose that sortilin may act as a link between Golgi proteins and the AP-5/SPG11/SPG15 complex. Together, our findings suggest that AP-5 functions in a novel sorting step out of late endosomes, acting as a backup pathway for retromer. This provides a mechanistic explanation for why mutations in AP-5/SPG11/SPG15 cause cells to accumulate aberrant endolysosomes, and highlights the role of endosome/lysosome dysfunction in the pathology of hereditary spastic paraplegia and other neurodegenerative disorders.",
             pubMedIds = "29381698",
             title = "Role of the AP-5 adaptor protein complex in late endosome-to-Golgi retrieval")
experimentData(hirst2018) <- miape


## Add feature metadata
fd <- read.csv("../extdata/journal.pbio.2004411.s003.tab4.csv",
               stringsAsFactors = FALSE)
suppressWarnings(fd2 <- dplyr::left_join(fData(hirst2018), fd))
rownames(fd2) <- fd2[, 1]
fData(hirst2018) <- fd2
fvarLabels(hirst2018)[8] <- "Hits"
fData(hirst2018)$Hits <- as.logical(fData(hirst2018)$Hits)

if (validObject(hirst2018))
    save(hirst2018, file = "../../data/hirst2018.rda",
         compress = "xz", compression_level = 9)
