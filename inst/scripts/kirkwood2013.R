library("MSnbase")
library("pRoloc")

i <- 2:35
kirkwood2013 <- readMSnSet2("../extdata/24043423.txt.gz",
                            sep = "\t", ecol = i, skip = 1)
sampleNames(kirkwood2013) <-
    sub("X", "Fraction", sampleNames(kirkwood2013))

fvarLabels(kirkwood2013)[1] <- "AC"
fvarLabels(kirkwood2013)[36:40] <- c("Peptide.Count",
                                     "Replicate.1.Peptide.Count",
                                     "Replicate.2.Peptide.Count",
                                     "Replicate.3.Peptide.Count",
                                     "PEP")
featureNames(kirkwood2013) <- fData(kirkwood2013)$AC

experimentData(kirkwood2013) <- 
    new("MIAPE",
        title = "Characterization of native protein complexes and protein isoform variation using size-fractionation-based quantitative proteomics.",
        name = "Kirkwood KJ1, Ahmad Y, Larance M and Lamond AI",
        abstract = "Proteins form a diverse array of complexes that mediate cellular function and regulation. A largely unexplored feature of such protein complexes is the selective participation of specific protein isoforms and/or post-translationally modified forms. In this study, we combined native size-exclusion chromatography (SEC) with high-throughput proteomic analysis to characterize soluble protein complexes isolated from human osteosarcoma (U2OS) cells. Using this approach, we have identified over 71,500 peptides and 1,600 phosphosites, corresponding to over 8,000 proteins, distributed across 40 SEC fractions. This represents >50% of the predicted U2OS cell proteome, identified with a mean peptide sequence coverage of 27% per protein. Three biological replicates were performed, allowing statistical evaluation of the data and demonstrating a high degree of reproducibility in the SEC fractionation procedure. Specific proteins were detected interacting with multiple independent complexes, as typified by the separation of distinct complexes for the MRFAP1-MORF4L1-MRGBP interaction network. The data also revealed protein isoforms and post-translational modifications that selectively associated with distinct subsets of protein complexes. Surprisingly, there was clear enrichment for specific Gene Ontology terms associated with differential size classes of protein complexes. This study demonstrates that combined SEC/MS analysis can be used for the system-wide annotation of protein complexes and to predict potential isoform-specific interactions. All of these SEC data on the native separation of protein complexes have been integrated within the Encyclopedia of Proteome Dynamics, an online, multidimensional data-sharing resource available to the community.",
        pubMedIds = "24043423",
        samples = list("U2OS cells (Human)"),
        other = list(
            Quantification = "Label-Free (Spectral Counting)",
            Fractionation = "Size Exclusion Chromatography",
            Note = "High redundancy of proteins due to the database search (several isoforms identified for a unique protein)"))


load("../extdata/protein-complex-markers.rda")

kirkwood2013 <- addMarkers(kirkwood2013, pcmrk)

stopifnot(validObject(kirkwood2013))

save(kirkwood2013, file = "../../data/kirkwood2013.rda")
