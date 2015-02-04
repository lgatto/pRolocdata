library("MSnbase")
library("pRoloc")

## getEcols("../extdata/25561571_Rep1.txt", "\t")
fabre2015r1 <- readMSnSet2("../extdata/25561571_Rep1.txt.gz",
                           ecol = 5:23, sep = "\t", fnames = "AC")
## getEcols("../extdata/25561571_Rep2.txt", "\t")
fabre2015r2 <- readMSnSet2("../extdata/25561571_Rep2.txt.gz",
                           ecol = 6:24, sep = "\t", fnames = "AC")

experimentData(fabre2015r1) <-
    experimentData(fabre2015r2) <-
        new("MIAPE",
            title = "Deciphering preferential interactions within supramolecular protein complexes: the proteasome case.",
            abstract = "In eukaryotic cells, intracellular protein breakdown is mainly performed by the ubiquitin-proteasome system. Proteasomes are supramolecular protein complexes formed by the association of multiple sub-complexes and interacting proteins. Therefore, they exhibit a very high heterogeneity whose function is still not well understood. Here, using a newly developed method based on the combination of affinity purification and protein correlation profiling associated with high-resolution mass spectrometry, we comprehensively characterized proteasome heterogeneity and identified previously unknown preferential associations within proteasome sub-complexes. In particular, we showed for the first time that the two main proteasome subtypes, standard proteasome and immunoproteasome, interact with a different subset of important regulators. This trend was observed in very diverse human cell types and was confirmed by changing the relative proportions of both 20S proteasome forms using interferon-gamma. The new method developed here constitutes an innovative and powerful strategy that could be broadly applied for unraveling the dynamic and heterogeneous nature of other biologically relevant supramolecular protein complexes.",
            name = "Bertrand Fabre et al.",
            pubMedIds = "25561571",
            samples  = list('U937 cells (Human)'),
            other = list(Quantification = 'Label-Free (iBAQ)',
                Fractionation = 'Glycerol Density Gradient Centrifugation'))

load("../extdata/protein-complex-markers.rda")

fabre2015r1 <- addMarkers(fabre2015r1, pcmrk)
fabre2015r2 <- addMarkers(fabre2015r2, pcmrk)

stopifnot(validObject(fabre2015r1))
stopifnot(validObject(fabre2015r2))

save(fabre2015r1, file = "../../data/fabre2015r1.rda")
save(fabre2015r2, file = "../../data/fabre2015r2.rda")
