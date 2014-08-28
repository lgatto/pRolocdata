library("pRolocdata")

f <- dir("../extdata", full.names = TRUE, pattern = "245589")

## Expression data
## getEcols(f[1], ",")
e <- 9:28
nikolovski2014 <- readMSnSet2(f[1], e-1, row.names = 1)
fData(nikolovski2014)$markers <- 
    as.character(fData(nikolovski2014)[, "Subcellular.localisation.marker..known.localisation."])
fData(nikolovski2014)$markers[fData(nikolovski2014)$markers == " "] <- "unknown"
fData(nikolovski2014)[, "Subcellular.localisation.marker..known.localisation."] <- NULL
pData(nikolovski2014) <- 
    data.frame(gradient = substr(sampleNames(nikolovski2014), 1, 1), 
               fraction = as.numeric(sub("[AB]_", "", sampleNames(nikolovski2014))),
               row.names = sampleNames(nikolovski2014))

experimentData(nikolovski2014) <-
    new("MIAPE",
        lab = "Cambridge Centre for Proteomics (CCP)",
        url = "http://proteomics.bio.cam.ac.uk/",
        name = "Kathryn S. Lilley",
        email = "k.s.lilley@bioc.cam.ac.uk",
        samples = list(
            species = "Arabidopsis thaliana",
            tissue = ""),
        title = "Label free protein quantification for plant Golgi protein localisation and abundance",
        abstract = "The proteomic composition of the Arabidopsis Golgi apparatus is currently reasonably well documented; however little is known about the relative abundances between different proteins within this compartment. Accurate quantitative information of Golgi resident proteins is of great importance: it facilitates a better understanding of the biochemical processes which take place within this organelle, especially those of different polysaccharide synthesis pathways. Golgi resident proteins are challenging to quantify since the abundance of this organelle is relatively low within the cell. In this study an organelle fractionation approach, targeting the Golgi apparatus, was combined with a label free quantitative mass spectrometry (MS), data-independent acquisition (DIA) method employing ion mobility separation known as LC-IMS-MSE (or HDMSE), to simultaneously localize proteins to the Golgi apparatus and assess their relative quantity. In total 102 Golgi localised proteins were quantified. These data provide new insight into Golgi apparatus organization and demonstrate that organelle fractionation in conjunction with label free quantitative MS is a powerful and relatively simple tool to access protein organelle localisation and their relative abundances. The findings presented open a unique view on the organization of the plant Golgi apparatus, leading towards novel hypotheses centered on the biochemical processes of this organelle. he proteomic composition of the Arabidopsis Golgi apparatus is currently reasonably well documented; however little is known about the relative abundances between different proteins within this compartment. Accurate quantitative information of Golgi resident proteins is of great importance: it facilitates a better understanding of the biochemical processes which take place within this organelle, especially those of different polysaccharide synthesis pathways. Golgi resident proteins are challenging to quantify since the abundance of this organelle is relatively low within the cell. In this study an organelle fractionation approach, targeting the Golgi apparatus, was combined with a label free quantitative mass spectrometry (MS), data-independent acquisition (DIA) method employing ion mobility separation known as LC-IMS-MSE (or HDMSE), to simultaneously localize proteins to the Golgi apparatus and assess their relative quantity. In total 102 Golgi localised proteins were quantified. These data provide new insight into Golgi apparatus organization and demonstrate that organelle fractionation in conjunction with label free quantitative MS is a powerful and relatively simple tool to access protein organelle localisation and their relative abundances. The findings presented open a unique view on the organization of the plant Golgi apparatus, leading towards novel hypotheses centered on the biochemical processes of this organelle.",
        pubMedIds = "25122472")



## Adding markers
markers <- read.csv(f[2], row.names = 1, stringsAsFactors = FALSE)
i <- match(rownames(markers), featureNames(nikolovski2014))
fData(nikolovski2014)$PMID <-
    fData(nikolovski2014)$Localisation.Method <- NA
fvarLabels(nikolovski2014)[8:9] <- 
    paste(fvarLabels(nikolovski2014)[8:9], "Markers", sep = ".")
fData(nikolovski2014)[i, c("Localisation.Method.Markers", "PMID.Markers")] <-
    markers[, c("Localisation.Method", "PMID")]

## Adding localisation results
res <- read.csv(f[3], row.names = 1)[, -(1:2)]
i <- match(rownames(res), featureNames(nikolovski2014))
fData(nikolovski2014)[, names(res)] <- NA
fData(nikolovski2014)[i, names(res)] <- res

stopifnot(validObject(nikolovski2014))

save(nikolovski2014, file = "../../data/nikolovski2014.rda",
     compress = "xz", compression_level = 9)
