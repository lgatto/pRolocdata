library("MSnbase")
library("readxl")
library("pRoloc")

## Spatial map
f <- "../extdata/1-s2.0-S2211124717311889-mmc3.xlsx"

x <- read_excel(f, sheet = 2)
i <- grep("LFQ intensity", names(x))
x <- x[, 1:max(i)]
itzhak2017 <- readMSnSet2(x, i)
sampleNames(itzhak2017) <- sub("LFQ intensity MC_", "", sampleNames(itzhak2017))

itzhak2017$map <- as.numeric(sub("_.+K", "",
                                 sub("Map", "", sampleNames(itzhak2017))))
itzhak2017$centrifugation <- sub("^.+_", "", sampleNames(itzhak2017))

miape <- new("MIAPE", pubMedIds = "28903049",
             title = "A Mass Spectrometry-Based Approach for Mapping Protein Subcellular Localization Reveals the Spatial Proteome of Mouse Primary Neurons",
             abstract = "We previously developed a mass spectrometry-based method, dynamic organellar maps, for the determination of protein subcellular localization and identification of translocation events in comparative experiments. The use of metabolic labeling for quantification (stable isotope labeling by amino acids in cell culture [SILAC]) renders the method best suited to cells grown in culture. Here, we have adapted the workflow to both label-free quantification (LFQ) and chemical labeling/multiplexing strategies (tandem mass tagging [TMT]). Both methods are highly effective for the generation of organellar maps and capture of protein translocations. Furthermore, application of label-free organellar mapping to acutely isolated mouse primary neurons provided subcellular localization and copy-number information for over 8,000 proteins, allowing a detailed analysis of organellar organization. Our study extends the scope of dynamic organellar maps to any cell type or tissue and also to high-throughput screening.")

itzhak2017@experimentData <- miape
itzhak2017 <- normalise(itzhak2017, method = "sum")
featureNames(itzhak2017) <- fData(itzhak2017)[, 2]

## Markers
x <- read_excel(f, sheet = 3)
i <- grep("nArea", names(x))
itzhak2017markers <- readMSnSet2(x, i)
fvarLabels(itzhak2017markers)[6] <- "markers"
featureNames(itzhak2017markers) <- fData(itzhak2017markers)[, 2]
sampleNames(itzhak2017markers) <- sub("nArea MC_", "",
                                      sampleNames(itzhak2017markers))
pData(itzhak2017markers) <- pData(itzhak2017)

mrk <- getMarkers(itzhak2017markers, verbose = FALSE)
itzhak2017 <- addMarkers(itzhak2017, mrk, verbose = FALSE)


if (validObject(itzhak2017))
    save(itzhak2017, file = "../../data/itzhak2017.rda",
         compress = "xz", compression_level = 9)

if (validObject(itzhak2017markers))
    save(itzhak2017markers, file = "../../data/itzhak2017markers.rda",
         compress = "xz", compression_level = 9)
