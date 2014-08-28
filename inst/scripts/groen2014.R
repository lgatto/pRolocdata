library("MSnbase")
library("pRoloc")

makeGroen <- function(csvfile, fractions, rep) {
  xx <- read.csv(csvfile, row.names = 1, stringsAsFactors = FALSE)
  eset <- as.matrix(xx[, c(8:13)])
  sn <- c("114/115", "114/116",  "114/117", "115/116", 
          "115/117", "116/117")
  sn <- paste(sn, rep, sep=".")
  colnames(eset) <- sn
  fd <- xx[, c(1:4)]
  names(fd)[1] <- "Protein.Accession"
  names(fd)[2] <- "Protein.Description"
  names(fd)[3] <- "Q.Value"
  names(fd)[4] <- "No.Unique.Peptides"
  pd <- data.frame(Fractions = fractions,
                   row.names = colnames(eset))
  colnames(pd) <- paste("Fractions", rep, sep=".")
  exp <- new("MIAPE",
             lab = "Cambridge Centre for Proteomics (CCP)",
             name = "Arnoud J Groen",
             contact = "Kathryn S. Lilley",
             email = "k.s.lilley@bioc.cam.ac.uk",
             samples = list(
               species = "Arabidopsis thaliana",
               tissue = "callus"
               ),
             title = "Identification of Trans Golgi Network proteins in Arabidopsis thaliana root tissue",
             abstract = "Knowledge of protein subcellular localisation assists in the elucidation of protein function and understanding of different biological mechanisms which occur at discrete subcellular niches.  Organelle-centric proteomics, enables localisation of thousands of proteins simultaneously. Although such techniques have enabled successful organelle protein catalogues, they rely on purification or significant enrichment of the organelle of interest, which is not achievable for many organelles. Incomplete separation of organelles leads to false discoveries, with erroneous assignments. Proteomics methods that measure the distribution patterns of specific organelle markers along density gradients are able to assign proteins of unknown localisation based on co-migration, when coupled with sophisticated computational tools to determine organelle specific co-clusters, without organelle purification. 
Here we apply multiple approaches to establish a high confidence dataset of Arabidopsis trans-Golgi network (TGN) proteins. We employ traditional immuno-isolations of the TGN and a probability based organelle proteomics approach LOPIT (Localisation of Organelle Protein by Isotope Tagging) using density centrifugation and semi-supervised machine learning methods. This combined approach facilitates a reduction in false assignments to the TGN. It also enables proteins present in more than one organelle and cargo proteins en route to other cellular destinations to be distinguished from protein whose steady state location favours the TGN.",               
             url = "http://proteomics.bio.cam.ac.uk",
             instrumentModel = "LTQ Orbitrap Velos",
             instrumentManufacturer = "ThermoScientific",
             analyser = "Orbitrap",
             detectorType = "Orbitrap",
             other = list(
               datasetInformation = "Expt 1 from Groen et al"))
  fd.Ann <- new("AnnotatedDataFrame", data = fd)
  fd.Ann@varMetadata[,1] <- c("TAIR Arabidopsis Thaliana (AT) number",
                            "Protein description",
                            "Protein Q value (unique peptides only)",
                            "Number of unique peptides")
  ans <- new("MSnSet",
             exprs = eset,
             experimentData = exp,
             featureData = fd.Ann,
             phenoData = new("AnnotatedDataFrame", data = pd))
  featureNames(ans) <- toupper(featureNames(ans))
  if (validObject(ans))
    return(ans)
}
## Read marker file
y <- read.csv("../extdata/Groen_MarkerList.csv", row.names=1)
mrk <- as.character(y[,1])
names(mrk) <- rownames(y)
## Expt 1
f <- c("2/4", "2/6", "2/8-10", "4/6", "4/8-10", "6/8-10")
groen2014r1 <- makeGroen("../extdata/Groen_Supplemental_Table3a.csv", f, "r1")
groen2014r1 <- addMarkers(groen2014r1, mrk, verbose=FALSE)
## Expt 2A
f <- c("2/4", "2/6", "2/8-11", "4/6", "4/8-11", "6/8-11")
groen2014r2 <- makeGroen("../extdata/Groen_Supplemental_Table3b.csv", f, "r2")
groen2014r2 <- addMarkers(groen2014r2, mrk, verbose=FALSE)
## Expt 2B
f <- c("3/4", "3/5", "3/7", "4/5", "4/7", "5/7")
groen2014r3 <- makeGroen("../extdata/Groen_Supplemental_Table3c.csv", f, "r3")
groen2014r3 <- addMarkers(groen2014r3, mrk, verbose=FALSE)
## Combined
cmb <- combine(groen2014r1, updateFvarLabels(groen2014r2))
cmb <- filterNA(cmb)
cmb <- combine(cmb, updateFvarLabels(groen2014r3))
cmb <- filterNA(cmb)
fData(cmb) <- fData(cmb)[, c(1,2,5)]
groen2014cmb <- cmb

save(groen2014r1, file="../../data/groen2014r1.RData",
     compress = "xz", compression_level = 9)
save(groen2014r2, file="../../data/groen2014r2.RData",
     compress = "xz", compression_level = 9)
save(groen2014r3, file="../../data/groen2014r3.RData",
     compress = "xz", compression_level = 9)
save(groen2014cmb, file="../../data/groen2014cmb.RData",
     compress = "xz", compression_level = 9)


