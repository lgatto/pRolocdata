library("MSnbase")

makeSet <- function(filename) {
  x1 <- read.csv(filename, stringsAsFactors = FALSE)
  es <- as.matrix(x1[, grep("_11", names(x1))])
  fd <- new("AnnotatedDataFrame",
            data = x1[, -grep("_11", names(x1))])
  obj <- new("MSnSet",
             exprs = es,
             featureData = fd)  
  pd <- data.frame(experiment = paste0("LOPIT", rep(1:4, each = 4)),
                   replicate = rep(rep(c("A", "B"), each = 4), 4),
                   tags = rep("iTRAQ4", 32),
                   reporter = rep(114:117, 8),
                   row.names = sampleNames(obj))
  phenoData(obj) <- new("AnnotatedDataFrame", data = pd)
  .experiment <- new("MIAPE",
                     lab = "Cambridge Centre for Proteomics (CCP)",
                     name = "Kathryn S. Lilley",
                     email = "k.s.lilley@bioc.cam.ac.uk",
                     samples = list(
                       species = "Arabidopsis thaliana",
                       tissue = "Callus"),
                     title = "Putative glycosyltransferases and other plant golgi apparatus proteins are revealed by {LOPIT} proteomics.",
                     abstract = "The Golgi apparatus is the central organelle in the secretory pathway and plays key roles in glycosylation, protein sorting, and secretion in plants. Enzymes involved in the biosynthesis of complex polysaccharides, glycoproteins, and glycolipids are located in this organelle, but the majority of them remain uncharacterized. Here, we studied the Arabidopsis (Arabidopsis thaliana) membrane proteome with a focus on the Golgi apparatus using localization of organelle proteins by isotope tagging. By applying multivariate data analysis to a combined data set of two new and two previously published localization of organelle proteins by isotope tagging experiments, we identified the subcellular localization of 1,110 proteins with high confidence. These include 197 Golgi apparatus proteins, 79 of which have not been localized previously by a high-confidence method, as well as the localization of 304 endoplasmic reticulum and 208 plasma membrane proteins. Comparison of the hydrophobic domains of the localized proteins showed that the single-span transmembrane domains have unique properties in each organelle. Many of the novel Golgi-localized proteins belong to uncharacterized protein families. Structure-based homology analysis identified 12 putative Golgi glycosyltransferase (GT) families that have no functionally characterized members and, therefore, are not yet assigned to a Carbohydrate-Active Enzymes database GT family. The substantial numbers of these putative GTs lead us to estimate that the true number of plant Golgi GTs might be one-third above those currently annotated. Other newly identified proteins are likely to be involved in the transport and interconversion of nucleotide sugar substrates as well as polysaccharide and protein modification.",
                     pubMedIds = "22923678")
  experimentData(obj) <- .experiment
  obj@processingData@files <- basename(filename)
  if (validObject(obj))
    return(obj)
}

f1 <- "../extdata/Nikolovski2012_SupplTable1.csv"
nikolovski2012 <- makeSet(f1)
featureNames(nikolovski2012) <- fData(nikolovski2012)$Protein.ID


f2 <- "../extdata/Nikolovski2012_SupplTable2.csv"
nikolovski2012imp <- makeSet(f2)
featureNames(nikolovski2012imp) <- fData(nikolovski2012imp)$Protein.ID
sel <- fData(nikolovski2012imp)$markers == ""
fData(nikolovski2012imp)$markers[sel] <- "unknown"

idx <- match(featureNames(nikolovski2012),
             featureNames(nikolovski2012imp))
fData(nikolovski2012)$markers <- fData(nikolovski2012imp)$markers[idx]

stopifnot(all.equal(sort(featureNames(nikolovski2012)),
                    sort(featureNames(nikolovski2012imp))))

stopifnot(validObject(nikolovski2012))
stopifnot(validObject(nikolovski2012imp))

save(nikolovski2012,
     file="../../data/nikolovski2012.RData",
     compress = "xz",
     compression_level = 9)
save(nikolovski2012imp,
     file="../../data/nikolovski2012imp.RData",
     compress = "xz",
     compression_level = 9)
