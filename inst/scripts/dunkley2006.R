library("MSnbase")

###########################################################
## From Dunkey et al. 2006 PMID:16618929
dunkley <- read.csv("../extdata/Dunkley2006.csv",row.names=1)
.exprs <- as.matrix(dunkley[,4:19])
.fData <- dunkley[,c(2:3)]
names(.fData) <-c("markers","assigned")
.fData$evidence <- rep("unknown",nrow(.fData))
.fData$evidence[grep("predicted",.fData$markers)] <- "predicted"
.fData$evidence[grep("known",.fData$markers)] <- "known"
.fData$evidence[grep("unknown",.fData$markers)] <- "unknown"
.fData$markers <- sub("known ","",.fData$markers)
.fData$markers <- as.factor(sub("predicted ","",.fData$markers))
.fData$method <- rep("PLSDA",nrow(.fData))
.fData$method[grep("\\*",.fData$assigned)] <- "manual"
.fData$assigned <- sub("\\*","",.fData$assigned)
tmp <- as.character(.fData$assigned)
tmp[tmp=="plasma membrane"] <- "PM"
tmp[tmp=="not classified"] <- "unknown"
.fData$assigned <- as.factor(tmp)
.fData$new <- rep("new",nrow(.fData))
.fData$new[.fData[,1]==.fData[,2]] <- "known"
.fData$new[.fData[,2]=="unknown"] <- "unknown"
.fData$new <- as.factor(.fData$new)
.fData$pd.2013 <- dunkley[,20]
.fData$pd.markers <- dunkley[,21]
.fData <- new("AnnotatedDataFrame",data=data.frame(as.matrix(.fData)))
.fData@varMetadata[,1] <- c("Original protein markers used in Dunkley et al 2006 PNAS paper",
                            "Protein localisation assigned by PLSDA and manually by Dunkley as described in Dunkley et al 2006 PNAS paper",
                            "Evidence for original proteins markers in column 'markers'",
                            "Method of assignment for each protein in column 'assigned'",
                            "Protein status following new assignment",
		            "PhenoDisco output as described in Breckels et al (2013) Journal of Proteomics. Accepted February 2013",	                          
                            "Updated protein markers (original markers see column 'markers' plus new protein markers found by phenoDisco and verified by Uniprot/literature as described in Breckels et al)")
.pData <- new("AnnotatedDataFrame",
              data=data.frame(membrane.prep=rep(1:2,each=8),
                fraction=rep(c(1,4,7,11,2,5,8,11),2),
                replicate=rep(rep(c("A","B"),each=4),2),
                row.names=names(dunkley)[4:19]))

.experiment <- new("MIAPE",
                   lab = "Cambridge Centre for Proteomics (CCP)",
                   name = "Kathryn S. Lilley",
                   email = "k.s.lilley@bioc.cam.ac.uk",
                   samples = list(
                     species = "Arabidopsis thaliana",
                     tissue = "Callus"),
                   title="Mapping the Arabidopsis organelle proteome.",
                   abstract = "A challenging task in the study of the secretory pathway is the identification and localization of new proteins to increase our understanding of the functions of different organelles. Previous proteomic studies of the endomembrane system have been hindered by contaminating proteins, making it impossible to assign proteins to organelles. Here we have used the localization of organelle proteins by the isotope tagging technique in conjunction with isotope tags for relative and absolute quantitation and 2D liquid chromatography for the simultaneous assignment of proteins to multiple subcellular compartments. With this approach, the density gradient distributions of 689 proteins from Arabidopsis thaliana were determined, enabling confident and simultaneous localization of 527 proteins to the endoplasmic reticulum, Golgi apparatus, vacuolar membrane, plasma membrane, or mitochondria and plastids. This parallel analysis of endomembrane components has enabled protein steady-state distributions to be determined. Consequently, genuine organelle residents have been distinguished from contaminating proteins and proteins in transit through the secretory pathway.",
                   pubMedIds = "16618929",
                   url = "http://www.bio.cam.ac.uk/proteomics/",
                   instrumentModel = "QSTAR",
                   instrumentManufacturer = "Applied Biosystems",
                   ionSource = "ESI",
                   analyser = "TOF",
                   detectorType = "PMT")               

.process <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".",sep=""),
                  paste("Normalised to sum of intensities.")),
                 normalised=TRUE,
                files="Dunkley2006.csv")

dunkley2006 <- new("MSnSet",
                   exprs = .exprs,
                   phenoData = .pData,
                   experimentData = .experiment,
                   featureData = .fData)
dunkley2006@processingData <- .process

if (validObject(dunkley2006))
  save(dunkley2006,file="../../data/dunkley2006.RData",
       compress = "xz", compression_level = 9)

