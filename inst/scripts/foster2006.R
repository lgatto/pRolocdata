library("MSnbase")
library("pRolocdata")
library("Vennerable")

###########################################################
## From Foster et al. 2006 PMID: 16615899

## To create peptide-level MSnSet, use 
## ../extdata/PIIS0092867406003692.mmc5-[high,low]Density.csv files
## and proceed as below

## Create protein-level MSnSet
pep2Prot <- read.csv("../extdata/PIIS0092867406003692.mmc2.csv.gz")
##pepInfo <- read.csv("../extdata/PIIS0092867406003692.mmc3.csv.gz")
protIntHD <- read.csv("../extdata/PIIS0092867406003692.mmc4-highDensity.csv.gz")
protIntLD <- read.csv("../extdata/PIIS0092867406003692.mmc4-lowDensity.csv.gz")
rownames(protIntHD) <- sub("\\.","",make.names(protIntHD[,"IPI"]))
rownames(protIntLD) <- sub("\\.","",make.names(protIntLD[,"IPI"]))

rnms <- unique(c(rownames(protIntHD),rownames(protIntLD)))
hfrac <- paste(grep("Fr",names(protIntHD),value=TRUE),"HD",sep="")
lfrac <- paste(grep("Fr",names(protIntLD),value=TRUE),"LD",sep="")
cnms <- c(hfrac,lfrac)

eset <- matrix(NA,length(rnms),length(cnms))
rownames(eset) <- rnms
colnames(eset) <- cnms

eset[rownames(protIntHD),hfrac] <- as.matrix(protIntHD[,grep("Fr",names(protIntHD))])
eset[rownames(protIntLD),lfrac] <- as.matrix(protIntLD[,grep("Fr",names(protIntLD))])

pdata <- data.frame(fractions=cnms,
                    density=c(rep("high",length(hfrac)),
                      rep("low",length(lfrac))),
                    num=sub("[H,L]D","",sub("Fr","",cnms)))
rownames(pdata) <- cnms
pdata <- new("AnnotatedDataFrame",data=pdata)

fdata <- matrix(NA,length(rnms),4)
rownames(fdata) <- rnms
colnames(fdata) <- c("UniProt","Name","NumPepHD","NumPepLD")
fdata[rownames(protIntHD),c("UniProt","Name","NumPepHD")] <- as.matrix(protIntHD[,c("UniProt","Name","NumPeptides")])
fdata[rownames(protIntLD),c("UniProt","Name","NumPepLD")] <- as.matrix(protIntLD[,c("UniProt","Name","NumPeptides")])
fdata <- data.frame(UniProt=fdata[,"UniProt"],
                    Name=fdata[,"Name"],
                    numPepHD=as.numeric(fdata[,"NumPepHD"]),
                    numPepLD=as.numeric(fdata[,"NumPepLD"]),
                    stringsAsFactors=FALSE)
## Training data (extracted manually)
fdata$markers <- "unknown"
fdata["IPI00269029","markers"] <- "Golgi"
fdata["IPI00118022","markers"] <- "PM"
fdata["IPI00453776","markers"] <- "EE"
fdata["IPI00223651","markers"] <- "TGN"
fdata["IPI00119618","markers"] <- "ER"
fdata["IPI00128071","markers"] <- "ERGDV"
fdata["IPI00113801","markers"] <- "mit"
fdata["IPI00113845","markers"] <- "PS" 
fdata["IPI00131845","markers"] <- "PS" 
fdata["IPI00131406","markers"] <- "PS" 

## Chi2 assignments
chi2vals <- read.csv("../extdata/PIIS0092867406003692.mmc6-localizations.csv.gz")
chi2loc <- read.csv("../extdata/PIIS0092867406003692.mmc6-refinedLoc.csv.gz")
rownames(chi2vals) <- sub(".","",make.names(chi2vals$IPI),fixed=TRUE)
rownames(chi2loc) <- sub(".","",make.names(chi2loc$IPI),fixed=TRUE)
loccols <- 6:14
chi2vals <- chi2vals[,loccols]
chi2loc <- chi2loc[,c(loccols,15)] ## also Notes col.
colnames(chi2vals) <- paste("chi2",
                            c("Mito-Fr4","ER-Fr11","Golgi-Fr13",
                              "ERGDV-Fr14","EE-Fr16","RE-TGN-Fr17",
                              "PM-Fr19","PS-Fr20","Cyto-Fr30"),
                            sep="-")

## some proteins are not present in our feature data
print(Venn(list(fdata=rownames(fdata),chi2vals=rownames(chi2vals))))
print(Venn(list(fdata=rownames(fdata),chi2loc=rownames(chi2loc))))
## let's keep only those for which we have quantitation
fdata <- cbind(fdata,chi2vals[rownames(fdata),])
fdata <- cbind(fdata,chi2loc[rownames(fdata),])
fdata <- new("AnnotatedDataFrame",data=fdata)

.experiment <- new("MIAPE",
                   lab="Center for Experimental BioInformatics (CEBI)",
                   samples = list(
                       species = "Mus musculus",
                       tissue = "Liver"), 
                   title="A mammalian organelle map by protein correlation profiling.",
                   abstract="Protein localization to membrane-enclosed organelles is a central feature of cellular organization. Using protein correlation profiling, we have mapped 1,404 proteins to ten subcellular locations in mouse liver, and these correspond with enzymatic assays, marker protein profiles, and confocal microscopy. These localizations allowed assessment of the specificity in published organellar proteomic inventories and demonstrate multiple locations for 39% of all organellar proteins. Integration of proteomic and genomic data enabled us to identify networks of coexpressed genes, cis-regulatory motifs, and putative transcriptional regulators involved in organelle biogenesis. Our analysis ties biochemistry, cell biology, and genomics into a common framework for organelle analysis.",
                   pubMedIds="16615899",
                   other = list(
                       MS = "LF",
                       spatexp = "PCP",
                       markers.fcol = "markers",
                       prediction.fcol = NA))

.process <- new("MSnProcess",
                processing=c(
                  paste("Loaded on ",date(),".",sep="")),
                normalised=TRUE,
                files=c("PIIS0092867406003692.mmc4-highDensity.csv.gz",
                  "PIIS0092867406003692.mmc4-lowDensity.csv.gz"))

foster2006 <- new("MSnSet",
                  exprs=eset,
                  featureData=fdata,
                  phenoData=pdata,
                  experimentData=.experiment,
                  processingData=.process)

stopifnot(pRolocdata:::valid.pRolocmetadata(pRolocmetadata(foster2006)))

if (validObject(foster2006))
  save(foster2006,file="../../data/foster2006.RData",
       compress = "xz", compression_level = 9)

