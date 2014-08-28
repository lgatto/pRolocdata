## Function to add Uniprot Tan IDs to the Tan datasets. See README for mnore details.
addTanIds <- function(object) {
  ## Load datasets from FlyMine
  ## csvF = the filtered dataset from FlyMine that contains *only* Uniprot Ids that are reviewed
  ## csvU = the unfiltered dataset from Flymine that conatins *all* Uniprot Ids per CG number
  ##        this includes multiple IDs per CG both unreviewed and reviewed entries
  csvF <- read.csv("../extdata/TanFlyMineFiltered.csv")
  csvU <- read.csv("../extdata/TanFlyMineUnfiltered.csv")
  
  ## Remove any duplicated Uniprot IDs in filtered set
  .rm <- which(duplicated(csvF[,5]))
  csvF <- csvF[-.rm, ]
  
  ## Remove any duplicated Uniprot IDs in unfiltered set
  .rm <- which(duplicated(csvU[,5]))
  csvU <- csvU[-.rm, ]
  
  ## Select unique from unfiltered
  .a <- which(isUnique(csvU[,1]))
  csvUa <- csvU[.a, ]
  csvUb <- csvU[-.a, ]
  
  ## Now select non-unique and concatenate these IDs per match e.g. if we have say
  ## CG10077 which corresponds to IDs Q7KU78, Q8MRJ6, Q8MZI3 and Q8T0I3 we wish to
  ## create concatenate these to one character "Q7KU78;Q8MRJ6;Q8MZI3;Q8T0I3"
  .b <- unique(csvUb[,1])
  nu <- data.frame(.b, vector("character", length(.b)), vector("character", length(.b)),
                   vector("character", length(.b)), vector("character", length(.b)))
  colnames(nu) <- c("CGnumberAll", "AcessionNoAll", "EntryName", "AcessionNo", "EntryName")
  nu[,2] <- as.character(nu[,2])
  nu[,3] <- as.character(nu[,3])
  nu[,4] <- as.character(nu[,4])
  nu[,5] <- as.character(nu[,5])
  .n <- sapply(.b, function(z) which(csvUb[,1]==z))
  for (i in 1:length(.n)) {
    if (length(.n[[i]])==2) {
      nu[i,2] <- paste(csvUb[.n[[i]][1], 4], csvUb[.n[[i]][2], 4], sep=";")
      nu[i,3] <- paste(csvUb[.n[[i]][1], 5], csvUb[.n[[i]][2], 5], sep=";") 
    } else {
      a <- paste(csvUb[.n[[i]][1], 4], csvUb[.n[[i]][2], 4], sep=";")
      b <- paste(csvUb[.n[[i]][1], 5], csvUb[.n[[i]][2], 5], sep=";")
      for (j in 3:length(.n[[i]])) {
        a <- paste(a, csvUb[.n[[i]][j], 4], sep=";")
        b <- paste(b, csvUb[.n[[i]][j], 5], sep=";")
      }
      nu[i,2] <- a
      nu[i,3] <- b
    }
    nu[i,4] <- as.character(csvUb[.n[[i]][1], 4])
    nu[i,5] <- as.character(csvUb[.n[[i]][1], 5])
  }
  nn <- data.frame(cbind(as.character(csvUa[,1]), as.character(csvUa[,4]), 
                         as.character(csvUa[,5]), as.character(csvUa[,4]), 
                         as.character(csvUa[,5])))
  colnames(nn) <- colnames(nu)
  
  ## Put together data.frame of all IDs
  allids <- rbind(nu, nn)
  .match <- sapply(csvF[,1], function(z) which(allids[,1]==as.character(z)))
  for (i in 1:length(.match)) {
    allids[.match[i], 4] <- as.character(csvF[i, 4])
    allids[.match[i], 5] <- as.character(csvF[i, 5])
  }
  
  ## Now add ids to MSnSet
  fData(object)$AccessionNo <- vector("character", nrow(object))
  fData(object)$EntryName <- vector("character", nrow(object))
  fData(object)$AccessionNoAll <- vector("character", nrow(object))
  fData(object)$EntryNameAll <- vector("character", nrow(object))
  .cg <- allids[,1]
  .a <- sapply(.cg, function(z) which(fData(object)$Protein.ID==z))
  for (i in 1:length(.a)) {
    if(length(.a[[i]]) > 0) {
      fData(object)$AccessionNo[.a[[i]]] <- as.character(allids[i, 4])
      fData(object)$EntryName[.a[[i]]] <- as.character(allids[i, 5])
      fData(object)$AccessionNoAll[.a[[i]]] <- as.character(allids[i, 2])
      fData(object)$EntryNameAll[.a[[i]]] <- as.character(allids[i, 3])
    }
  }
  missing <- which(fData(object)$AccessionNo=="")
  for (i in 1:length(missing)) {
    fData(object)$AccessionNo[missing[i]] <- paste("NO_ID_", i, sep="")
  }
  return(object)
}

