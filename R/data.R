pRolocdata <- function() 
  data(package = "pRolocdata")


##' Extracts relevant metadata from an \code{MSnSet} instance. See
##' \code{README.md} for a description and explanation of the metadata
##' fields.
##'
##' @title Extract pRoloc metadata
##' @param x A \code{pRolocdata} data.
##' @return An instance of class \code{pRolocmetadata}.
##' @author Laurent Gatto
##' @aliases print.pRolocmetadata
##' @examples
##' library("pRolocdata")
##' data(dunkley2006)
##' data(dunkley2006)
##' pRolocmetadata(dunkley2006)
pRolocmetadata <- function(x) {
    ans <- list(Species = experimentData(x)@samples$species,
                Tissue = experimentData(x)@samples$tissue,
                CellLine = ifelse(experimentData(x)@samples$tissue == "Cell", 
                    experimentData(x)@samples$cellLine, NA),
                PMID = pubMedIds(x),
                MS = otherInfo(experimentData(x))$MS,
                Experiment = otherInfo(experimentData(x))$spatexp,
                MarkerCol = otherInfo(experimentData(x))$markers.fcol,
                PredictionCol = otherInfo(experimentData(x))$prediction.fcol)
    class(ans) <- c("list", "pRolocmetadata")
    ans
}

print.pRolocmetadata <- function(x, ...) {
    cat("pRoloc experiment metadata:\n")
    nx <- names(x)
    for (i in nx)
        cat(paste0(" ", i, ": ", x[[i]], "\n"))
}

valid.pRolocmetadata <- function(x) {
    stopifnot(inherits(x, "pRolocmetadata"))
    !any(sapply(x, is.null))
}

##' @title Create a SpatialMaps Account
##' @param x A \code{pRolocdata} data.
##' @return An instance of class \code{pRolocmetadata}.
##' @aliases print.pRolocmetadata
##' @examples
##' library("pRolocdata")
createAccount <- function(password = "prompt", email = "prompt") {
  projectAPI <- "AIzaSyC7iVp_D4iCAOl1e6ymW9TB7aC9E8tbjD4"
  email <- readline(prompt = "Enter Email: ")
  AuthUrl <- paste0("https://www.googleapis.com/identitytoolkit/v3/relyingparty/signupNewUser?key=", projectAPI)
  userData <- POST(url = AuthUrl, body = list("email" = email, "password" = password), encode = "json")
  return(content(userData))
  print("Your SpatialMaps Account was successfully created!")
  print("We send you a mail - Please verify your email.")
}

##' @title Reset the SpatialMaps account password
##' @param x A \code{pRolocdata} data.
##' @return An instance of class \code{pRolocmetadata}.
##' @aliases print.pRolocmetadata
##' @examples
##' library("pRolocdata")
resetPassword <- function(email){
  projectAPI <- "AIzaSyC7iVp_D4iCAOl1e6ymW9TB7aC9E8tbjD4"
  AuthUrl <- paste0("https://www.googleapis.com/identitytoolkit/v3/relyingparty/getOobConfirmationCode?key=", projectAPI)
  userData <- POST(url = AuthUrl, body = list("email" = email, "requestType" = "PASSWORD_RESET"), encode = "json")
  if ("error" %in% names(content(userData))) {
    warning(paste0("User email ", email, " was not found in the database"))
  } else {
    print(paste0("Password reset email was send to ", email))
  }
}

##' @title Login & Create Securitytoken
##' @param x A \code{pRolocdata} data.
##' @return An instance of class \code{pRolocmetadata}.
##' @aliases print.pRolocmetadata
##' @examples
##' library("pRolocdata")
login <- function(pgp) {
  if (password == "prompt" && email == "prompt"){
    email <- readline(prompt = "Email: ")
    password <- readline(prompt = "Password: ")
    print("Connecting to SpatialMaps:")
  } else if(password != "prompt" && email != "prompt"){
    print("Connecting to SpatialMaps:")
  } else {
    warning("please provide your email and password")
  }
  
  projectAPI = "AIzaSyC7iVp_D4iCAOl1e6ymW9TB7aC9E8tbjD4"
  AuthUrl <- paste0("https://www.googleapis.com/identitytoolkit/v3/relyingparty/verifyPassword?key=", projectAPI)
  userData <- httr::POST(url = AuthUrl, body = list("email" = email, "password" = password), encode = "json")
  return(content(userData))
}

##' @title Download Datasets from SpatialMaps
##' @param x A \code{pRolocdata} data.
##' @return An instance of class \code{pRolocmetadata}.
##' @aliases print.pRolocmetadata
##' @examples
##' library("pRolocdata")
download <- function(dataset, randomKey, pgp, password="none") {
    dbURL <- "https://spatialmap-1b08e.firebaseio.com"
    path <- paste0("/objects/", dataset)
    #retrieving data
    data <- GET(paste0(dbURL,path,".json"))
    retrievedData <- httr::content(data,"text")
    tempPath2 <- tempfile()
    writeBin(base64_dec(fromJSON(retrievedData)), tempPath2)
    x <- readRDS(tempPath2)
    assign(toString(as.name(dataset)), x, envir = .GlobalEnv)
    return(paste0(dataset, " was transfered"))
}

##' @title Upload Datasets from SpatialMaps
##' @param x A \code{pRolocdata} data.
##' @return An instance of class \code{pRolocmetadata}.
##' @aliases print.pRolocmetadata
##' @examples
##' library("pRolocdata")
upload <- function(dataset, name){
  dbURL <- "https://spatialmap-1b08e.firebaseio.com"
  #pRolocMetaData
  pRolocMeta <- pRolocMetaFrame(eval(as.name(dataset)), varName = name)
  Response <- POST(paste0(dbURL,"/meta",".json"), body = toJSON(pRolocMeta, auto_unbox = TRUE))
  #pRolocRawData
  pRolocRaw <- pRolocRawData(eval(as.name(dataset)))
  PUT(paste0(dbURL,"/raw/",content(Response),".json"), body = toJSON(pRolocRaw, auto_unbox = TRUE))
  #pRolocData
  pRolocDataVar <- pRolocFData(eval(as.name(dataset)))
  PUT(paste0(dbURL,"/data/",content(Response),".json"), body = toJSON(pRolocDataVar, auto_unbox = TRUE))
  #success message
  print(paste0(name, " got transfered to firebase."))
}

##' @title Upload Datasets from SpatialMaps
##' @param x A \code{pRolocdata} data.
##' @return An instance of class \code{pRolocmetadata}.
##' @aliases print.pRolocmetadata
##' @examples
##' library("pRolocdata")
createColors <- function(object){
  markers <- fData(object)$markers
  uniqueMarkers <- unique(markers)
  markerVec <- c()
  for(i in 1:length(uniqueMarkers)){
    markerColor <- ifelse(uniqueMarkers[i] == "unknown", getStockcol(), getStockcol()[i])
    markerVec <- c(markerVec, markerColor)
  }
  colorTable <- data.frame(uniqueMarkers, markerVec, stringsAsFactors = FALSE)
  
  colorAssigment <- unlist(sapply(markers, function(x) colorTable$markerVec[which(colorTable$uniqueMarkers == x)]))
  return(colorAssigment)
}

##' @title Upload Datasets from SpatialMaps
##' @param x A \code{pRolocdata} data.
##' @return An instance of class \code{pRolocmetadata}.
##' @aliases print.pRolocmetadata
##' @examples
##' library("pRolocdata")
#extract data from MSnSet object
pRolocRawData <- function(object){
  #convert object to base64
  tempPath <- tempfile()
  saveRDS(object, file = tempPath)
  binarySet <- readBin(tempPath, what = "raw", n = 50000000)
  base64Set <- jsonlite::base64_enc(binarySet)
  #adding key by assigning to data.frame
  pRolocList <- list("base64Set" =  base64Set)
  return(pRolocList)
}

##' @title Upload Datasets from SpatialMaps
##' @param x A \code{pRolocdata} data.
##' @return An instance of class \code{pRolocmetadata}.
##' @aliases print.pRolocmetadata
##' @examples
##' library("pRolocdata")
pRolocFData <- function(object){
  pcaData <- as.data.frame(plot2D(object, plot = FALSE))
  
  fScatter <- data.frame("PCA1" = pcaData[[1]], 
                        "PCA2" = pcaData[[2]], 
                        "Colors" = createColors(object))
  fSetData <- fData(object)
  
  for (i in 1:length((fSetData))){
    if (i == 1){
      p <- data.frame(fSetData[[i]])
    } else {
      p <- data.frame(p, fSetData[[i]])
    }
  }
  
  #filtering forbidden keys
  originalNames <- names(fSetData)
  originalNames <- gsub("\\$","-", originalNames)
  originalNames <- gsub("\\#","-", originalNames)
  originalNames <- gsub("\\]","-", originalNames)
  originalNames <- gsub("\\[","-", originalNames)
  originalNames <- gsub("\\/","-", originalNames)
  originalNames <- gsub("\\.","-", originalNames)
  names(p) <- originalNames
  
  p <- cbind(p, data.frame("id" = row.names(fSetData)))
  fSet <- cbind(fScatter,p)
  
  exprsSet <- exprs(object)
  exprsSet <- cbind(exprsSet, data.frame("id" = row.names(exprsSet)))
  row.names(exprsSet) <- NULL
  
  #filtering forbidden keys
  originalNames2 <- names(exprsSet)
  originalNames2 <- gsub("\\$","-", originalNames2)
  originalNames2 <- gsub("\\#","-", originalNames2)
  originalNames2 <- gsub("\\]","-", originalNames2)
  originalNames2 <- gsub("\\[","-", originalNames2)
  originalNames2 <- gsub("\\/","-", originalNames2)
  originalNames2 <- gsub("\\.","-", originalNames2)
  names(exprsSet) <- originalNames2
  
  pRolocList <- list("fSet" = fSet, "exprsSet" = exprsSet)
  return(pRolocList)
}

##' @title Upload Datasets from SpatialMaps
##' @param x A \code{pRolocdata} data.
##' @return An instance of class \code{pRolocmetadata}.
##' @aliases print.pRolocmetadata
##' @examples
##' library("pRolocdata")
pRolocMetaFrame <- function(object, varName){
  #meta
  #varName <- "varName"
  title <-  object@experimentData@title
  author <- object@experimentData@name
  email <- object@experimentData@email
  contact <- object@experimentData@contact
  dataStamp <- object@experimentData@dateStamp
  abstract <- object@experimentData@abstract
  lab <- object@experimentData@lab
  pubMedIds <- object@experimentData@pubMedIds
  
  tissue <- object@experimentData@samples$tissue
  cellLine <- object@experimentData@samples$cellLine
  species <- object@experimentData@samples$species
  operator <- object@experimentData@samples$operator
  
  markerClasses <- toString(pRoloc::getMarkerClasses(object))
  featureNames <- toString(featureNames(object))
  
  #List generation
  pRolocList <- list("varName" = varName, 
                    "title" = title,
                    "author" = author, 
                    "email" = email, 
                    "contact" = contact, 
                    "dataStamp" = dataStamp, 
                    "abstract" = abstract, 
                    "lab" = lab, 
                    "pubMedIds" = pubMedIds,
                    "tissue" = tissue,
                    "cellLine" = cellLine,
                    "species" = species,
                    "operator" = operator,
                    "markerClasses" = markerClasses,
                    "featureNames" = featureNames
  )
  
  return(pRolocList)
}

##' @title Upload Datasets from SpatialMaps
##' @param x A \code{pRolocdata} data.
##' @return An instance of class \code{pRolocmetadata}.
##' @aliases print.pRolocmetadata
##' @examples
##' library("pRolocdata")
update <- function(dataset, name, email, password, randomKey) {
    dbURL <- "https://spatialmap-1b08e.firebaseio.com"
    #pRolocMetaData
    pRolocMeta <- pRolocMetaFrame(eval(as.name(dataset)), varName = name)
    PUT(paste0(dbURL,"/meta/",randomKey,".json"), body = toJSON(pRolocMeta, auto_unbox = TRUE))
    #pRolocRawData
    pRolocRaw <- pRolocRawData(eval(as.name(dataset)))
    PUT(paste0(dbURL,"/raw/",randomKey,".json"), body = toJSON(pRolocRaw, auto_unbox = TRUE))
    #pRolocData
    pRolocDataVar <- pRolocFData(eval(as.name(dataset)))
    PUT(paste0(dbURL,"/data/",randomKey,".json"), body = toJSON(pRolocDataVar, auto_unbox = TRUE))
    #success message
    print(paste0(name, " got transfered to firebase."))
}

##' @title SpatialMap database backup function 
##' @param outputFile The output file name that should include the .json identifier {string}. File is written to the working directory. 
##' @param secretKey The database secret key {string}. Can be optained from options -> project config -> accounts -> database secrets.
##' @return An instance of class \code{pRolocmetadata}.
dataBackup <- function(outputFile="SpatialMaps.json", secretKeyPGP = NULL){
    if (is.null(secretKeyPGP)) secretKeyPGP <- readline(prompt = "Enter SecretKey: ")
    print("Fetching Data")
    dbURL <- "https://spatialmap-1b08e.firebaseio.com"
    urlPath = paste0(dbURL, "/.json?auth=", secretKeyPGP)
    curl_download(url = urlPath,
                  destfile = outputFile,
                  quiet = FALSE)
    print(paste0("Backup file created: ", outputFile))
}

