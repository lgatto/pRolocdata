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
##' @param password The SpatialMaps account password.
##' @param email Any email that can receive the confirmation mail.  
##' @return Message of success or failure.
createAccount <- function(password = "prompt", email = "prompt") {
  projectAPI <- fromJSON("keys.json")$projectAPI
  email <- readline(prompt = "Enter Email: ")
  AuthUrl <- paste0("https://www.googleapis.com/identitytoolkit/v3/relyingparty/signupNewUser?key=", projectAPI)
  userData <- POST(url = AuthUrl, body = list("email" = email, "password" = password), encode = "json")
  return(content(userData))
  print("Your SpatialMaps Account was successfully created!")
  print("We send you a mail - Please verify your email.")
}

##' @title API key store
##' @param Key The name of the API key to be called.
##' @description An internal function to store.
##' @param email Any email that can receive the confirmation mail. 
##' @return Message of success or failure.
createAccount <- function(password = "prompt", email = "prompt") {
  projectAPI <- fromJSON("keys.json")$projectAPI
  email <- readline(prompt = "Enter Email: ")
  AuthUrl <- paste0("https://www.googleapis.com/identitytoolkit/v3/relyingparty/signupNewUser?key=", projectAPI)
  userData <- POST(url = AuthUrl, body = list("email" = email, "password" = password), encode = "json")
  return(content(userData))
  print("Your SpatialMaps Account was successfully created!")
  print("We send you a mail - Please verify your email.")
}

##' @title Reset the SpatialMaps account password.
##' @param email The SpatialMaps account email. 
##' @return Returns success or failure warning.
resetPassword <- function(email){
  projectAPI <- fromJSON("keys.json")$projectAPI
  AuthUrl <- paste0("https://www.googleapis.com/identitytoolkit/v3/relyingparty/getOobConfirmationCode?key=", projectAPI)
  userData <- POST(url = AuthUrl, body = list("email" = email, "requestType" = "PASSWORD_RESET"), encode = "json")
  if ("error" %in% names(content(userData))) {
    warning(paste0("User email ", email, " was not found in the database"))
  } else {
    print(paste0("Password reset email was send to ", email))
  }
}

##' @title Login & Create Security token.
##' @param x A \code{pRolocdata} data.
##' @return An instance of class \code{pRolocmetadata}.
##' @aliases print.pRolocmetadata
login <- function(password = "prompt", email = "prompt", pgp = "none"){
  if (password == "prompt" && email == "prompt") {
    email <- readline(prompt = "Email: ")
    password <- readline(prompt = "Password: ")
    print("Connecting to SpatialMaps:")
  } else if (password != "prompt" && email != "prompt") {
    print("Connecting to SpatialMaps:")
  } else {
    warning("please provide your email and password")
  }
  
  projectAPI = fromJSON("keys.json")$projectAPI
  AuthUrl <- paste0("https://www.googleapis.com/identitytoolkit/v3/relyingparty/verifyPassword?key=", projectAPI)
  userData <- httr::POST(url = AuthUrl, body = list("email" = email, "password" = password), encode = "json")
  return(content(userData))
}

##' @title Download Datasets from SpatialMaps.
##' @param x A \code{pRolocdata} data.
##' @return An instance of class \code{pRolocmetadata}.
download <- function(dataset, randomKey, pgp, password="none") {
    dbURL <- dbURL <- fromJSON("keys.json")$dbURL
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
##' @param dataset The proloc object.
##' @param name A string to add the name of the dataset.
##' @return Reports the successfull transfer and outputs the random ID.
upload <- function(dataset, name, token){
  dbURL <- fromJSON("keys.json")$dbURL
  #pRolocMetaData
  pRolocMeta <- pRolocMetaFrame(eval(as.name(dataset)), varName = name, token = token)
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
##' @param dataset The proloc object.
##' @param name A string to add the name of the dataset.
##' @return Reports the successfull transfer and outputs the random ID.
tokenCheck <- function(token){
  projectAPI <- fromJSON("keys.json")$projectAPI
  dbURL <- dbURL <- fromJSON("keys.json")$dbURL
  requestURL <- paste0("https://www.googleapis.com/identitytoolkit/v3/relyingparty/getAccountInfo?key=", projectAPI)
  requestData = httr::POST(url = requestURL,  body = list("idToken" = token), encode = "json")
  return(httr::content(requestData)$users[[1]]$localId)
}

##' @title Append colors to the fSet table.
##' @param object The pRoloc object.
##' @return Appends the color column to the fSet.
createColors <- function(object){
  markers <- fData(object)$markers
  uniqueMarkers <- unique(markers)
  markerVec <- c()
  for (i in 1:length(uniqueMarkers)) {
    markerColor <- ifelse(uniqueMarkers[i] == "unknown", getStockcol(), getStockcol()[i])
    markerVec <- c(markerVec, markerColor)
  }
  colorTable <- data.frame(uniqueMarkers, markerVec, stringsAsFactors = FALSE)
  
  colorAssigment <- unlist(sapply(markers, function(x) colorTable$markerVec[which(colorTable$uniqueMarkers == x)]))
  return(colorAssigment)
}

##' @title Transform MSnSet object to binary base64 string.
##' @param object The MsnSet object.
##' @description Converts the S4 object to a base64 string to prevent any column reordering or other modifications. 
##' @return A list containing the base64 encrypted MSnSet s4 object.
pRolocRawData <- function(object){
  tempPath <- tempfile()
  saveRDS(object, file = tempPath)
  binarySet <- readBin(tempPath, what = "raw", n = 50000000)
  base64Set <- jsonlite::base64_enc(binarySet)
  #Adding the database key by assigning the base64 string to a list
  pRolocList <- list("base64Set" =  base64Set)
  return(pRolocList)
}

##' @title Upload Datasets from SpatialMaps
##' @param x A \code{pRolocdata} data.
##' @return An instance of class \code{pRolocmetadata}.
pRolocFData <- function(object){
  pcaData <- as.data.frame(plot2D(object, plot = FALSE))
  fScatter <- data.frame("PCA1" = pcaData[[1]], 
                        "PCA2" = pcaData[[2]], 
                        "Colors" = createColors(object))
  fSetData <- fData(object)
  
  for (i in 1:length((fSetData))) {
    if (i == 1) {
      p <- data.frame(fSetData[[i]])
    } else {
      p <- data.frame(p, fSetData[[i]])
    }
  }
  
  #filtering forbidden key names
  names(p) <- filterKeys(names(fSetData))
  
  p <- cbind(p, data.frame("id" = row.names(fSetData)))
  fSet <- cbind(fScatter,p)
  exprsSet <- exprs(object)
  exprsSet <- cbind(exprsSet, data.frame("id" = row.names(exprsSet)))
  row.names(exprsSet) <- NULL
  
  #filtering forbidden key names
  names(exprsSet) <- filterKeys(names(exprsSet))
  
  pRolocList <- list("fSet" = fSet, "exprsSet" = exprsSet)
  return(pRolocList)
}
##' @title Filter forbidden json keys
##' @param x The names of the dataset, called with names(x)
##' @return returns correctedNames and reports changes in the naming
filterKeys <- function(x){
  modifiedNames <- x
  forbiddenKeys <- c("\\$","\\#","\\]","\\[","\\/","\\.")
  for (i in forbiddenKeys) modifiedNames <- gsub(i,"-", modifiedNames)
  if (!identical(x, modifiedNames)) 
    warning("One or more column keys were renamed for the SpatialMaps usage. Note: Your raw dataset stays the same.")
  return(modifiedNames)
}

##' @title Meta data creation
##' @param object The MSnSet object.
##' @param varname The name of the object.
##' @description Creates the /Meta data entry for SpatialMaps
##' @return An instance of class \code{pRolocmetadata}.
pRolocMetaFrame <- function(object, varName, token){
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
  UID <- tokenCheck(token = token)
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
                    "UID" = UID,
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

##' @title Update existing SpatialMaps dataset. 
##' @param dataset A MSnSet object.
##' @param name Name of the MSnSet dataset.
##' @param email Email of the user.
##' @param password Password of the user
##' @param randomKey The MSnSet random key. This key is provided either in on the 
##' SpatialMaps website or returned after each upload.
##' @return Returns success message. 
update <- function(dataset, name, email, password, randomKey) {
    dbURL <- dbURL <- fromJSON("keys.json")$dbURL
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
    dbURL <- fromJSON("keys.json")$dbURL
    urlPath = paste0(dbURL, "/.json?auth=", secretKeyPGP)
    curl_download(url = urlPath,
                  destfile = outputFile,
                  quiet = FALSE)
    print(paste0("Backup file created: ", outputFile))
}