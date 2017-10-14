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


##' @title SpatialMap Account Creation
##' @description Uploads require a user login. The createAccount function allows to create the account within R with valid password and email. In case a spatialMaps account was already created on the platform, no new account is neccesssary. 
##' @param email Any email that can receive the confirmation mail. 
##' @param password Any password with a minimum of complexity. 
##' @return Message of success or failure.
##' @export
createAccount <- function(password = "prompt", email = "prompt") {
  projectAPI <- fromJSON("keys.json")$projectAPI
  email <- readline(prompt = "Enter Email: ")
  AuthUrl <- paste0("https://www.googleapis.com/identitytoolkit/v3/relyingparty/signupNewUser?key=", projectAPI)
  userData <- POST(url = AuthUrl, body = list("email" = email, "password" = password), encode = "json")
  return(httr::content(userData))
  print("Your SpatialMaps Account was successfully created!")
  print("We send you a mail - Please verify your email.")
}

##' @title spatialMaps account password reset
##' @description The password reset function in case the password was lost. The function triggers a password reset email sent to the input email.
##' @param email The SpatialMaps account email. 
##' @return Returns success or failure warning.
##' @export
resetPassword <- function(email){
  projectAPI <- fromJSON("keys.json")$projectAPI
  AuthUrl <- paste0("https://www.googleapis.com/identitytoolkit/v3/relyingparty/getOobConfirmationCode?key=", projectAPI)
  userData <- POST(url = AuthUrl, body = list("email" = email, "requestType" = "PASSWORD_RESET"), encode = "json")
  if ("error" %in% names(httr::content(userData))) {
    warning(paste0("User email ", email, " was not found in the database"))
  } else {
    print(paste0("Password reset email was send to ", email))
  }
}

##' @title Login & Create Security token.
##' @description The login function provides all information to interact with the SpatialMaps database. Its main functioanlity is to receive the JWT token required to upload MSnSets.
##' @param password The SpatialMaps account password.
##' @param email The SpatialMaps account email
##' @return returns list with various account information
##' @export 
login <- function(email = "prompt", password = "prompt", simple = TRUE){
  if (password == "prompt" && email == "prompt") {
    email <- readline(prompt = "Email: ")
    password <- readline(prompt = "Password: ")
    print("Connecting to SpatialMaps:")
  } else if (password != "prompt" && email != "prompt") {
    print("Connecting to SpatialMaps:")
  } else {
    warning("please provide your email and password")
  }
  
  projectAPI = jsonlite::fromJSON("keys.json")$projectAPI
  AuthUrl <- paste0("https://www.googleapis.com/identitytoolkit/v3/relyingparty/verifyPassword?key=", projectAPI)
  userData <- httr::POST(url = AuthUrl, body = list("email" = email, "password" = password), encode = "json")
  if (simple) { 
    return(httr::content(userData)$idToken)
  } else {
    return(httr::content(userData))
  }
}

##' @title Download Datasets from SpatialMaps.
##' @description SpatialMaps download function. The dataset name can be gathered from the SpatialMaps platform.
##' @param dataset The shortlink received from SpatialMaps
##' @return returns the requested SpatialMaps s4 class. 
##' @export
download <- function(dataset) {
    dbURL <- jsonlite::fromJSON("keys.json")$dbURL
    path <- paste0("/raw/", dataset)
    #retrieving data
    data <- httr::GET(paste0(dbURL,path,".json"))
    retrievedData <- httr::content(data,"text")
    tempPath2 <- tempfile()
    writeBin(base64_dec(fromJSON(retrievedData)$base64Set), tempPath2)
    return(readRDS(tempPath2))
}

##' @title Upload Datasets from SpatialMaps
##' @param dataset The proloc object.
##' @param name A string to add the name of the dataset.
##' @param token The JWT token created with the login() function.
##' @return Reports the successfull transfer and outputs the random ID.
##' @export
upload <- function(dataset, name, token, public = TRUE){
  dbURL <- fromJSON("keys.json")$dbURL
  #pRolocMetaData
  pRolocMeta <- pRolocMetaFrame(eval(as.name(dataset)), varName = name, token = token, public = public)
  Response <- POST(paste0(dbURL,"/meta",".json"), body = toJSON(pRolocMeta, auto_unbox = TRUE))
  #pRolocRawData
  pRolocRaw <- pRolocRawData(eval(as.name(dataset)))
  PUT(paste0(dbURL,"/raw/",httr::content(Response),".json"), body = toJSON(pRolocRaw, auto_unbox = TRUE))
  #pRolocData
  pRolocDataVar <- pRolocFData(eval(as.name(dataset)))
  PUT(paste0(dbURL,"/data/",httr::content(Response),".json"), body = toJSON(pRolocDataVar, auto_unbox = TRUE))
  #success message
  print(paste0(name, " got transfered to firebase."))
}


##' @title Token verification
##' @description Verifies the authenticity of the generated token and returns the uniquer user identifier.
##' @param token The JWT token created with the login() function.
##' @return Reports the successfull transfer and outputs the random ID.
tokenCheck <- function(token){
  projectAPI <- fromJSON("keys.json")$projectAPI
  dbURL <- dbURL <- fromJSON("keys.json")$dbURL
  requestURL <- paste0("https://www.googleapis.com/identitytoolkit/v3/relyingparty/getAccountInfo?key=", projectAPI)
  requestData = httr::POST(url = requestURL,  body = list("idToken" = token), encode = "json")
  return(httr::content(requestData)$users[[1]]$localId)
}

##' @title Append colors to the fSet table.
##' @description A helper function to pass pRoloc conform colors to the fData table.
##' @param object The pRoloc object.
##' @return Appends the color column to the fSet.
createColors <- function(object){
  markers <- fData(object)$markers
  uniqueMarkers <- unique(markers)
  markerVec <- c()
  for (i in 1:length(uniqueMarkers)) {
    markerColor <- ifelse(uniqueMarkers[i] == "unknown", pRoloc::getStockcol(), pRoloc::getStockcol()[i])
    markerVec <- c(markerVec, markerColor)
  }
  colorTable <- data.frame(uniqueMarkers, markerVec, stringsAsFactors = FALSE)
  
  colorAssigment <- unlist(sapply(markers, function(x) colorTable$markerVec[which(colorTable$uniqueMarkers == x)]))
  return(colorAssigment)
}

##' @title s4 class to base64 binary
##' @description Transforms MSnSet object to binary base64 string. This step is required to keep store the s4 class without interfering with the JSON formatting requirements of firebase and to prevent any column reordering or other modifications.
##' @param object The MsnSet object.
##' @return A base64 encrypted MSnSet object list.
pRolocRawData <- function(object){
  tempPath <- tempfile()
  saveRDS(object, file = tempPath)
  binarySet <- readBin(tempPath, what = "raw", n = 50000000)
  base64Set <- jsonlite::base64_enc(binarySet)
  #Adding the database key by assigning the base64 string to a list
  pRolocList <- list("base64Set" =  base64Set)
  return(pRolocList)
}

##' @title fData list aggregation.
##' @description Helper function to reformat fData to DB suitable format.
##' @param object The MSnSet.
pRolocFData <- function(object){
  pcaData = as.data.frame(pRoloc::plot2D(object, plot = FALSE))
  
  fScatter = data.frame("PCA1" = pcaData[[1]], 
                        "PCA2" = pcaData[[2]], 
                        "Colors" = createColors(object))
  fSetData = fData(object)
  
  for (i in 1:length((fSetData))) {
    if (i == 1) {
      p = data.frame(fSetData[[i]])
    } else {
      p = data.frame(p, fSetData[[i]])
    }
  }
  
  #filtering forbidden keys
  names(p) <- filterKeys(names(fSetData))

  exprsSet = exprs(object)
  exprsSet = cbind(exprsSet, data.frame("id" = row.names(exprsSet)))
  row.names(exprsSet) = NULL
  
  names(exprsSet) <- filterKeys(names(exprsSet))
  fSet = cbind(fScatter,p,exprsSet)
  
  pRolocList = list("fSet" = fSet)
  return(pRolocList)
}


##' @title Filter forbidden json keys
##' @description The firebase realtime database probhibts certain synthax relevant column names. This function helps to overcome potential errors and warns the user if one or more column names needed to be renamed. 
##' @param x The column names of the dataset.
##' @return returns correctedNames and reports changes in the naming.
filterKeys <- function(x){
  modifiedNames <- x
  forbiddenKeys <- c("\\$","\\#","\\]","\\[","\\/","\\.")
  for (i in forbiddenKeys) modifiedNames <- gsub(i,"-", modifiedNames)
  if (!identical(x, modifiedNames)) 
    print("One or more column keys were renamed for the SpatialMaps usage. Note: Your raw dataset stays the same.")
  return(modifiedNames)
}

##' @title Metadata creation
##' @param object The MSnSet object.
##' @param varname The name of the object.
##' @param token The JWT token created with the login() function.
##' @description Creates the /Meta data entry for SpatialMaps
pRolocMetaFrame <- function(object, varName, token, public){
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
                    "public" = public,
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
##' @description Other than the upload function, the update function takes a specified random key as path (other than creating a new one). The main value of the update function is not to create duplicates in the database.
##' @param dataset A MSnSet object.
##' @param name Name of the MSnSet dataset.
##' @param token The JWT token created with the login() function.
##' @param randomKey The MSnSet random key. This key is provided either in on the 
##' SpatialMaps website or returned after each upload.
##' @return Returns success message. 
##' @export
update <- function(dataset, name, token, randomKey) {
    dbURL <- dbURL <- fromJSON("keys.json")$dbURL
    #pRolocMetaData
    pRolocMeta <- pRolocMetaFrame(eval(as.name(dataset)), varName = name, token = token)
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
##' @description Backup of all SpatialMaps database entries
##' @param outputFile The output file name that should include the .json identifier {string}. File is written into the working directory.
##' @param secret The database secret key. Can be optained from options -> project config -> accounts -> database secrets.
##' @return An instance of class \code{pRolocmetadata}.
##' @export
dataBackup <- function(outputFile="SpatialMaps.json", secret = NULL){
    if (is.null(secret)) secret <- readline(prompt = "Enter secret: ")
    print("Fetching Data")
    dbURL <- fromJSON("keys.json")$dbURL
    urlPath = paste0(dbURL, "/.json?auth=", secret)
    curl_download(url = urlPath,
                  destfile = outputFile,
                  quiet = FALSE)
    print(paste0("Backup file created: ", outputFile))
}