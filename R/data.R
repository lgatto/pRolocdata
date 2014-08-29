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
                Type = otherInfo(experimentData(x))$type,
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
