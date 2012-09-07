.onAttach <- function(libname, pkgname) {
  x <- apply(pRolocdata(), 1, function(.d) paste(.d, collapse = " from "))
  n <- length(x)
  x <- paste0(" ", 1:n, ": ", x, "\n")  
  msg <- c(paste0("\nThis is pRolocdata version ",
                 packageVersion("pRolocdata"), ".\n"),
           "Available data sets:\n",
           x)           
  packageStartupMessage(msg)  
}
