.onAttach <- function(libname, pkgname) {
  
  x <- apply(pRolocdata(), 1, function(.d) paste(.d, collapse = ": "))
  n <- length(x)
  x <- paste0(" ", 1:n, ": ", x, "\n")  
  msg <- c(paste("\nThis is pRolocdata version",
                 packageVersion("pRolocdata"), ".\n"),
           "Available data sets:\n",
           x)           
  packageStartupMessage(msg)  
}
