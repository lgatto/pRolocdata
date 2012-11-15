.onAttach <- function(libname, pkgname) {
  msg <- paste0("\nThis is pRolocdata version ",
                packageVersion("pRolocdata"), ".\n",
                "Use 'pRolocdata()' to list available data sets.")
  packageStartupMessage(msg)  
}
