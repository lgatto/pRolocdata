library("MSnbase")
library("pRoloc")
library("dplyr")

e <- 4:17

lopimsSyn2 <-
    readMSnSet2("../extdata/lopimsSyn2.csv",
                ecol = e, fnames = 1) %>%
    addMarkers("../extdata/lopims_markers.csv")
save(lopimsSyn2, file = "../../data/lopimsSyn2.rda")


lopimsSyn1 <-
    readMSnSet2("../extdata/lopimsSyn1.csv",
                ecol = e, fnames = 1) %>%
    addMarkers("../extdata/lopims_markers.csv")
save(lopimsSyn1, file = "../../data/lopimsSyn1.rda")

lopimsSyn2_0frags <-
    readMSnSet2("../extdata/lopimsSyn2_0frags.csv",
                ecol = e, fnames = 1) %>%
    addMarkers("../extdata/lopims_markers.csv")
save(lopimsSyn2_0frags, file = "../../data/lopimsSyn2_0frags.rda")
