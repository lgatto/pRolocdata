data(mulvey2015)
data(mulvey2015norm)

library("SummarizedExperiment")
mulvey2015_se <- as(mulvey2015, "SummarizedExperiment")
mulvey2015norm_se <- as(mulvey2015norm, "SummarizedExperiment")

save(mulvey2015_se, file = "../../data/mulvey2015_se.rda", compress = "xz", compression_level = 9)
save(mulvey2015norm_se, file = "../../data/mulvey2015norm_se.rda", compress = "xz", compression_level = 9)
