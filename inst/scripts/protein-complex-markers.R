tmp <- read.delim("../extdata/Protein_list.txt", header=TRUE,
                  as.is = TRUE)
n <- sapply(MSnbase:::utils.ssv2list(tmp$AC), length)

pcmrkdf <- data.frame(complexe = rep(tmp[, 1], time = n),
                      accession = MSnbase:::utils.ssv2vec(tmp$AC))
table(duplicated(pcmrkdf[, 2]))
pcmrkdf <- pcmrkdf[!duplicated(pcmrkdf[, 2]), ]

pcmrk <- as.character(pcmrkdf[, 1])
names(pcmrk) <- pcmrkdf[, 2]

save(pcmrk, file = "../extdata/protein-complex-markers.rda")
