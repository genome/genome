#!/usr/bin/env Rscript

library(edgeR)
library(getopt)

filterCountsFile <- function(inputFile, countsPerMillion, numSamples, outputFile) {
    counts <- read.table(inputFile, head=T, row.names=1)
    keep <- rowSums(cpm(counts)>countsPerMillion) >= numSamples
    output <- counts[keep,]
    write.table(output, outputFile,  quote=FALSE, sep="\t",
        row.names=TRUE, col.names=NA)
}

parseOptions <- function(optSpec) {
    opts <- getopt(optSpec)

    if (!is.null(opts$help)) {
        cat(getopt(optSpec, usage=TRUE))
        q(status=0)
    }

    required <- optSpec[optSpec[,3] == "1", 1]
    for (req in required) {
        if (is.null(opts[[req]])) {
            cat(paste("Required option", req, "missing!\n"))
            q(status=1)
        }
    }

    opts
}

optSpec = matrix(c(
    "help", "h", 0, NA,
    "input-file", "i", 1, "character",
    "counts-per-million", "c", 1, "integer",
    "num-samples", "n", 1, "integer",
    "output-file", "o", 1, "character"
    ), byrow=TRUE, ncol=4)

opts <- parseOptions(optSpec)

dge <- filterCountsFile(opts$`input-file`, opts$`counts-per-million`, opts$`num-samples`, opts$`output-file`)
