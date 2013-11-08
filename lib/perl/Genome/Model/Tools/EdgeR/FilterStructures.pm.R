
#!/usr/bin/env Rscript

library(edgeR)
library(getopt)

filterStructures <- function(inputFile, outputFile, quant, min_count) {
    counts <- read.table(inputFile, head=T, row.names=1)
    q <- apply(counts, 1, quantile, quant)
    filtered <- counts[q >= min_count, ]
    write.table(filtered, outputFile, quote=FALSE, sep="\t",
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
    "output-file", "o", 1, "character",
    "quantile", "q", 1, "double",
    "min-count", "c", 1, "integer"
    ), byrow=TRUE, ncol=4)

opts <- parseOptions(optSpec)

filterStructures(
    opts$`input-file`,
    opts$`output-file`,
    opts$quantile,
    opts$`min-count`
    )
