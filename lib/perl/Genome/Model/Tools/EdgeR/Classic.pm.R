#!/usr/bin/env Rscript

library(edgeR)
library(getopt)

classicAnalysis <- function(inputFile, groups, outputFile, pvalue) {
    counts <- read.table(inputFile, head=T, row.names=1)
    dge <- DGEList(counts=counts, group=groups)
    dge <- calcNormFactors(dge)
    dge <- estimateCommonDisp(dge)
    dge <- estimateTagwiseDisp(dge)
    et <- exactTest(dge)
    de <- decideTestsDGE(et, p=pvalue)
    out <- cbind(et$table, de)
    colnames(out)[4] <- "test.result"
    out <- out[order(out$p.value), ]

    write.table(out, outputFile, quote=FALSE, sep="\t",
        row.names=TRUE, col.names=NA)
    dge
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
    "groups", "G", 1, "character",
    "output-file", "o", 1, "character",
    "pvalue", "p", 1, "double"
    ), byrow=TRUE, ncol=4)

opts <- parseOptions(optSpec)

groups <- factor(unlist(strsplit(opts$groups, split=",")))
dge <- classicAnalysis(opts$`input-file`, groups, opts$`output-file`, opts$`pvalue`)
