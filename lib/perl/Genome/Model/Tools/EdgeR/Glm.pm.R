#!/usr/bin/env Rscript

library(edgeR)
library(getopt)

glmAnalysis <- function(inputFile, groups, subjects, outputFile, pvalue) {
    counts <- read.table(inputFile, head=T, row.names=1)
    dge <- DGEList(counts=counts, genes=rownames(counts))
    dge <- calcNormFactors(dge)
    design <- model.matrix(~subjects+groups)
    rownames(design) <- colnames(dge)
    dge <- estimateGLMCommonDisp(dge, design)
    dge <- estimateGLMTrendedDisp(dge, design)
    dge <- estimateGLMTagwiseDisp(dge, design)
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit)
    de <- decideTestsDGE(lrt, p=pvalue)
    top <- topTags(lrt, n=length(lrt$table$PValue))
    out <- cbind(top$table, de)

    colnames(out)[length(out)] <- "test.result"
    out <- out[order(out$PValue), ]

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
    "pvalue", "p", 1, "double",
    "subjects", "s", 1, "character"
    ), byrow=TRUE, ncol=4)

opts <- parseOptions(optSpec)

groups <- factor(unlist(strsplit(opts$groups, split=",")))
subjects <- factor(unlist(strsplit(opts$subjects, split=",")))
dge <- glmAnalysis(opts$`input-file`, groups, subjects, opts$`output-file`, opts$`pvalue`)
