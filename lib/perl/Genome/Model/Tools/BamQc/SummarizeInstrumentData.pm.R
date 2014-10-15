#!/usr/bin/env Rscript

library(plyr)
library(gplots)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
    stop("Give input file, output pdf")
}

input_file <- args[1]
output_pdf <- args[2]

load_isize_data <- function(path) {
    x <- read.table(path, head=TRUE, skip=10)
    colnames(x) <- c("insert_size", "count")
    x$count <- x$count / sum(x$count)
    x
}

load_sample <- function(name, data) {
    ds <- data[data$sample_name == name,]
    rv <- NULL
    for (i in 1:nrow(ds)) {
        path <- Sys.glob(sprintf("%s/*insert_size_metrics", ds$bamqc_path[i]))
        lib <- ds$library_name[i]
        id <- ds$id[i]
        x <- load_isize_data(path)
        x$lib <- lib
        x$id <- id
        x$flow_cell_id <- ds$flow_cell_id[i]
        x$lane <- ds$lane[i]
        rv <- rbind(rv, x)
    }
    rv$lib <- factor(rv$lib)
    o <- order(rv$lib, rv$id)
    rv$id <- factor(rv$id, levels=unique(rv$id[o]))
    rv$flow_cell_id <- paste(rv$flow_cell_id, rv$lane)
    rv
}

isize_plot <- function(df) {
    ggplot(df) +
        geom_bar(mapping=aes(x=insert_size, y=count, color=lib, fill=lib),
            stat="identity", position="stack") +
        labs(x="Insert size", y="Density") +
        facet_grid(id ~ ., scales="free_y") +
        scale_color_brewer(palette="Set1", name="Library") +
        scale_fill_brewer(palette="Set1", name="Library") +
        opts(panel.background=theme_blank(), panel.grid=theme_blank(),
            legend.position="top", legend.direction="vertical")
}

data <- read.table(input_file,
                   head=TRUE, stringsAsFactors=FALSE)
data$duplicate_percent <- as.numeric(as.character(data$duplicate_percent))

dup.smp <- ddply(data, .(sample_name), function (x) {
        dr <- sum(x$clusters * x$duplicate_percent) / sum(x$clusters)
        data.frame(duplication=dr, common_name=unique(x$common_name), ndata=nrow(x))
    })

dup.lib <- ddply(data, .(library_name), function (x) {
        dr <- sum(x$clusters * x$duplicate_percent) / sum(x$clusters)
        data.frame(duplication=dr, common_name=unique(x$common_name), ndata=nrow(x))
    })

dup.smp$sample_name <- sprintf("%s (%d idata)", dup.smp$sample_name, dup.smp$ndata)
dup.lib$library_name <- sprintf("%s (%d idata)", dup.lib$library_name, dup.lib$ndata)

print(dup.smp)
print(dup.lib)

pdf(output_pdf, width=11, height=8)

dup.smp$sample_name <- factor(
    dup.smp$sample_name,
    levels=dup.smp$sample_name[order(dup.smp$duplication)]
    )

dup.lib$library_name <- factor(
    dup.lib$library_name,
    levels=dup.lib$library_name[order(dup.lib$duplication)]
    )

maxdup <- min(100, max(data$duplicate_percent + 10))

text_hjust <- 1.2

ggplot(dup.smp, aes(sample_name, duplication)) +
    geom_bar(mapping=aes(color=common_name), stat="identity", position="stack") +
    coord_flip(ylim=c(0, maxdup)) +
    geom_text(aes(label=sprintf("%.2f%%", duplication), color=common_name), hjust=text_hjust) +
    opts(panel.background=theme_blank(), panel.grid=theme_blank()) +
    opts(title="Duplication by Sample") +
    labs(x="Duplicate percent", y="Sample name") +
    scale_color_brewer(palette="Set2")

ggplot(dup.lib, aes(library_name, duplication, color=common_name)) +
    geom_bar(stat="identity", position="stack") +
    coord_flip(ylim=c(0, maxdup)) +
    geom_text(aes(label=sprintf("%.2f%%", duplication)), hjust=text_hjust) +
    opts(panel.background=theme_blank(), panel.grid=theme_blank()) +
    opts(title="Duplication by Library") +
    labs(x="Duplicate percent", y="Sample name") +
    scale_color_brewer(palette="Set2")

ggplot(data, aes(duplicate_percent, fill=common_name)) +
    geom_histogram() +
    facet_grid(common_name ~ .) +
    labs(x="Duplicate percent") +
    opts(legend.position="top", legend.direction="vertical", title="Duplication by Tissue Type") +
    scale_fill_brewer(palette="Set2")

ggplot(data, aes(clusters, duplicate_percent, color=common_name)) +
    geom_point() +
    labs(x="# Clusters", y="Duplicate percent") +
    opts(title="Duplication by # Clusters") +
    opts(legend.position="top", legend.direction="vertical") +
    scale_color_brewer(palette="Set2")

samples <- unique(data$sample_name)
for (s in samples) {
    cat(sprintf("Processing sample %s...\n", s))
    ds <- data[data$sample_name == s,]
    ds <- ds[order(ds$library_name, ds$flow_cell_id, ds$lane),]
    libs <- unique(ds$library_name)
    df <- load_sample(s, data)
    title <- sprintf("Insert Size: Sample %s (%s)", s, unique(ds$common_name))
    p <- isize_plot(df) + opts(title=title)
    tbl <- ds[, c("library_name", "id", "flow_cell_id", "lane", "clusters", "duplicate_percent")]
    textplot(tbl, show.rownames=FALSE)
    title(title)
    print(p)
}

dev.off()
