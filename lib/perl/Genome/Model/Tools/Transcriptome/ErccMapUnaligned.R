#!/usr/bin/env Rscript

library(ggplot2)
library(optparse)

option_parser <- OptionParser(option_list=list(
    make_option(c("-f", "--filename")),
    make_option(c("--epsilon"), type="double", default=0.001)
))

options  <- parse_args(option_parser)
filename <- options[[1]]
epsilon  <- options[[2]]


data = read.delim(filename)

data$count = as.numeric(data$count)

data$log_concentration = log2(data$concentration)
data$log_count         = log2(data$count + 1)


pdf('concentration-histogram.pdf')
ggplot(data, aes(x=log_concentration)) + geom_histogram()
dev.off()

pdf('count-histogram.pdf')
ggplot(data, aes(x=log_count)) + geom_histogram()
dev.off()

count_model <- lm(log_count ~ log_concentration, data=data)
count_r_squared = summary(count_model)[['r.squared']]

pdf('count-vs-concentration.pdf')
ggplot(data, aes(x=log_concentration, y=log_count)
       ) + geom_point(aes(color=subgroup)
       ) + geom_smooth(method=lm
       ) + annotate('text', 5, -3,
                    label=paste("R^2 =", count_r_squared, sep=' '))
dev.off()
