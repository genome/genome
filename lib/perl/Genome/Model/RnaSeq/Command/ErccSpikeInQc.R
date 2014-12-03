#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(TRUE)

filename <- args[1]
epsilon  <- 0.001

data = read.delim(filename)

data$FPKM  = as.numeric(data$FPKM)
data$count = as.numeric(data$count)

data$log_concentration = log2(data$concentration)
data$log_fpkm          = log2(data$FPKM + epsilon)
data$log_count         = log2(data$count + 1)


pdf('concentration-histogram.pdf')
ggplot(data, aes(x=log_concentration)) + geom_histogram()
dev.off()

pdf('fpkm-histogram.pdf')
ggplot(data, aes(x=log_fpkm)) + geom_histogram()
dev.off()

pdf('count-histogram.pdf')
ggplot(data, aes(x=log_count)) + geom_histogram()
dev.off()

fpkm_model <- lm(log_fpkm ~ log_concentration, data=data)
fpkm_r_squared = summary(fpkm_model)[['r.squared']]

pdf('fpkm-vs-concentration.pdf')
ggplot(data, aes(x=log_concentration, y=log_fpkm)
       ) + geom_point(aes(color=subgroup)
       ) + geom_smooth(method=lm
       ) + annotate('text', 5, -9,
                    label=paste("R^2 =", fpkm_r_squared, sep=' '))
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
