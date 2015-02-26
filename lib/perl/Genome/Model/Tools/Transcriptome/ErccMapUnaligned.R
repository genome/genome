#!/usr/bin/env Rscript

library(plyr, quietly=TRUE)
library(reshape, warn.conflicts=FALSE, quiet=TRUE)
library(ggplot2, quietly=TRUE)
library(getopt)

opt <- getopt(
matrix(
  c('filename', 'f', 2, "character",
  'epsilon' ,   'e', 2, "double",
  'help'    ,   'h', 0, "logical"),
  byrow=TRUE, ncol=4 )
);

#help was asked for.
if ( !is.null(opt$help) ) {
  #get the script name (only works when invoked with Rscript).
  self <- commandArgs()[1];
  #print a friendly message and exit with a non-zero error code
  cat(paste("Usage: ",self," [-[-help|h]] [-[-filename|f] <output.pdf>] [-[-epsilon|e] <0.01>]\n",sep=""));
  q(status=1);
}

#set some reasonable defaults for the options
 if ( is.null(opt$filename ) ) { opt$filename <- 'output.pdf' }
 if ( is.null(opt$epsilon  ) ) { opt$epsilon <- 0.01          }

filename <- opt$filename
epsilon  <- opt$epsion

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

# signal success and exit
q(status=0)
