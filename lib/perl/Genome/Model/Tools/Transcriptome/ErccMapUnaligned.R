#!/usr/bin/env Rscript

library(plyr, quietly=TRUE)
library(reshape, warn.conflicts=FALSE, quiet=TRUE)
library(ggplot2, quietly=TRUE)
library(getopt)

opt <- getopt(
  matrix(
    c('data'    ,   'd', 1, "character",
      'output'  ,   'o', 2, "character",
      'epsilon' ,   'e', 2, "double",
      'help'    ,   'h', 0, "logical"),
    byrow=TRUE, ncol=4
  )
);

#help was asked for.
if ( !is.null(opt$help) ) {
  #get the script name (only works when invoked with Rscript).
  self <- commandArgs()[1];
  #print a friendly message and exit with a non-zero error code
  cat(paste("Usage: ",self," [-[-help|h]] [-[-data|d] <data.tsv>] [-[-output|o] <output.pdf>] [-[-epsilon|e] <0.01>]\n",sep=""));
  q(status=1);
}

if ( is.null(opt$data) ) {
  #get the script name (only works when invoked with Rscript).
  cat(paste("Please specify an input data file via the '--data=<data.tsv>' option! \n",sep=""));
  q(status=1);
}

#set some reasonable defaults for the options
 if ( is.null(opt$output  ) ) { opt$output  <- 'output.pdf' }
 if ( is.null(opt$epsilon ) ) { opt$epsilon <- 0.01         }

inputTSV <- opt$data
epsilon  <- opt$epsion
pdfFile  <- opt$output

data = read.delim(inputTSV)

data$count = as.numeric(data$count)

data$log_concentration = log2(data$concentration)
data$log_count         = log2(data$count + 1)


pdf(pdfFile)
ggplot(data, aes(x=log_concentration)) + geom_histogram()

ggplot(data, aes(x=log_count)) + geom_histogram()

count_model <- lm(log_count ~ log_concentration, data=data)
count_r_squared = summary(count_model)[['r.squared']]

ggplot(data, aes(x=log_concentration, y=log_count)
       ) + geom_point(aes(color=subgroup)
       ) + geom_smooth(method=lm
       ) + annotate('text', 5, -3,
                    label=paste("R^2 =", count_r_squared, sep=' '))
dev.off()

r2 <- paste("R^2 :", count_r_squared, sep=" ")
cat(paste("\n", "===========", r2, "===========", "\n", sep="\n"))
cat(paste("See", pdfFile, "for the analysis plots", "\n", sep=" "))

# signal success and exit
q(status=0)
