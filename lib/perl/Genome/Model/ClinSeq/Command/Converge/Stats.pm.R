#!/usr/bin/env Rscript

usage <- function() {
  print("Rscript Stats.pm.R /tmp/all.stats.tsv /tmp/all.plots.pdf")
}

plot_stats <- function(stats_f, plot_f) {
  stats<-read.table(stats_f, head = T)
  #Assume first six cols are "Question  Data_Type Analysis_Type Statistic_Type  Extra_Description File_Source"
  #followed by the sample data.
  stats.short<-stats[, c(1, 7:ncol(stats))]
  pdf(plot_f)
  for (i in 1:nrow(stats.short)) {
    stats.short.m <- melt(stats.short[i, ], id.vars = "Question");
    #skip cases with all NA's
    if(!any(complete.cases(stats.short.m))) {
      next
    }
    print(ggplot(na.omit(stats.short.m)) + geom_bar(aes(x = variable, y = value)) +
      xlab("Sample") + ylab(stats.short.m$Question) + opts(axis.text.x = theme_text(angle = 90)))
  }
  dev.off()
}

main <- function()
{
  args <- commandArgs(trailingOnly = TRUE)
  if(length(args) != 2) {
    usage()
    quit()
  }
  library(ggplot2)
  stats_f = args[1]
  plot_f = args[2]
  plot_stats(stats_f, plot_f)
}

main()
