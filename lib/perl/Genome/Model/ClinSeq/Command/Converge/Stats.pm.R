#!/usr/bin/env Rscript

plot_stats <- function(stats_f, plot_f) {
  stats<-read.table(stats_f, head = T)
  stats.short<-stats[, c(1, 7:ncol(stats))]
  pdf(plot_f)
  for (i in 1:nrow(stats.short)) {
    print(i);
    stats.short.m <- melt(stats.short[i, ]);
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
