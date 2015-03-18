#!/usr/bin/env Rscript

# L I B R A R I E S ###########################################################
library(plyr, quietly=TRUE)
library(reshape, warn.conflicts=FALSE, quiet=TRUE)
library(ggplot2, quietly=TRUE)
library(getopt)

# C O N F I G S ###############################################################
# allow output to be 100 columns wide
options(width=as.integer(100))

# F U N C T I O N S ###########################################################
process.opts <- function() {
  opts <- getopt(
    matrix(
      c('data'    ,   'd', 1, "character",
        'output'  ,   'o', 2, "character",
        'epsilon' ,   'e', 2, "double",
        'help'    ,   'h', 0, "logical"),
      byrow=TRUE, ncol=4
    )
  );

  #help was asked for.
  if ( !is.null(opts$help) ) {
    #get the script name (only works when invoked with Rscript).
    self <- commandArgs()[1];
    #print a friendly message and exit with a non-zero error code
    cat(paste("Usage: ",self," [-[-help|h]] [-[-data|d] <data.tsv>] [-[-output|o] <output.pdf>] [-[-epsilon|e] <0.01>]\n",sep=""));
    q(status=1);
  }

  if ( is.null(opts$data) ) {
    #get the script name (only works when invoked with Rscript).
    cat(paste("Please specify an input data file via the '--data=<data.tsv>' option! \n",sep=""));
    q(status=1);
  }

  #set some reasonable defaults for the options
  if ( is.null(opts$output  ) ) { opts$output  <- 'output.pdf' }
  if ( is.null(opts$epsilon ) ) { opts$epsilon <- 0.01         }

  #inputTSV <- opts$data
  #epsilon  <- opts$epsilon
  #pdfFile  <- opts$output
  return(opts)
}

linearity.plot <- function(df) {
  df$count <- as.numeric(df$count)

  df$logMix1  <- log2(df$Mix1)
  df$logMix2  <- log2(df$Mix2)
  df$logCount <- log2(df$count + 1)


  count_model <- lm(logCount ~ logMix1, data=df)
  count_r_squared = summary(count_model)[['r.squared']]
  correlation_label <- paste("R^2 =", count_r_squared, sep=' ')

  p <- ggplot(df, aes(x=logMix1, y=logCount)) +
          geom_point(aes(color=subgroup)) +
          geom_smooth(method=lm) +
          annotate('text', 5, -3, label=correlation_label)

  p <- p +
       opts(title="log2(alignment counts) vs log2(Mix 1 concentrations)")

  # adjust axis labels
  p <- p + ylab("log2(Count)") + xlab("log2(Mix1)")

  print(p)

  return(count_r_squared)
}

make.MixFrame <- function(mixType, df) {
  mix <- aggregate(df[c(mixType, "count")],
            by=list(Group=df$subgroup),
            FUN=sum)

  colnames(mix)[3] <- "AlignmentCounts"
  mix$ExpectedMixRatios <- mix[,mixType] / sum(mix[,mixType])
  mix$ExpectedCounts <- mix$ExpectedMixRatios * sum(mix$AlignmentCounts)
  colnames(mix)[2] <- paste("Concentration", mixType, sep="")

  return(mix)
}

mix.plot <- function(MixDF, MixType) {
  MeltDF <- suppressMessages(melt(MixDF))
  df <- MeltDF[MeltDF$variable == "AlignmentCounts" |
               MeltDF$variable == "ExpectedCounts", ]
  # Reset the number factors on subsetted data frame
  df$variable <- factor(df$variable)

  p <- ggplot(df, aes(x=Group, y=value, fill=variable))
  p <- p + geom_bar(position="dodge", color="black") 
  title <- paste("Observed & Expected Alignment Counts (",
                 MixType, ")",
                 sep="")

  # adjust plot title
  p <- p + opts(title=title)

  # adjust y axis label
  p <- p + ylab("Count")

  # adjust legend placement
  p <- p + opts(legend.position="top", legend.direction="horizontal")

  # suppress legend title
  p <- p + scale_fill_discrete(name="")

  print(p)
}

test.mixture <- function(Mix1DF, Mix2DF) {
  hypo.mix1 <- TRUE
  hypo.mix2 <- TRUE
  p_value_threshold <- 0.05

#  cat(paste("Running Chi-squared test for Mix 1 usage\n"))
  test.mix1 <- chisq.test(Mix1DF$AlignmentCounts/10000, p=Mix1DF$ExpectedMixRatios)
#  print(test.mix1)

  if (test.mix1[["p.value"]] < p_value_threshold) {
    hypo.mix1 <- FALSE
  }

#  cat(paste("Running Chi-squared test for Mix 2 usage\n"))
  test.mix2 <- chisq.test(Mix2DF$AlignmentCounts/10000, p=Mix2DF$ExpectedMixRatios)
#  print(test.mix2)

  if (test.mix2[["p.value"]] < p_value_threshold) {
    hypo.mix2 <- FALSE
  }

  if ( (hypo.mix1 == FALSE && hypo.mix2 == FALSE) ||
       (hypo.mix1 == TRUE && hypo.mix2 == TRUE) ) {
    cat(paste("==> Couldn't ascertain whether mix 1 or mix 2 was used! <==\n"))
    cat(paste("\n"))
    cat(paste("==== Chi-squared test for Mix 1 usage results ==== \n"))
    print(test.mix1)
    cat(paste("==== Chi-squared test for Mix 2 usage results ==== \n"))
    print(test.mix2)
    return
  }

  if (hypo.mix1 == TRUE) {
    cat(paste("\n===> Mix 1 was used <=== \n"))
  }

  if (hypo.mix2 == TRUE) {
    cat(paste("\n===> Mix 2 was used <=== \n"))
  }
}

readTSV <- function(inputTSV) {
  df <- read.table(inputTSV, header=TRUE)

  # let's rename the columns so they're easier to use in this script
  # these were the original column names in the TSV:
  # 1. 'Re-sort ID'
  # 2. 'ERCC ID'
  # 3. 'subgroup'
  # 4. 'Mix 1 concentration (attomoles/ul)'
  # 5. 'Mix 2 concentration (attomoles/ul)'
  # 6. 'label'
  # 7. 'count'

  newnames <- c(
      "ReSortID",
      "ERCCID",
      "subgroup",
      "Mix1",
      "Mix2",
      "label",
      "count"
  )

  colnames(df) <- newnames
  df$subgroup <- factor(df$subgroup)
  return(df)
}

main <- function() {
  opts <- process.opts()

  df <- readTSV(opts$data)

  mix1DF <- make.MixFrame("Mix1", df)
  mix2DF <- make.MixFrame("Mix2", df)

  cat(paste("\n===> Mix1 Summary Stats <=== \n"))
  print(mix1DF)
  cat(paste("\n===> Mix2 Summary Stats <=== \n"))
  print(mix2DF)

  if (sum(mix1DF$AlignmentCounts) > 0 &&
      sum(mix2DF$AlignmentCounts) > 0) {
    test.mixture(mix1DF, mix2DF)
  }
  else {
    cat(paste("===> No ERCC Spike-Ins Found! <===\n"))
    q(status=0)
  }

  pdf(opts$output)
  r2 <- linearity.plot(df)
  mix.plot(mix1DF, "Mix1")
  mix.plot(mix2DF, "Mix2")

  sink(file="/dev/null")
  dev.off()
  sink()


  r2_text <- paste("R^2 :", r2, sep=" ")
  cat(paste("===>", r2_text, "<===", "\n", sep=" "))

  cat(paste("\nSee", opts$output, "for the analysis plots", "\n", sep=" "))
}

# M A I N #####################################################################
main()

# signal success and exit
q(status=0)
