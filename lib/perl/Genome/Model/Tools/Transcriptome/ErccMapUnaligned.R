#!/usr/bin/env Rscript

library(plyr, quietly=TRUE)
library(reshape, warn.conflicts=FALSE, quiet=TRUE)
library(ggplot2, quietly=TRUE)
library(getopt)

# allow output to be 100 columns wide
options(width=as.integer(100))

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

#inputTSV <- opt$data
#epsilon  <- opt$epsilon
#pdfFile  <- opt$output

#################################
make.LinearityPlot <- function(df) {
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
  p <- p + opts(title="log2(alignment counts) vs log2(Mix 1 concentrations)")
  print(p)

  return(count_r_squared)
}

make.MixFrame <- function(mixType, df) {
  mix <- aggregate(df[c(mixType, "count")],
            by=list(Group=df$subgroup),
            FUN=sum)

  epsilon <- 0.0001

  colnames(mix)[3] <- "AlignmentCounts"
  mix$Probability <- mix[,mixType] / sum(mix[,mixType])
  mix$ExpectedCounts <- mix$Probability * sum(mix$AlignmentCounts)

  # account for 0 count cases ( log2(0) is undefined )
  tmpCounts <- ifelse(mix$AlignmentCounts == 0, epsilon, mix$AlignmentCounts)
  mix$log2AlignmentCounts <- log2(tmpCounts)

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
  test.mix1 <- chisq.test(Mix1DF$log2AlignmentCounts, p=Mix1DF$Probability)

  if (test.mix1[["p.value"]] < p_value_threshold) {
    hypo.mix1 <- FALSE
  }

#  cat(paste("Running Chi-squared test for Mix 2 usage\n"))
  test.mix2 <- chisq.test(Mix2DF$log2AlignmentCounts, p=Mix2DF$Probability)

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

df <- readTSV(opt$data)

mix1DF <- make.MixFrame("Mix1", df)
mix2DF <- make.MixFrame("Mix2", df)

cat(paste("\n===> Mix1 Summary Stats <=== \n"))
mix1DF
cat(paste("\n===> Mix2 Summary Stats <=== \n"))
mix2DF

test.mixture(mix1DF, mix2DF)

pdf(opt$output)
r2 <- make.LinearityPlot(df)
mix.plot(mix1DF, "Mix1")
mix.plot(mix2DF, "Mix2")

sink(file="/dev/null")
dev.off()
sink()


r2_text <- paste("R^2 :", r2, sep=" ")
cat(paste("===>", r2_text, "<===", "\n", sep=" "))

cat(paste("\nSee", opt$output, "for the analysis plots", "\n", sep=" "))

# signal success and exit
q(status=0)

#################################

#data = read.delim(inputTSV)
#
#data$count = as.numeric(data$count)
#
#data$log_concentration = log2(data$concentration)
#data$log_count         = log2(data$count + 1)
#
#
#pdf(pdfFile)
#ggplot(data, aes(x=log_concentration)) + geom_histogram()
#
#ggplot(data, aes(x=log_count)) + geom_histogram()
#
#count_model <- lm(log_count ~ log_concentration, data=data)
#count_r_squared = summary(count_model)[['r.squared']]
#
#ggplot(data, aes(x=log_concentration, y=log_count)
#       ) + geom_point(aes(color=subgroup)
#       ) + geom_smooth(method=lm
#       ) + annotate('text', 5, -3,
#                    label=paste("R^2 =", count_r_squared, sep=' '))
#dev.off()
#
#r2 <- paste("R^2 :", count_r_squared, sep=" ")
#cat(paste("\n", "===========", r2, "===========", "\n", sep="\n"))
#cat(paste("See", pdfFile, "for the analysis plots", "\n", sep=" "))
#
## signal success and exit
#q(status=0)

#> aggregate(df["Mix1"], by=list(Group=df$subgroup), FUN=sum)
#  Group     Mix1
#1     A 41406.01
#2     B 20703.01
#3     C 20703.01
#4     D 20703.01
#> aggregate(df["Mix2"], by=list(Group=df$subgroup), FUN=sum)
#  Group     Mix2
#1     A 10351.50
#2     B 20703.01
#3     C 31054.51
#4     D 41406.01
#> 41406.01 + 20703.01 + 20703.01 + 20703.01
#[1] 103515
#> 10351.50 + 20703.01 + 31054.51 + 41406.01
#[1] 103515
#> class(aggregate(df["Mix2"], by=list(Group=df$subgroup), FUN=sum))
#[1] "data.frame"
#> mix1 <- aggregate(df["Mix1"], by=list(Group=df$subgroup), FUN=sum)
#> mix2 <- aggregate(df["Mix2"], by=list(Group=df$subgroup), FUN=sum)
#> mix1$Mix1 / sum(mix1$Mix1)
#[1] 0.4 0.2 0.2 0.2
#> mix2$Mix2 / sum(mix2$Mix2)
#[1] 0.1 0.2 0.3 0.4
#>
#> mix1$Probability <- mix1$Mix1 / sum(mix1$Mix1)
#> mix2$Probability <- mix2$Mix2 / sum(mix2$Mix2)
#> mix1
#  Group     Mix1 Probability
#1     A 41406.01         0.4
#2     B 20703.01         0.2
#3     C 20703.01         0.2
#4     D 20703.01         0.2
#> mix2
#  Group     Mix2 Probability
#1     A 10351.50         0.1
#2     B 20703.01         0.2
#3     C 31054.51         0.3
#4     D 41406.01         0.4
#> aggregate(df["count"], by=list(AlignmentCounts=df$subgroup), FUN=sum)
#  AlignmentCounts  count
#1               A 278722
#2               B 170081
#3               C 138932
#4               D 148281
#> aggregate(df["Mix2", "count"], by=list(Group=df$subgroup), FUN=sum)
#Error in aggregate.data.frame(as.data.frame(x), ...) :
#  arguments must have same length
#> aggregate(df[c("Mix2", "count")], by=list(Group=df$subgroup), FUN=sum)
#  Group     Mix2  count
#1     A 10351.50 278722
#2     B 20703.01 170081
#3     C 31054.51 138932
#4     D 41406.01 148281
#> aggregate(df[c("Mix1", "count")], by=list(Group=df$subgroup), FUN=sum)
#  Group     Mix1  count
#1     A 41406.01 278722
#2     B 20703.01 170081
#3     C 20703.01 138932
#4     D 20703.01 148281
#> mix1 <- aggregate(df[c("Mix1", "count")], by=list(Group=df$subgroup), FUN=sum)
#
#> mix2 <- aggregate(df[c("Mix2", "count")], by=list(Group=df$subgroup), FUN=sum)
#
#> mix1$Probability <- mix1$Mix1 / sum(mix1$Mix1)
#> mix2$Probability <- mix2$Mix2 / sum(mix2$Mix2)
#> mix1
#  Group     Mix1  count Probability
#1     A 41406.01 278722         0.4
#2     B 20703.01 170081         0.2
#3     C 20703.01 138932         0.2
#4     D 20703.01 148281         0.2
#> sum(mix1$count)
#[1] 736016
#> sum(df$count)
#[1] 736016
#> mix1$Probability * sum(mix1$count)
#[1] 294406.4 147203.2 147203.2 147203.2
#> sum(mix1$Probability * sum(mix1$count))
#[1] 736016
#> mix1$ExpectedCounts <- mix1$Probability * sum(mix1$count)
#> mix1
#  Group     Mix1  count Probability ExpectedCounts
#1     A 41406.01 278722         0.4       294406.4
#2     B 20703.01 170081         0.2       147203.2
#3     C 20703.01 138932         0.2       147203.2
#4     D 20703.01 148281         0.2       147203.2
#> mix1$ExpectedCounts <- mix2$Probability * sum(mix2$count)
#> mix1$ExpectedCounts <- mix1$Probability * sum(mix1$count)
#> mix2$ExpectedCounts <- mix2$Probability * sum(mix2$count)
#> mix1
#  Group     Mix1  count Probability ExpectedCounts
#1     A 41406.01 278722         0.4       294406.4
#2     B 20703.01 170081         0.2       147203.2
#3     C 20703.01 138932         0.2       147203.2
#4     D 20703.01 148281         0.2       147203.2
#> mix2
#  Group     Mix2  count Probability ExpectedCounts
#1     A 10351.50 278722         0.1        73601.6
#2     B 20703.01 170081         0.2       147203.2
#3     C 31054.51 138932         0.3       220804.8
#4     D 41406.01 148281         0.4       294406.4
#> chisq.test(mix1$count, p=mix1$Probability)
#
#        Chi-squared test for given probabilities
#
#data:  mix1$count
#X-squared = 4863.81, df = 3, p-value < 2.2e-16
#
#> chisq.test(mix2$count, p=mix2$Probability)
#
#        Chi-squared test for given probabilities
#
#data:  mix2$count
#X-squared = 678091.5, df = 3, p-value < 2.2e-16
#
#> chisq.test(mix2$count, p=mix2$Probability, correct=FALSE)
#
#        Chi-squared test for given probabilities
#
#data:  mix2$count
#X-squared = 678091.5, df = 3, p-value < 2.2e-16
#
#> chisq.test(mix1$count, p=mix1$Probability, correct=FALSE)
#
#        Chi-squared test for given probabilities
#
#data:  mix1$count
#X-squared = 4863.81, df = 3, p-value < 2.2e-16
#
#> ((278722-294406.4)^2/294406.4) + ((170081-147203.2)^2/147203.2) + ((138932-147203.2)^2/147203.2) + ((148281-147203.2)^2/147203.2)
#[1] 4863.81
#> mix1$count - mix1$ExpectedCounts
#[1] -15684.4  22877.8  -8271.2   1077.8
#> mix2$count - mix2$ExpectedCounts
#[1]  205120.4   22877.8  -81872.8 -146125.4
#> log(mix1$count) - log(mix1$ExpectedCounts)
#[1] -0.054746256  0.144460849 -0.057829340  0.007295177
#> log(mix2$count) - log(mix2$ExpectedCounts)
#[1]  1.3315481  0.1444608 -0.4632944 -0.6858520
#> log(mix1$count) - log(mix1$ExpectedCounts)
#> logMix1Counts <- log(mix1$count)
#> logMix2Counts <- log(mix2$count)
#> chisq.test(logMix1Counts, p=mix1$Probability)
#
#        Chi-squared test for given probabilities
#
#data:  logMix1Counts
#X-squared = 3.9819, df = 3, p-value = 0.2634
#
#> chisq.test(logMix2Counts, p=mix2$Probability)
#
#        Chi-squared test for given probabilities
#
#data:  logMix2Counts
#X-squared = 16.2073, df = 3, p-value = 0.001028
#
#Warning message:
#In chisq.test(logMix2Counts, p = mix2$Probability) :
#  Chi-squared approximation may be incorrect
#> logMix1Counts <- log2(mix1$count)
#> logMix2Counts <- log2(mix2$count)
#> chisq.test(logMix1Counts, p=mix1$Probability)
#
#        Chi-squared test for given probabilities
#
#data:  logMix1Counts
#X-squared = 5.7447, df = 3, p-value = 0.1247
#
#> chisq.test(logMix2Counts, p=mix2$Probability)
#
#        Chi-squared test for given probabilities
#
#data:  logMix2Counts
#X-squared = 23.3822, df = 3, p-value = 3.361e-05
#
#>
#> df3 <- melt(mix2)
#Using Group as id variables
#> p2 <- ggplot(df3[df3$variable == "count" | df3$variable == "ExpectedCounts",], aes(x=Group, y=value, fill=variable))
#> p2 <- p2 + geom_bar(position="dodge", color="black") + scale_fill_brewer(palette="Pastel1")
#> p2
