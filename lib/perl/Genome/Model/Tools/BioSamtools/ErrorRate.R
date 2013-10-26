# Read in a TSV file of Error Rate metrics

generatePlots <- function(inputFile,readEnd,countPlotFile,ratePlotFile,rateDistFile) {
  fullFile <- readTable(inputFile)
  readEndErrorRate <- getReadEnd(fullFile,readEnd)
  summaryErrorRate <- getSummaryErrorRate(readEndErrorRate)
  positionErrorRate <- getPositionErrorRate(readEndErrorRate)
  makeCountPlot(positionErrorRate,countPlotFile)
  makeRatePlot(positionErrorRate,ratePlotFile)
  makeRateDist(positionErrorRate,rateDistFile)
}

readTable <- function(filename) {
  fullFile <- read.delim(filename,header=TRUE,comment.char="#")
  return(fullFile)
}

getReadEnd <- function(fullFile,readEnd) {
  readEndErrorRate <- fullFile[fullFile$read_end == readEnd, ]
}

getSummaryErrorRate <- function(readEndErrorRate){
  # The last line is the SUMMARY
  summaryErrorRate <- readEndErrorRate[nrow(readEndErrorRate),]
  return(summaryErrorRate)
}

# Remove the last line to get all Positions
getPositionErrorRate <- function(readEndErrorRate) {
  positionErrorRate <- readEndErrorRate[-nrow(readEndErrorRate),]
  return(positionErrorRate)
}

makeCountPlot <- function(positionErrorRate,countPlotFile) {
  # These are the column names of the raw counts of the occurence of each type of base in comparison to the reference
  countColumnNames <- c("total","match","error","mismatch","ambiguous","insertion","deletion")

  # Subset the table to get only the above columns
  countErrorRate <- positionErrorRate[countColumnNames]

  # Assign colors to the columns of interest
  plotColors <- c("black","green","red","pink","brown","yellow","blue")

  # Open a device to print png image
  png(countPlotFile, bg="transparent", width=800, height=600)

  # Plot the first columna and each additional column
  plot(countErrorRate$total,type="l",log="y",col=plotColors[1],xlab="Position (bp)",ylab="log(Count)",main="Error Per Position")
  lines(countErrorRate$match,col=plotColors[2])
  lines(countErrorRate$error,col=plotColors[3])
  lines(countErrorRate$mismatch,col=plotColors[4])
  lines(countErrorRate$ambiguous,col=plotColors[5])
  lines(countErrorRate$insertion,col=plotColors[6])
  lines(countErrorRate$deletion,col=plotColors[7])
  legend("topright", names(countErrorRate),col=plotColors,cex=0.8,lty=1:3, lwd=2, bty="n");

  # Close the device we are printing the raw count plot to
  dev.off()
}

makeRatePlot <- function(positionErrorRate,ratePlotFile,readEnd=0) {
  # New color scheme for the ratio plot
  plotColors <- c("green","red","pink","brown","yellow","blue")

  # subset all the postion metrics to the ratio/rates
  ratioColumnNames = c("total","error_rate","mismatch_rate","ambiguous_rate","insertion_rate","deletion_rate")
  ratioErrorRate <- positionErrorRate[ratioColumnNames]

  # Start writing the positional error rate plot
  png(ratePlotFile, bg="transparent", width=800, height=600)
  par(mar = c(5, 4, 4, 4) + 0.3)
  plot(ratioErrorRate$error_rate,type="l",col=plotColors[2],xlab="Position (bp)",ylab="Error Rate",main=paste("Positional Error Ratio - Read",readEnd,sep=" "))
  lines(ratioErrorRate$mismatch_rate,col=plotColors[3])
  lines(ratioErrorRate$ambiguous_rate,col=plotColors[4])
  lines(ratioErrorRate$insertion_rate,col=plotColors[5])
  lines(ratioErrorRate$deletion_rate,col=plotColors[6])

  par(new=TRUE)
  plot(ratioErrorRate$total,type="l",col=plotColors[1],axes=FALSE,bty="n",xlab="",ylab="")
  axis(side=4,at=pretty(range(ratioErrorRate$total)))
  mtext("Total Sequence",side=4,line=3)
  legend("topright", names(ratioErrorRate),col=plotColors,cex=0.8,lty=1:3, lwd=2, bty="n");
  dev.off()
}

makeRateDist <- function(positionErrorRate,rateDistFile) {
  ratioColumnNames = c("total","error_rate","mismatch_rate","ambiguous_rate","insertion_rate","deletion_rate")
  ratioErrorRate <- positionErrorRate[ratioColumnNames]
  plotColors <- c("red","pink","brown","yellow","blue")
  png(rateDistFile, bg="transparent", width=800, height=600)
  boxPlotErrorRate <- ratioErrorRate[,-1]
  boxplot(boxPlotErrorRate,col=plotColors,xlab="Type",ylab="Ratio",main="Distribution of Error Ratio")
  dev.off()
}
  
