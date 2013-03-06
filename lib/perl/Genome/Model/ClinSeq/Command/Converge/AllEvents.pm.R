#!/usr/bin/env Rscript
#Written by Malachi Griffith

#args = (commandArgs(TRUE))
#datadir = args[1];        #used to load input files and dump results files
#sample_count = args[2]    #number of samples considered

#if (length(args) < 2){
#  message_text1 = "Required arguments missing: ./AllEvents.pm.R /tmp/converge_all_events 8"
#  stop(message_text1)
#}
#datadir = "/tmp/converge_all_events"
#sample_count = 8

#Load required libraries
library(plotrix)

#Load in the matrix data
setwd(datadir)
legend_data = read.table("events_legend.txt", header=TRUE, sep="\t", as.is=c(1:4), comment.char = "")
sample_data_all = read.table("events_final.tsv", header=TRUE, sep="\t", as.is=c(1:(sample_count+2)), na.strings=c("-"))
sample_data_labels = read.table("events_final_labels.tsv", header=TRUE, sep="\t", as.is=c(1:(sample_count+2)), na.strings=c("-"))
sample_data_numerical = read.table("events_final_numerical.tsv", header=TRUE, sep="\t", as.is=c(1:2))

#Get the sample names
sample_names = names(sample_data)[3:((3+sample_count)-1)]

#Get the events matrix
events = sample_data_numerical[,sample_names]

#Create a color matrix by mapping each event index to a color from the legend file
events_color_matrix = events
for (i in 1:length(events[1,])){
  for (j in 1:length(legend_data[,1])){
    x = which(events[,i] == legend_data[j,"index"])
    if (length(x) > 0){
      events_color_matrix[x,i] = legend_data[j,"color"]
    }
  }
}

#color2D.matplot(events, cellcolors=color_array, axes=FALSE, xlab=NA, ylab=NA, show.legend=FALSE, nslices=nmax, show.values=TRUE, vcex=0.4, main=title)

main_title = "Events by sample by gene"
color2D.matplot(events, cellcolors=as.matrix(events_color_matrix), axes=FALSE, xlab=NA, ylab=NA, show.legend=FALSE, show.values=FALSE, vcex=0.4, main=main_title)
axis(side=1, at=c(0.5:(dim(events)[2]-0.5)), labels=sample_names, las=2, cex.axis=0.5)
axis(side=2, at=c((dim(events)[1]-0.5):0.5), labels=sample_data_all[,"ensg_name"], las=2, cex.axis=0.65, font=3)

plot(x=1:10, y=1:10, type="n", xaxt="n", yaxt="n", ylab=NA, xlab=NA)
color.legend(5,2,6,8, legend_data[,"event_type"], rect.col=legend_data[,"color"], align="rb", gradient="y")











