#!/usr/bin/env Rscript
#Written by Malachi Griffith

args = (commandArgs(TRUE))
datadir = args[1];                    #used to load input files and dump results files
sample_count = as.numeric(args[2])    #number of samples considered
max_genes = as.numeric(args[3])       #max number of genes/rows that will be allowed in the heat map

if (length(args) < 3){
  message_text1 = "Required arguments missing: ./AllEvents.pm.R /tmp/converge_all_events 8 50"
  stop(message_text1)
}

#datadir = "/tmp/converge_all_events"
#sample_count = 8
#max_genes = 50

#Load required libraries
library(plotrix)

#Load in the matrix data
setwd(datadir)
legend_data = read.table("events_legend.txt", header=TRUE, sep="\t", as.is=c(1:4), comment.char = "")
sample_data_all = read.table("events_final.tsv", header=TRUE, sep="\t", as.is=c(1:(sample_count+2)), na.strings=c("-"))
sample_data_labels = read.table("events_final_labels.tsv", header=TRUE, sep="\t", as.is=c(1:(sample_count+2)), na.strings=c("-"))
sample_data_numerical = read.table("events_final_numerical.tsv", header=TRUE, sep="\t", as.is=c(1:2))

#Reorder all data objects according the total number of samples affected
o = order(sample_data_all[,"grand_subject_count"], decreasing=TRUE)
sample_data_all = sample_data_all[o,]
sample_data_labels = sample_data_labels[o,]
sample_data_numerical = sample_data_numerical[o,]

#If there are more than $max_genes row present, limit to only the top $max_genes
if (dim(sample_data_all)[1] > max_genes){
  sample_data_all = sample_data_all[1:max_genes,]
  sample_data_labels = sample_data_labels[1:max_genes,]
  sample_data_numerical = sample_data_numerical[1:max_genes,]
}

#Get the sample names
sample_names = names(sample_data_all)[3:((3+sample_count)-1)]


########### Plot Version 1 ###########
#Get the events matrix
#Create a color matrix by mapping each event index to a color from the legend file
events = sample_data_numerical[,sample_names]
events_color_matrix = events
for (i in 1:length(events[1,])){
  for (j in 1:length(legend_data[,1])){
    x = which(events[,i] == legend_data[j,"index"])
    if (length(x) > 0){
      events_color_matrix[x,i] = legend_data[j,"color"]
    }
  }
}

#Create the heatmap and legend as a two page PDF
pdf(file="heatmap_ordered_by_recurrence.pdf")
par(mar=c(5, 5.5, 4, 2) +0.1) #c(bottom, left, top, right) - make some extra room for long gene names
main_title = "Events by sample by gene (ordered by recurrence)"

if (max_genes > 100){
  color2D.matplot(events, cellcolors=as.matrix(events_color_matrix), axes=FALSE, xlab=NA, ylab="Genes", show.legend=FALSE, show.values=FALSE, vcex=0.4, main=main_title, border=NA)
  axis(side=1, at=c(0.5:(dim(events)[2]-0.5)), labels=sample_names, las=2, cex.axis=0.75)
}else{
  color2D.matplot(events, cellcolors=as.matrix(events_color_matrix), axes=FALSE, xlab=NA, ylab=NA, show.legend=FALSE, show.values=FALSE, vcex=0.4, main=main_title)
  axis(side=1, at=c(0.5:(dim(events)[2]-0.5)), labels=sample_names, las=2, cex.axis=0.75)
  axis(side=2, at=c((dim(events)[1]-0.5):0.5), labels=sample_data_all[,"ensg_name"], las=2, cex.axis=0.65, font=3)
}

plot(x=1:10, y=1:10, type="n", xaxt="n", yaxt="n", ylab=NA, xlab=NA)
color.legend(5,2,6,8, legend_data[,"event_type"], rect.col=legend_data[,"color"], align="rb", gradient="y")
dev.off()

########### Plot Version 2 ########### 
#Recreate the heatmap ordered by gene name
o = order(sample_data_all[,"ensg_name"], decreasing=FALSE)
sample_data_all = sample_data_all[o,]
sample_data_labels = sample_data_labels[o,]
sample_data_numerical = sample_data_numerical[o,]

events = sample_data_numerical[,sample_names]
events_color_matrix = events
for (i in 1:length(events[1,])){
  for (j in 1:length(legend_data[,1])){
    x = which(events[,i] == legend_data[j,"index"])
    if (length(x) > 0){
      events_color_matrix[x,i] = legend_data[j,"color"]
    }
  }
}

pdf(file="heatmap_ordered_by_gene_name.pdf")
par(mar=c(5, 5.5, 4, 2) +0.1) #c(bottom, left, top, right) - make some extra room for long gene names
main_title = "Events by sample by gene (ordered by gene name)"

if (max_genes > 100){
  color2D.matplot(events, cellcolors=as.matrix(events_color_matrix), axes=FALSE, xlab=NA, ylab="Genes", show.legend=FALSE, show.values=FALSE, vcex=0.4, main=main_title, border=NA)
  axis(side=1, at=c(0.5:(dim(events)[2]-0.5)), labels=sample_names, las=2, cex.axis=0.75)
}else{
  color2D.matplot(events, cellcolors=as.matrix(events_color_matrix), axes=FALSE, xlab=NA, ylab=NA, show.legend=FALSE, show.values=FALSE, vcex=0.4, main=main_title)
  axis(side=1, at=c(0.5:(dim(events)[2]-0.5)), labels=sample_names, las=2, cex.axis=0.75)
  axis(side=2, at=c((dim(events)[1]-0.5):0.5), labels=sample_data_all[,"ensg_name"], las=2, cex.axis=0.65, font=3)
}

plot(x=1:10, y=1:10, type="n", xaxt="n", yaxt="n", ylab=NA, xlab=NA)
color.legend(5,2,6,8, legend_data[,"event_type"], rect.col=legend_data[,"color"], align="rb", gradient="y")
dev.off()


#Dump out the list of genes actually plotted
write.table(sample_data_labels, file="genes_displayed.tsv", quote = FALSE, sep="\t", row.names=FALSE)


########### Plot Version 3 ###########
#Recreate the heatmap ordered by gene name and by patient mutation occurence
o = order(sample_data_all[,"grand_subject_count"], decreasing=TRUE)
sample_data_all = sample_data_all[o,]
sample_data_labels = sample_data_labels[o,]
sample_data_numerical = sample_data_numerical[o,]

events = sample_data_numerical[,sample_names]
events_color_matrix = events
for (i in 1:length(events[1,])){
  for (j in 1:length(legend_data[,1])){
    x = which(events[,i] == legend_data[j,"index"])
    if (length(x) > 0){
      events_color_matrix[x,i] = legend_data[j,"color"]
    }
  }
}

#Convert all numbers to 1 to allow hierarchical sorting by patient columns
sample_data_numerical2 = sample_data_numerical
for (i in 3:length(sample_data_numerical2[1,])){
  j = which(sample_data_numerical2[,i] > 0)
  if (length(j) > 0){
    sample_data_numerical2[j,i] = 1
  }
}

dd = sample_data_numerical2
dt = t(dd)
dts = dt[do.call(order, c(lapply(3:NCOL(dt), function(i) dt[, i]), decreasing=TRUE)), ] 
dds = t(dts)
sample_names_sorted = colnames(dds)[3:NCOL(dds)]

pdf(file="heatmap_ordered_by_gene_name_and_patient.pdf")
par(mar=c(5, 5.5, 4, 2) +0.1) #c(bottom, left, top, right) - make some extra room for long gene names
main_title = "Events by sample by gene (hierarchical ordering by recurrence)"

if (max_genes > 100){  
  color2D.matplot(events[,sample_names_sorted], cellcolors=as.matrix(events_color_matrix[,sample_names_sorted]), axes=FALSE, xlab=NA, ylab="Genes", show.legend=FALSE, show.values=FALSE, vcex=0.4, main=main_title, border=NA)
  axis(side=1, at=c(0.5:(dim(events)[2]-0.5)), labels=sample_names_sorted, las=2, cex.axis=0.75)
}else{
  color2D.matplot(events[,sample_names_sorted], cellcolors=as.matrix(events_color_matrix[,sample_names_sorted]), axes=FALSE, xlab=NA, ylab=NA, show.legend=FALSE, show.values=FALSE, vcex=0.4, main=main_title)
  axis(side=1, at=c(0.5:(dim(events)[2]-0.5)), labels=sample_names_sorted, las=2, cex.axis=0.75)
  axis(side=2, at=c((dim(events)[1]-0.5):0.5), labels=sample_data_all[,"ensg_name"], las=2, cex.axis=0.65, font=3)
}

plot(x=1:10, y=1:10, type="n", xaxt="n", yaxt="n", ylab=NA, xlab=NA)
color.legend(5,2,6,8, legend_data[,"event_type"], rect.col=legend_data[,"color"], align="rb", gradient="y")
dev.off()




