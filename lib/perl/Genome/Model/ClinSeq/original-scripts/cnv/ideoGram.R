#!/usr/bin/env Rscript
#Written by Malachi Griffith

#Load user specified arguments
args = (commandArgs(TRUE))
ideogram_file = args[1];          #Path to a file containing chromosome ideogram values
target_chr = args[2];
chr_start = args[3];
chr_end = args[4];

if (length(args) < 4){
  message_text = "Required arguments missing for ideoGram.R - (e.g. ideoGram.R  /gsc/scripts/opt/genome/db/tgi/misc-annotation/human/build37-20130113/centromere.csv chr1 0 0) [ideogram data, chromosome, start, end]"
  stop(message_text)
  message(message_text)
}

#Example input
#ideogram_file = "/gscmnt/sata132/techd/mgriffit/reference_annotations/hg19/ideogram/ChrBandIdeogram.tsv"
#target_chr = "chr1"
#chr_start = 0
#chr_end = 0

#Load input data
ideo_data = read.table(ideogram_file, sep="\t", header=FALSE, comment.char="#", as.is=c(1,4,5))
names(ideo_data) = c("chrom","chromStart", "chromEnd", "name", "gieStain")

drawIdeogram = function(ideo_data, target_chr){

  main_title = paste("Ideogram for ", target_chr, sep = "")

  #Get the ideogram value for the current chromosome
  i = which(ideo_data[,"chrom"] == target_chr)
  ideo = ideo_data[i,]

  #Assign colors to each ideogram feature according to giemsa staining value
  ideo[,"color"] = "white"
  ii = which(ideo[,"gieStain"] == "gneg"); if (length(ii) > 0){ideo[ii, "color"] = "white"}
  ii = which(ideo[,"gieStain"] == "gpos25"); if (length(ii) > 0){ideo[ii, "color"] = "grey75"}
  ii = which(ideo[,"gieStain"] == "gpos50");  if (length(ii) > 0){ideo[ii, "color"] = "grey50"}
  ii = which(ideo[,"gieStain"] == "gpos75"); if (length(ii) > 0){ideo[ii, "color"] = "grey25"}
  ii = which(ideo[,"gieStain"] == "gpos100"); if (length(ii) > 0){ideo[ii, "color"] = "black"}
  ii = which(ideo[,"gieStain"] == "gvar"); if (length(ii) > 0){ideo[ii, "color"] = "firebrick1"}
  ii = which(ideo[,"gieStain"] == "stalk"); if (length(ii) > 0){ideo[ii, "color"] = "firebrick2"}
  ii = which(ideo[,"gieStain"] == "acen"); if (length(ii) > 0){ideo[ii, "color"] = "firebrick"}

  #Define some y coords
  ideo_bottom = 1
  ideo_middle = 1.5
  ideo_top = 2
  ideo_text = 2.95
  ideo_text_cex = 0.75

  #Get the midpoints of each ideogram feature
  ideo[,"chromMid"] = ideo[,"chromStart"] + ((ideo[,"chromEnd"] - ideo[,"chromStart"])/2)

  #Identify the non-centromeric ideogram features
  nc = which(!ideo[,"gieStain"] == "acen")

  #Create a blank plot.  If a chr range was specified, limit the ylim accordingly
  x = c(ideo[,"chromStart"], ideo[,"chromEnd"])
  y = rep(1, length(x))
  y_lower = 0
  y_upper = 3.75
  x_lower = min(x)
  x_upper = max(x)
  if (chr_start > 0){
    x_lower = chr_start
    x_upper = chr_end
  }
  plot(x=x, y=y, col="black", xlim=c(x_lower, x_upper), ylim=c(y_lower, y_upper), main=main_title, xlab="", ylab="", xaxt="n", yaxt="n", type="n", bty="n")

  icount_nc = length(ideo[nc,"chromStart"])
  rect(xleft=ideo[nc,"chromStart"], ybottom=rep(ideo_bottom, icount_nc), xright=ideo[nc,"chromEnd"], ytop=rep(ideo_top, icount_nc), col=ideo[nc,"color"])

  #Draw centromeres as triangles (using polygon)
  #Assume there are only two centromere features per chromosome
  cn = which(ideo[,"gieStain"] == "acen")
  polygon(x=c(ideo[cn[1],"chromStart"], ideo[cn[1],"chromStart"], ideo[cn[1],"chromEnd"]), y=c(ideo_bottom, ideo_top, ideo_middle), col=ideo[cn[1],"color"])
  polygon(x=c(ideo[cn[2],"chromStart"], ideo[cn[2],"chromEnd"], ideo[cn[2],"chromEnd"]), y=c(ideo_middle, ideo_top, ideo_bottom), col=ideo[cn[2],"color"])

  #Add labels for each ideogram feature
  text(x=ideo[,"chromMid"], y=ideo_text, labels=ideo[,"name"], srt=45, cex=ideo_text_cex)
}

filename=paste(target_chr, ".ideogram.jpeg", sep="")
jpeg(filename=filename, quality=100, width = 1500, height = 300)
drawIdeogram (ideo_data, target_chr)
dev.off()

