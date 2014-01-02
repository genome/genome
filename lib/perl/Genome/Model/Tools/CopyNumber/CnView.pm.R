#!/usr/bin/env Rscript
#Written by Malachi Griffith

#This script is used in conjunction with 'intersectCnvGeneList.pl'
#This script produces graphs to visualize CNV results by chromosome or portion thereof
#The input CNV data is a 'cnvs.hq' file from the TGI somatic variation pipeline (i.e. CNV diff value between Tumor and Normal over 10 kb windows)
#The user supplies a category name for the types of gene names to be highlighted on the plot.  
#The user must also supply a file containing these gene names and associated mean CNV value over the region of each gene
#If the gene passes a default cutoff the name of the gene will be place on the plot at the appropriate location


#Load user specified arguments
args = (commandArgs(TRUE))
category_name = args[1];          #Name of gene group being considered (e.g. 'Cancer Gene Census')
cnv_file = args[2];               #Path to file containing raw CNV diff values from the somatic variation pipeline
gene_file = args[3];              #Path to file containing gene records and associated mean CNV diff values for each gene
segments_file = args[4];          #Path to file containing 
ideogram_file = args[5];          #Path to a file containing chromosome ideogram values
outdir = args[6];                 #Directory where output files will be generated
chr = args[7];                    #Chr to be processed or 'ALL' for all chromosomes
chr_start = as.numeric(args[8]);  #If a single Chr is specified, a sub-region of the Chr can also be plotted - start position for this zoom in
chr_end = as.numeric(args[9]);    #If a single Chr is specified, a sub-region of the Chr can also be plotted - end position for this zoom in
image_type = args[10];

if (length(args) < 10){
  message_text = "Required arguments missing for CNView.R"
  stop(message_text)
  message(message_text)
}

#Debug values
#category_name = "Cancer Genes" 
#cnv_file = "/gscmnt/ams1183/info/model_data/2877353501/build112474967/variants/cnvs.hq"
#gene_file = "/gscmnt/sata132/techd/mgriffit/braf_resistance/cnvs/CNView_CancerGenes/CNView_CancerGenes_chr22_900000-2100000.tsv"
#ideogram_file = "/gscmnt/sata132/techd/mgriffit/reference_annotations/hg19/ideogram/ChrBandIdeogram.tsv"
#outdir = "/gscmnt/sata132/techd/mgriffit/braf_resistance/cnvs/CNView_CancerGenes/"
#chr = "chr22"
#chr_start = 900000
#chr_end = 2100000
#image_type = "none"

#FUNCTION SECTION

#############################################################################################################################################
#plotChrCNV - Function to plot CNV values by position and annotate with gene names and other relevant info                                  #
#############################################################################################################################################
plotChrCNV = function(target_chr, direction){
  #Create plotting space
  lay = layout(matrix(c(1,2,3,4), ncol=2, byrow=F), heights=c(0.75,2), widths=c(1,0.3), FALSE) #layout.show(lay)

  #Plot 1 - Draw the ideogram at the top
  #Adjust margins -> c(bottom, left, top, right)
  par(mar=c(0, 4, 1.5, 0))
  drawIdeogram(ideo_data, target_chr)
  par(mar=c(3, 4, 1.5, 0))

  print(paste(direction, "_", target_chr, sep=""))
  #Get the chromosome or region of interest
  #target_chr = "chr17"
  i = which(cnvs[,"CHR"] == target_chr)

  #Get the ideogram value for the current chromosome
  ideo_thischr = ideo_data[which(ideo_data[,"chrom"] == target_chr), ]

  #Define some display genes.  Those in the input gene list on the target chr
  #They should also be below the lowest cut or above the highest cut
  gi_up = which(genes[,"chr"] == target_chr & (genes[,"mean_cnv_diff"] >= cut5 | genes[,"cnvhmm_status"]=="Gain"))
  gi_down = which(genes[,"chr"] == target_chr & (genes[,"mean_cnv_diff"] <= cut2 | genes[,"cnvhmm_status"]=="Loss"))

  #Assign colors according to the magnitude of the DIFF
  cnvs[,"COLOR"] = "grey50"
  ii = which(cnvs[,"DIFF"] <= cut1 & cnvs[,"DIFF"] > cut2); if (length(ii) > 0){cnvs[ii, "COLOR"] = losscolor1}
  ii = which(cnvs[,"DIFF"] <= cut2 & cnvs[,"DIFF"] > cut3); if (length(ii) > 0){cnvs[ii, "COLOR"] = losscolor2}
  ii = which(cnvs[,"DIFF"] <= cut3); if (length(ii) > 0){cnvs[ii, "COLOR"] = losscolor3}
  ii = which(cnvs[,"DIFF"] >= cut4 & cnvs[,"DIFF"] < cut5); if (length(ii) > 0){cnvs[ii, "COLOR"] = gaincolor1}
  ii = which(cnvs[,"DIFF"] >= cut5 & cnvs[,"DIFF"] < cut6); if (length(ii) > 0){cnvs[ii, "COLOR"] = gaincolor2}
  ii = which(cnvs[,"DIFF"] >= cut6); if (length(ii) > 0){cnvs[ii, "COLOR"] = gaincolor3}

  #Grab the CNV and POS values for that region.  Divide into two groups to help with skewed scale
  total_windows=length(i)
  display_cut_up=-2
  display_cut_down=1
  cnvs_chr=cnvs[i,]
  gain_i = which(cnvs_chr[,"DIFF"] > display_cut_up)
  loss_i = which(cnvs_chr[,"DIFF"] < display_cut_down)
  xg=cnvs_chr[gain_i,"MID"]
  yg=cnvs_chr[gain_i,"DIFF"]
  cg=cnvs_chr[gain_i,"COLOR"]
  xl=cnvs_chr[loss_i,"MID"]
  yl=cnvs_chr[loss_i,"DIFF"]
  cl=cnvs_chr[loss_i,"COLOR"]

  #Grab the CNV HMM segments
  j = which(segments[,"CHR"] == target_chr)
  segments_chr=segments[j,]
  gain_j = which(segments_chr[,"Status"] == "Gain")
  loss_j = which(segments_chr[,"Status"] == "Loss")

  yu1 = max(yg)
  yu2 = max(yg)
  yu3 = max(yg)
  if (length(gi_up) > 0){
    #If there are going to be text gene labels, make room for them
    yu1 = max(yg) + (max(yg)+display_cut_down)*0.3
    yu2 = max(yg) + (max(yg)+display_cut_down)*0.7
    yu3 = max(yg) + (max(yg)+display_cut_down)*1.1
  }
  yd1 = -2
  yd2 = -2
  yd3 = -2
  if(length(gi_down) > 0){
    #If there are going to be text gene labels, make room for them
    yd1 = (abs(min(yl))+(2*0.4)) *(-1)
    yd2 = (abs(min(yl))+(2*1.0)) *(-1)
    yd3 = (abs(min(yl))+(2*1.4)) *(-1)
  }

  #Get the limits for the X axis
  x_ideo_chr = c(ideo_thischr[,"chromStart"], ideo_thischr[,"chromEnd"])
  xlim_lower = min(x_ideo_chr)
  xlim_upper = max(x_ideo_chr)

  #Axis labels
  xlabel=paste("Position (bp) on ", target_chr, sep="")
  ylabel="CNV Difference"

  #Plot 2 - Plot the values for the gains
  if (direction == "gains"){
    mainlabel="Gains"
    ylim_lower=display_cut_up; if (ylim_lower > -2){ylim_lower = -2}
    ylim_upper=yu3; if (ylim_upper < 2){ylim_upper = 2}
    plot(x=xg, y=yg, pch=16, col=cg, ylim=c(ylim_lower, ylim_upper), xlim = c(xlim_lower, xlim_upper), xlab="", ylab=ylabel, main=mainlabel)
    abline(h=0, lty=3, lwd=0.5, col="black")
    abline(h=cut4, lty=2, lwd=0.5, col=gaincolor1)
    abline(h=cut5, lty=2, lwd=0.5, col=gaincolor2)
    abline(h=cut6, lty=2, lwd=0.5, col=gaincolor3)
    #If a hardcap is being utilized... mark it
    if (max(yg) >= hard_cap_upper){
      abline(h=hard_cap_upper, lty=3, lwd=0.5, col=gaincolor3)
    }
    cex_text=0.75
    if (length(gi_up) > 10){cex_text=0.6}
    if(length(gi_up) > 0){ 
      odds_up=seq(1, length(gi_up), 2)
      x0_odd=genes[gi_up[odds_up],"Mid"]; x1_odd=x0_odd; y0_odd=rep(cut5, length(x0_odd)); y1_odd=rep(yu1, length(x0_odd))
      segments(x0=x0_odd, x1=x1_odd, y0=y0_odd, y1=y1_odd, col=gaincolor2, lty=2, lwd=0.5)
      text(x=genes[gi_up[odds_up],"Mid"], y=yu1, labels=genes[gi_up[odds_up],"mapped_gene_name"], srt=45, cex=cex_text, col=gaincolor4)
    }
    if(length(gi_up) > 1){ 
      evens_up=seq(2, length(gi_up), 2)
      x0_even=genes[gi_up[evens_up],"Mid"]; x1_even=x0_even; y0_even=rep(cut5, length(x0_even)); y1_even=rep(yu2, length(x0_even))
      segments(x0=x0_even, x1=x1_even, y0=y0_even, y1=y1_even, col=gaincolor2, lty=2, lwd=0.5)
      text(x=genes[gi_up[evens_up],"Mid"], y=yu2, labels=genes[gi_up[evens_up],"mapped_gene_name"], srt=45, cex=cex_text, col=gaincolor4)
    }
    if (length(gain_j)>0){
      segx1g = segments_chr[gain_j,"START"]
      segx2g = segments_chr[gain_j,"END"]
      segyg = segments_chr[gain_j,"Adjusted_CN_DIFF"]
      segments(x0=segx1g, x1=segx2g, y0=segyg, y1=segyg, col=gaincolor4, lty=1, lwd=2)
    }
  }

  #Plot 3 - Plot the values for the losses
  if (direction == "losses"){
    mainlabel="Losses"
    ylim_lower=yd3; if (ylim_lower > -2){ylim_lower = -2}
    ylim_upper=display_cut_down;
    plot(x=xl, y=yl, pch=16, col=cl, ylim=c(ylim_lower,ylim_upper), xlim = c(xlim_lower, xlim_upper), xlab=xlabel, ylab=ylabel, main=mainlabel)
    abline(h=0, lty=3, lwd=0.5, col="black")
    abline(h=cut1, lty=2, lwd=0.5, col=losscolor1)
    abline(h=cut2, lty=2, lwd=0.5, col=losscolor2)
    abline(h=cut3, lty=2, lwd=0.5, col=losscolor3)
    #If a hardcap is being utilized... mark it
    if (min(yl) <= hard_cap_lower){
      abline(h=hard_cap_lower, lty=3, lwd=0.5, col=losscolor3)
    }
    cex_text=0.75

    #Add text labels
    if (length(gi_down) > 10){cex_text=0.6}
    if (length(gi_down) > 0){
      odds_down=seq(1, length(gi_down), 2)
      x0_odd=genes[gi_down[odds_down],"Mid"]; x1_odd=x0_odd; y0_odd=rep(cut2, length(x0_odd)); y1_odd=rep(yd1, length(x0_odd))
      segments(x0=x0_odd, x1=x1_odd, y0=y0_odd, y1=y1_odd, col=losscolor1, lty=2, lwd=0.5)
      text(x=genes[gi_down[odds_down],"Mid"], y=yd1, labels=genes[gi_down[odds_down],"mapped_gene_name"], srt=45, cex=cex_text, col=losscolor4)
    }
    if (length(gi_down) > 1){
      evens_down=seq(2, length(gi_down), 2)
      x0_even=genes[gi_down[evens_down],"Mid"]; x1_even=x0_even; y0_even=rep(cut2, length(x0_even)); y1_even=rep(yd2, length(x0_even))
      segments(x0=x0_even, x1=x1_even, y0=y0_even, y1=y1_even, col=losscolor1, lty=2, lwd=0.5)
      text(x=genes[gi_down[evens_down],"Mid"], y=yd2, labels=genes[gi_down[evens_down],"mapped_gene_name"], srt=45, cex=cex_text, col=losscolor4)
    }
    if (length(loss_j)>0){
      segx1l = segments_chr[loss_j,"START"]
      segx2l = segments_chr[loss_j,"END"]
      segyl = segments_chr[loss_j,"Adjusted_CN_DIFF"]
      segments(x0=segx1l, x1=segx2l, y0=segyl, y1=segyl, col=losscolor4, lty=1, lwd=2)
    }
  }

  #Plot 4 - Plot the legend for the ideogram
  ideo_legend_cols = c("white","grey75","grey50","grey25","black","firebrick1","firebrick2","firebrick")
  ideo_legend_text = c("gneg","gpos25","gpos50","gpos75","gpos100","gvar","stalk","acen")
  plot.new()
  par(mar=c(0,0,0,0))
  legend("center", col=ideo_legend_cols, pch=15, legend=ideo_legend_text, title="Giemsa staining", cex=0.9)

  #Plot 5 - Plot the legend for the gains
  if (direction == "gains"){
    red_count = length(which(cnvs_chr[,"COLOR"] == gaincolor3))
    orange_count = length(which(cnvs_chr[,"COLOR"] == gaincolor2))
    yellow_count = length(which(cnvs_chr[,"COLOR"] == gaincolor1))

    plot.new()
    par(mar=c(0,0,0,0))
    gain_legend_text=c(paste("Gain > ", cut6, " (n = ", prettyNum(red_count, big.mark=",", scientific=FALSE), ")", sep=""), 
                       paste("Gain > ", cut5, " (n = ", prettyNum(orange_count, big.mark=",", scientific=FALSE), ")", sep=""), 
                       paste("Gain > ", cut4, " (n = ", prettyNum(yellow_count, big.mark=",", scientific=FALSE), ")", sep=""),
                       "HMM Gain Segs")
    gain_legend_cols=c(gaincolor3, gaincolor2, gaincolor1, gaincolor4)
    gain_legend_title=paste("Gain (Windows = ", prettyNum(total_windows, big.mark=",", scientific=FALSE), ")", sep="")
    legend("center", col=gain_legend_cols, pch=c(16,16,16,NA), lty=c(2,2,2,1), lwd=c(1,1,1,2), legend=gain_legend_text, title=gain_legend_title, cex=0.9)
  }

  #Plot 6 - Plot the legend for the losses
  if (direction == "losses"){
    lightblue_count = length(which(cnvs_chr[,"COLOR"] == losscolor1))
    blue_count = length(which(cnvs_chr[,"COLOR"] == losscolor2))
    darkblue_count = length(which(cnvs_chr[,"COLOR"] == losscolor3))
    plot.new()
    par(mar=c(0,0,0,0))
    loss_legend_text=c(paste("Loss < ", cut1, " (n = ", prettyNum(lightblue_count, big.mark=",", scientific=FALSE), ")", sep=""), 
                       paste("Loss < ", cut2, " (n = ", prettyNum(blue_count, big.mark=",", scientific=FALSE), ")", sep=""), 
                       paste("Loss < ", cut3, " (n = ", prettyNum(darkblue_count, big.mark=",", scientific=FALSE), ")", sep=""),
                       "HMM Loss Segs")
    loss_legend_cols=c(losscolor1, losscolor2, losscolor3, losscolor4)
    loss_legend_title=paste("Loss (Windows = ", prettyNum(total_windows, big.mark=",", scientific=FALSE), ")", sep="")
    legend("center", col=loss_legend_cols, pch=c(16,16,16,NA), lty=c(2,2,2,1), lwd=c(1,1,1,2), legend=loss_legend_text, title=loss_legend_title, cex=0.9)
  }
}

#############################################################################################################################################
#plotChrCNV_Compact - Function to plot CNV values by position in a compact way that allows all chromosomes to be visualized at once         #
#############################################################################################################################################

#Create a version of this function that produces a more compact display allowing many to be plotted on the same page
plotChrCNV_Compact = function(target_chr, type){
	
  #Get the chromosome or region of interest
  i = which(cnvs[,"CHR"] == target_chr)
  gi_up = which(genes[,"chr"] == target_chr & (genes[,"mean_cnv_diff"] >= cut5 | genes[,"cnvhmm_status"]=="Gain"))
  gi_down = which(genes[,"chr"] == target_chr & (genes[,"mean_cnv_diff"] <= cut2 | genes[,"cnvhmm_status"]=="Loss"))
  cnvs[,"COLOR"] = "grey50"
  ii = which(cnvs[,"DIFF"] <= cut1 & cnvs[,"DIFF"] > cut2); if (length(ii) > 0){cnvs[ii, "COLOR"] = losscolor1}
  ii = which(cnvs[,"DIFF"] <= cut2 & cnvs[,"DIFF"] > cut3); if (length(ii) > 0){cnvs[ii, "COLOR"] = losscolor2}
  ii = which(cnvs[,"DIFF"] <= cut3); if (length(ii) > 0){cnvs[ii, "COLOR"] = losscolor3}
  ii = which(cnvs[,"DIFF"] >= cut4 & cnvs[,"DIFF"] < cut5); if (length(ii) > 0){cnvs[ii, "COLOR"] = gaincolor1}
  ii = which(cnvs[,"DIFF"] >= cut5 & cnvs[,"DIFF"] < cut6); if (length(ii) > 0){cnvs[ii, "COLOR"] = gaincolor2}
  ii = which(cnvs[,"DIFF"] >= cut6); if (length(ii) > 0){cnvs[ii, "COLOR"] = gaincolor3}
  display_cut_up=-2
  display_cut_down=1
  cnvs_chr=cnvs[i,]
  gain_i = which(cnvs_chr[,"DIFF"] > display_cut_up)
  loss_i = which(cnvs_chr[,"DIFF"] < display_cut_down)
  xg=cnvs_chr[gain_i,"MID"]
  yg=cnvs_chr[gain_i,"DIFF"]
  cg=cnvs_chr[gain_i,"COLOR"]
  xl=cnvs_chr[loss_i,"MID"]
  yl=cnvs_chr[loss_i,"DIFF"]
  cl=cnvs_chr[loss_i,"COLOR"]

  #Get the ideogram value for the current chromosome
  ideo_thischr = ideo_data[which(ideo_data[,"chrom"] == target_chr), ]

  #Get the limits for the X axis
  x_ideo_chr = c(ideo_thischr[,"chromStart"], ideo_thischr[,"chromEnd"])
  xlim_lower = min(x_ideo_chr)
  xlim_upper = max(x_ideo_chr)

  #Grab the CNV HMM segments
  j = which(segments[,"CHR"] == target_chr)
  segments_chr=segments[j,]
  gain_j = which(segments_chr[,"Status"] == "Gain")
  loss_j = which(segments_chr[,"Status"] == "Loss")

  if (type == "GAIN"){
    ylim_lower = min(yg); if (ylim_lower > -2){ylim_lower = -2}
    ylim_upper = max(yg); if (ylim_upper < 2){ylim_upper = 2}
    plot(x=xg, y=yg, pch=16, col=cg, xlab="", ylab="CNV", main=target_chr, ylim=c(ylim_lower, ylim_upper), xlim=c(xlim_lower, xlim_upper))
  }else if (type == "LOSS"){
    ylim_lower = min(yl); if (ylim_lower > -2){ylim_lower = -2}
    ylim_upper = 1
    plot(x=xl, y=yl, pch=16, col=cl, xlab="", ylab="CNV", main=target_chr, ylim=c(ylim_lower, ylim_upper), xlim=c(xlim_lower, xlim_upper))
  }else if (type == "BOTH"){
    ylim_lower = -2
    ylim_upper = 6
    yl[which(yl > 5.5)] = 5.5
    yl[which(yl < -2)] = -2
    plot(x=xl, y=yl, pch=16, col=cl, xlab="", ylab="CNV", main=target_chr, ylim=c(ylim_lower, ylim_upper), xlim=c(xlim_lower, xlim_upper))
  }
  if (length(gain_j)>0){
    segx1g = segments_chr[gain_j,"START"]
    segx2g = segments_chr[gain_j,"END"]
    segyg = segments_chr[gain_j,"Adjusted_CN_DIFF"]
    segments(x0=segx1g, x1=segx2g, y0=segyg, y1=segyg, col=gaincolor4, lty=1, lwd=2)
  }
  if (length(loss_j)>0){
    segx1l = segments_chr[loss_j,"START"]
    segx2l = segments_chr[loss_j,"END"]
    segyl = segments_chr[loss_j,"Adjusted_CN_DIFF"]
    segments(x0=segx1l, x1=segx2l, y0=segyl, y1=segyl, col=losscolor4, lty=1, lwd=2)
  }
  abline(h=0, lty=2, lwd=1, col="black")
  abline(h=cut4, lty=2, lwd=0.5, col=gaincolor1)
  abline(h=cut5, lty=2, lwd=0.5, col=gaincolor2)
  abline(h=cut6, lty=2, lwd=0.5, col=gaincolor3)
  abline(h=cut1, lty=2, lwd=0.5, col=losscolor1)
  abline(h=cut2, lty=2, lwd=0.5, col=losscolor2)
  abline(h=cut3, lty=2, lwd=0.5, col=losscolor3)
}

#############################################################################################################################################
#drawIdeogram - Draw a chromsome ideogram for the coordinates being considered - A chromosome or piece of a chromosome                      #
#############################################################################################################################################
drawIdeogram = function(ideo_data, target_chr){

  main_title = paste("Copy number differences for ", target_chr, " (", category_name, " labeled)", sep = "")

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


#############################################################################################################################################
#openImageFile - Deal with open file handles for different image types                                                                      #
#############################################################################################################################################
#Create a function to deal with opening each possible image file desired
openImageFile = function(name, type, image_width, image_height){
  #Adjust file name if a chromosome range was specified
  if (chr_start > 0){
    name = paste(name, "_", format(chr_start, scientific=FALSE), "-", format(chr_end, scientific=FALSE), sep="")
  }
  if(type == "pdf"){
    filename = paste(name, ".pdf", sep="")
    pdf(filename, width=image_width, height=image_height)
  }else if(type == "jpeg"){
    filename = paste(name, ".jpeg", sep="")
    jpeg(filename=filename, width=image_width, height=image_height, units="in", pointsize=12, quality=100, bg = "white", res=150)
  }else if (type == "png"){
    filename = paste(name, ".png", sep="")
    png(filename=filename, width=image_width, height=image_height, units="in", pointsize=12, bg = "white", res=150)
  }else if (type == "bmp"){
    filename = paste(name, ".bmp", sep="")
    bmp(filename=filename, width=image_width, height=image_height, units="in", pointsize=12, bg = "white", res=150)
  }else if (type == "tiff"){
    filename = paste(name, ".tiff", sep="")
    tiff(filename=filename, width=image_width, height=image_height, units="in", pointsize=12, bg = "white", res=300, compression="none")
  }
}


#SCRIPT SECTION

#Load data
cnvs=read.table(cnv_file, comment.char="#", header=TRUE)
segments=read.table(segments_file, sep="\t", as.is=c(1,11), header=TRUE)
genes=read.table(gene_file, sep="\t", header=TRUE, as.is=c(1:5,8,10))
ideo_data = read.table(ideogram_file, sep="\t", header=FALSE, comment.char="#", as.is=c(1,4,5))
names(ideo_data) = c("chrom","chromStart", "chromEnd", "name", "gieStain")

#Define the image file type to be produced (pdf, jpeg, png, bmp, tiff, none)
image_width_1 = 12
image_height_1 = 5
image_width_2 = 14
image_height_2 = 9

#Define some hard caps for upper and lower display limits
hard_cap_upper=15  #Tried 20 for a while. Trying 15 now
hard_cap_lower=-2

#Change output dir
setwd(outdir)

#List of chrs to process if user specified 'ALL'
chr_list = chr
if (chr == "ALL"){
  chr_list=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")
}

#Add "chr" to the chromosome names in the cnv and genes objects
cnvs[,"CHR"]=paste("chr", cnvs[,"CHR"], sep="")
genes[,"chr"]=paste("chr", genes[,"chr"], sep="")

#Make sure segments file has at least one row of data, otherwise following commands will choke, need to be skipped  
if (length(rownames(segments))>0){
  #Add "chr" to the chromosome names in the segments object
  segments[,"CHR"]=paste("chr", segments[,"CHR"], sep="")
  #Calculate the adjusted CN difference for segments
  segments[,"Adjusted_CN_DIFF"]=segments[,"Adjusted_CN1"]-segments[,"Adjusted_CN2"]
}else{ #Just add colnames to empty matrix to prevent further downstream errors
  segments=as.data.frame(setNames(replicate(length(colnames(segments))+1,numeric(0), simplify = F), c(colnames(segments),"Adjusted_CN_DIFF")))
}

#If the user specified a subrange of a single chromosome, remove all other data
if (length(chr_list) == 1 & chr_start > 0){
  i = which ((cnvs[,"CHR"] == chr) & (cnvs[,"POS"] >= chr_start) & (cnvs[,"POS"] <= chr_end))
  cnvs=cnvs[i,]
}

#Reorder genes by start position
o=order(genes[,"start"])
genes=genes[o,]

#Calculate the midpoint for CNV window coords
window_size=cnvs[2,"POS"]-cnvs[1,"POS"]
cnvs[,"MID"] = cnvs[,"POS"]+(window_size/2)

#Calculate the midpoint for gene coords
genes[,"Mid"] = genes[,"start"]+((genes[,"end"]-genes[,"start"])/2)

#Reset values smaller than -2 to be -2 (both copies deleted)
if (length(which(cnvs[,"DIFF"] < hard_cap_lower)) > 0){
  cnvs[which(cnvs[,"DIFF"] < hard_cap_lower),"DIFF"] = hard_cap_lower
}

#Reset values larger than 20 to be 20 (arbitrary - for display purposes).  Single outliers, usually false positives near the centromeres obscure the data by jacking up the scale...
if (length(which(cnvs[,"DIFF"] > hard_cap_upper)) > 0){
  cnvs[which(cnvs[,"DIFF"] > hard_cap_upper),"DIFF"] = hard_cap_upper
  segments[which(segments[,"Adjusted_CN_DIFF"] > hard_cap_upper),"Adjusted_CN_DIFF"] = hard_cap_upper
}

#TODO: Apply a further more aggressive hard-capping near centromeres where outliers in the gain direction are common?
#for (chr_name in chr_list){
#
#
#}


#Define some display cutoffs
cut1=-0.5
cut2=-0.75
cut3=-0.9
cut4=1
cut5=2
cut6=5

#Define gain/loss colors
gaincolor1="#FDCC8A"
gaincolor2="#FC8D59"
gaincolor3="#E34A33"
gaincolor4="#B30000"
losscolor1="#BDC9E1"
losscolor2="#74A9CF"
losscolor3="#2B8CBE"
losscolor4="#045A8D"

#Generate a figure for each chromosome
for (chr_name in chr_list){
  if (image_type == "none"){
    plotChrCNV(chr_name, "gains")
    plotChrCNV(chr_name, "losses")
  }else{
    #Gains
    file_name = paste("Gains_", chr_name, sep="")
    openImageFile(file_name, image_type, image_width_1, image_height_1)
    plotChrCNV(chr_name, "gains")
    dev.off()

    #Losses
    file_name = paste("Losses_", chr_name, sep="")
    openImageFile(file_name, image_type, image_width_1, image_height_1)
    plotChrCNV(chr_name, "losses")
    dev.off()
  }
}

#Adjust margins -> c(bottom, left, top, right)
margins=c(2, 3.75, 1.75, 1.75)+0.1

#Set up the legend
legend_text=c(paste("Copy number gain > ", cut6, sep=""),
              paste("Copy number gain > ", cut5, sep=""),
              paste("Copy number gain > ", cut4, sep=""),
              paste("Copy number loss < ", cut1, sep=""),
              paste("Copy number loss < ", cut2, sep=""),
              paste("Copy number loss < ", cut3, sep=""),
              "HMM Gain Seg",
              "HMM Loss Seg")
legend_cex=1

#Generate a figure for all chromosome gains on a single page - only when processing all chromosomes though...
if (length(chr_list) > 1){
  if (image_type != "none"){
    openImageFile("Gains_AllChrs", image_type, image_width_2, image_height_2)
  }
  par(mfrow=c(6,4), mar=margins)
  print("All chromosomes - Gains")
  for (chr_name in chr_list){
    plotChrCNV_Compact(chr_name, "GAIN")
  }
  plot.new()
  par(mar=c(0,0,0,0))
  legend("center", col=c(gaincolor3, gaincolor2, gaincolor1,losscolor1, losscolor2, losscolor3, gaincolor4, losscolor4), pch=c(16,16,16,16,16,16,NA,NA), lty=c(2,2,2,2,2,2,1,1), lwd=c(1,1,1,1,1,1,2,2), legend=legend_text, title="Cutoffs", cex=legend_cex)
  if (image_type != "none"){
    dev.off()
  }
}

#Generate a figure for all chromosome losses on a single page
if (length(chr_list) > 1){
  if (image_type != "none"){
    openImageFile("Losses_AllChrs", image_type, image_width_2, image_height_2)
  }
  par(mfrow=c(6,4), mar=margins)
  print("All chromosomes - Losses")
  for (chr_name in chr_list){
    plotChrCNV_Compact(chr_name, "LOSS")
  }
  plot.new()
  par(mar=c(0,0,0,0))
  legend("center", col=c(gaincolor3, gaincolor2, gaincolor1,losscolor1, losscolor2, losscolor3, gaincolor4, losscolor4), pch=c(16,16,16,16,16,16,NA,NA), lty=c(2,2,2,2,2,2,1,1), lwd=c(1,1,1,1,1,1,2,2), legend=legend_text, title="Cutoffs", cex=legend_cex)
  if (image_type != "none"){
    dev.off()
  }
}

if (length(chr_list) > 1){
  if (image_type != "none"){
    openImageFile("Both_AllChrs", image_type, image_width_2, image_height_2)
  }
  par(mfrow=c(6,4), mar=margins)
  print("All chromosomes - Gains and Losses")
  for (chr_name in chr_list){
    plotChrCNV_Compact(chr_name, "BOTH")
  }
  plot.new()
  par(mar=c(0,0,0,0))
  legend("center", col=c(gaincolor3, gaincolor2, gaincolor1,losscolor1, losscolor2, losscolor3, gaincolor4, losscolor4), pch=c(16,16,16,16,16,16,NA,NA), lty=c(2,2,2,2,2,2,1,1), lwd=c(1,1,1,1,1,1,2,2), legend=legend_text, title="Cutoffs", cex=legend_cex)
  if (image_type != "none"){
    dev.off()
  }
}




