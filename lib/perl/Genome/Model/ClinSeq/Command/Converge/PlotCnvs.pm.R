#From Chris Miller, Charles Lu and Zach Skidmore, assumes sample is Male by default
#if(!gtools %in% installed.packages()) install.packages(gtools)
#package needs to be installed, this might do it automatically but if not, install the package gtools.
require(gtools)

usage <- function() {
  print("Rscript PlotCnvs.pm.R microarray_cnv wgs_cnv exome_cnv op_plot_file sample_name");
  print("The *cnv_files are of the form 'chr\tstart\tend\tn_markers(NA if unavail)\ttumor_CN', without headers");
}

plotEm <- function(files, sample_name){
  chrNames=TRUE
  for(file in files){    
    infile = NULL;
    id = "";
    if(!(file=="header")){
      infile=file
      print(file)
      cnvfile_suffix = ".cnvs.txt.cn"
      name_prefix = paste(sample_name, ".", sep = "")
      id = gsub(cnvfile_suffix, "", basename(file))
      id = gsub(name_prefix, "", id)
      id = gsub("wgs", "WGS", id)
      id = gsub("exome", "Exome", id)
      id = gsub("microarrays", "Microarray", id)
    }
    plotSegments(filename=infile, 
                 showNorm=TRUE, multiplePlot=TRUE, drawLabel=chrNames,
                 title=chrNames, sampleName=id)
    chrNames=FALSE;
  }
}


addLegend <- function(){
  plot(0, 0, xlim=c(0,100), ylim=c(0,100),xaxt="n",yaxt="n", bty="n",
      bg = "white", type = "n")
  rect(seq(40,60,0.2),50,seq(40.2,60.2,0.2),70,col=getColor(seq(0,4,0.04)),
      border = 0, xaxt="n")
  text(40,50,"0",pos=1)
  text(45,50,"1",pos=1)
  text(50,50,"2",pos=1)
  text(55,50,"3",pos=1)
  text(60,50,"4+",pos=1)
  text(50,100,sample_name, pos=1, cex = 2)
}


##---------------------------------------------------------
## create the Female entrypoints table
##
createFemaleEntrypoints <- function(){
  p = data.frame(0, 250000000, 2)
  colnames(p) = c("chr", "length", "ploidy")
  p = rbind(c(1, 249250621, 2), p)
  p = rbind(c(2, 243199373, 2), p)
  p = rbind(c(3, 198022430, 2), p)
  p = rbind(c(4, 191154276, 2), p)
  p = rbind(c(5, 180915260, 2), p)
  p = rbind(c(6, 171115067, 2), p)
  p = rbind(c(7, 159138663, 2), p)
  p = rbind(c(8, 146364022, 2), p)
  p = rbind(c(9, 141213431, 2), p)
  p = rbind(c(10, 135534747, 2), p)
  p = rbind(c(11, 135006516, 2), p)
  p = rbind(c(12, 133851895, 2), p)
  p = rbind(c(13, 115169878, 2), p)
  p = rbind(c(14, 107349540, 2), p)
  p = rbind(c(15, 102531392, 2), p)
  p = rbind(c(16, 90354753, 2), p)
  p = rbind(c(17, 81195210, 2), p)
  p = rbind(c(18, 78077248, 2), p)
  p = rbind(c(20, 63025520, 2), p)
  p = rbind(c(19, 59128983, 2), p)
  p = rbind(c(21, 48129895, 2), p)
  p = rbind(c(22, 51304566, 2), p)
  p = rbind(c("X", 155270560, 2), p)
  colnames(p) = c("chr", "length", "ploidy")
  p$ploidy<-as.numeric(p$ploidy)
  p$length<-as.numeric(p$length)
  index <- mixedsort(p$chr)
  p <- p[match(index, p$chr),]
  return(p)
}

##---------------------------------------------------------
## create the Male entrypoints table
##
createMaleEntrypoints <- function(){
  p = data.frame(0, 250000000, 2)
  p = rbind(c(1, 249250621, 2), p)
  p = rbind(c(2, 243199373, 2), p)
  p = rbind(c(3, 198022430, 2), p)
  p = rbind(c(4, 191154276, 2), p)
  p = rbind(c(5, 180915260, 2), p)
  p = rbind(c(6, 171115067, 2), p)
  p = rbind(c(7, 159138663, 2), p)
  p = rbind(c(8, 146364022, 2), p)
  p = rbind(c(9, 141213431, 2), p)
  p = rbind(c(10, 135534747, 2), p)
  p = rbind(c(11, 135006516, 2), p)
  p = rbind(c(12, 133851895, 2), p)
  p = rbind(c(13, 115169878, 2), p)
  p = rbind(c(14, 107349540, 2), p)
  p = rbind(c(15, 102531392, 2), p)
  p = rbind(c(16, 90354753, 2), p)
  p = rbind(c(17, 81195210, 2), p)
  p = rbind(c(18, 78077248, 2), p)
  p = rbind(c(19, 59128983, 2), p)
  p = rbind(c(20, 63025520, 2), p)
  p = rbind(c(21, 48129895, 2), p)
  p = rbind(c(22, 51304566, 2), p)
  p = rbind(c("X", 155270560, 1), p)
  p = rbind(c("Y", 59373566, 1), p)
  colnames(p) = c("chr", "length", "ploidy")
  p$ploidy<-as.numeric(p$ploidy)
  p$length<-as.numeric(p$length)
  index <- mixedsort(p$chr)
  p <- p[match(index, p$chr),]
  return(p)
}

##---------------------------------------------------------
## read in the entrypoints file
##
readEntrypoints <- function(file){
  p=read.table(file,sep="\t",quote="",comment.char="#",
    colClasses=c("character","numeric","numeric"))
  #p = cbind(data.frame(0,300000000,2),p)
   index <- mixedsort(p$V1)
   p <- p[match(index, p$V1),]
   return(p)
}

##---------------------------------------------------------
## add offsets to entrypoints file
##
addOffsets <- function(df){
  ##starts with chr,len,ploidy
  offsets = c()
  sum = 0
  for(i in 1:length(df[,1])){
    offsets = c(offsets,sum)
    sum = sum + df[i,2]
  }
  df[,4] <- offsets
  return(df)
}

#color gradient
getColor <- function(val){ 
  val[val>4] = 4
  val[val<0] = 0
  len = length(val)
  ## add in the extremes so that our palette is centered properly
  val = c(val,0,4)
  colFunc <- colorRampPalette(c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6"),
      bias=1)
  return(colFunc(100)[as.numeric(cut(val,breaks=100))[1:len]])
}


## main function - plot the segments
plotSegments <- function(chr="ALL", filename, ymax=0, ymin=1,
                         highlights=NULL, lowRes=FALSE, lowResMin=NULL,
                         lowResMax=NULL, showNorm=FALSE, baseline=2,
                         gainThresh=2.5, lossThresh=1.5, annotationsTop=NULL,
                         annotationsBottom = NULL, plotTitle="",
                         gainColor="red", lossColor="blue", ylabel="",
                         xmin=NULL, xmax=NULL, label_size=0.6, multiplePlot=FALSE, 
                        drawLabel=TRUE, percentagePlot=FALSE, title=FALSE, sampleName=""){
  ## add options for plotting just a smaller region - TODO
  xlim = NULL
  ## read in the segments, handle case of empty file
  segs = data.frame(V1=numeric(0),V2=numeric(0),V3=numeric(0),V4=numeric(0),V5=numeric(0))  
  print(filename)
  if(!is.null(filename)){
    if(file.info(filename)$size > 0){
      segs=read.table(filename,comment.char="#")
    }
  }
  ## read in the entrypoints
  entrypoints = addOffsets(createMaleEntrypoints())
  entrypoints = addOffsets(entrypoints)
  names(entrypoints) = c("chr","length","ploidy","offset")
  ## if we have regions to highlight, read them in too
  hlRegions = NULL;
  if(!(is.null(highlights))){
    hlRegions=read.table(highlights)
  }
  ## if we have annotations, read them in too
  annTopRegions=NULL
  if(!(is.null(annotationsTop))){
    annTopRegions=read.table(annotationsTop,comment.char="#")
  }
  annBtmRegions=NULL
  if(!(is.null(annotationsBottom))){
    annBtmRegions=read.table(annotationsBottom,comment.char="#")
  }

  ##validate that we have entrypoints for all of our chromosomes
  chrnames = names(table(segs$V1))
  if(!is.null(chrnames)){
    for(i in 1:length(chrnames)){
  ##raise an error if entrypoints and chrs don't match
      if(length(which(entrypoints$chr==chrnames[i])) < 1){
        cat("\nERROR - no entrypoint found for chromosome ",chrnames[i]," found in segs file\n")
          cat("maybe you meant to use the male entrypoints?\n")
          q(save="no")
      }
    }
  }

  ##---------------------------------------------------------
  ## function to expand the size of features
  ## so that they exceed the minimum pixel size on small
  ## plots
  makeVisible <- function(segs,minSize=lr.min){
    pos = which((segs[,3]-segs[,2]) < minSize)
    mid = ((segs[pos,3]-segs[pos,2])/2)+segs[pos,2]
    segs[pos,2] = mid - (minSize/2)
    segs[pos,3] = mid + (minSize/2)
    return(segs)
  }

##---------------------------------------------------------
## add annotations to the top/btm of the plot
## with lines to the peaks
addAnnos <- function(annos, segs, top=TRUE, offset=FALSE, chr=NULL, leftEdge=50000000){
    if(!(is.null(chr))){
      annos=annos[which(annos[,1] == chr),]
    }
    if(length(annos[,1]) < 1){
      return(0)
    }
    
    ypos = ymin*0.8
    if (top){
      ypos = ymax*0.8
    }
    
    for(i in 1:length(annos[,1])){
      st=as.numeric(annos[i,2])
      sp=as.numeric(annos[i,3])
      mid=(sp-st)+st
      midNoOffset = mid

      if(offset){
        offsetNum=as.numeric(entrypoints[which(entrypoints$chr==annos[i,1]),4])
        st = st + offsetNum
        sp = sp + offsetNum
        mid = mid + offsetNum
      }
      
      ##get the height of the peak at this position (if it exists)
      ptop = baseline
      peakNum=which((segs[,1] == as.character(annos[i,1])) &
        (segs[,2] <= midNoOffset) & (segs[,3] >= midNoOffset))

      if(length(peakNum > 0)){
        if(top){
          ptop = max(segs[peakNum,5])+((ymax-as.numeric(baseline))*0.05)
        } else {
          ptop = min(segs[peakNum,5])+((ymin-as.numeric(baseline))*0.05)
        }
      }

      ##adjust at left edge of plot
      ## todo - make this a percentage of the plot, rather than a
      ## set distance
      mid2 = mid
      if(mid < leftEdge){
        mid2 = leftEdge
      }
      
      if(is.null(annos[i,5])){
        annos[i,5] = 0
      } else {
        if (is.na(annos[i,5])){
          annos[i,5] = 0
        }
      }
      text(mid2,(ypos+annos[i,5]),annos[i,4],cex=0.5,font=3)
      lines(c(mid,mid2),c(ptop,(ypos+annos[i,5])*.90))
    }
}

  
  ## if we haven't set a ymax/ymin, set it to be just
  ## higher than the maximum peaks
  if(is.null(ymax)){
    ymax=max(segs[,5])*1.1
  }
  if(is.null(ymin)){
    ymin=min(segs[,5])*1.1
  }

  ## set the xlim to the width of the genome
  if(is.null(xlim)){
    xlim=c(1,sum(entrypoints$length))
  }else{
    a = segs[which(((segs$V3 >= xlim[1]) & (segs$V3 <= xlim[2])) | ((segs$V2 >= xlim[1]) & (segs$V2 <= xlim[2]))),]
  }

  ## outline the plot
  if (percentagePlot == TRUE){
    plot(0, 0, xlim=xlim, ylim=c(ymin,ymax), pch=".", ylab=ylabel, xlab="",
         xaxt="n", yaxt="n", cex.lab=1, las=2, cex.axis=label_size)
    axis(side = 2, las=2, at = c(1,0.75, 0.5,0.25,0,-0.25,-0.5,-0.75,-1),
         cex.axis=label_size, labels=c("100%",  "75%", "50%", "25%", 0, "25%", "50%", "75%", "100%"))
  } else if (multiplePlot == TRUE){
      plot(0, 0, xlim=xlim, ylim=c(ymin, ymax), pch=".", ylab=ylabel, xlab="",
           xaxt = "n", yaxt="n", cex.lab=1, las=2, cex.axis=label_size*0.8, xaxs="i", yaxs="i")
  } else {
    plot(0, 0, xlim=xlim, ylim=c(ymin,ymax), pch=".", ylab=ylabel, xlab="",
         xaxt="n", cex.lab=1, cex.axis=label_size)
  }

  if(title==TRUE){
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey50")
  }
  title(ylab=ylabel,line=2,cex.lab=label_size)

  ## add the title
  if(!(is.null(plotTitle))){
    if (label_size < 0.8) {
      mtext(plotTitle,padj=-.5,cex=0.8,side=3);
    }
    else {
      mtext(plotTitle,padj=-.5,cex=label_size,side=3);
    }
  }

  ## draw baselines
  abline(h=baseline, col="grey50")

  offsets = as.numeric(entrypoints[,4])
  offsets = append(offsets,sum(entrypoints$length))

  ## draw chromosome labels if needed
  if (drawLabel){
    for(i in 2:(length(offsets)-1)){
      text((offsets[i]+offsets[i+1])/2, 0.5, labels= gsub("chr","",entrypoints[i,1]), cex=1,col="white")
    }    
  }
  
  ## draw highlight regions, if specified
  if(!(is.null(highlights))){      
    chrNames = names(table(hlRegions[,1]))
    for(i in 1:length(chrNames)){
      ##don't plot things that we aren't considering
      if(chrNames[i] %in% entrypoints$chr){   
        ## offset is equal to whatever chromosome we're on
        offset = as.numeric(entrypoints[which(entrypoints$chr==chrNames[i]),4])
        reg = which(hlRegions[,1] == chrNames[i])
        if(length(reg) > 0){          
          hlRegions[reg,2] = hlRegions[reg,2] + offset
          hlRegions[reg,3] = hlRegions[reg,3] + offset
        }
      }
      ## hlRegions[which(hlRegions[,1] == chrNames[i]),3] = hlRegions[which(hlRegions[,1] == chrNames[i]),3] + offset
    }
    
    if(lowRes){
      lrpos = which((hlRegions[,3]-hlRegions[,2]) > lowResMin)
      hlRegions[lrpos,] = makeVisible(hlRegions[lrpos,],lowResMax)
    }
    rect(hlRegions[,2], ymin*2, hlRegions[,3], ymax*2, col="gold",lwd=0,lty="blank")
  }



## function to actually draw the segments
drawSegs <- function(segs,color="black"){
  chrNames = names(table(segs[,1]))
  for(i in 1:length(chrNames)){
  ## offset is equal to whatever chromosome we're on
    offset = as.numeric(entrypoints[which(entrypoints$chr==chrNames[i]),4])
    toAdj = which(segs[,1] == chrNames[i])
    if(length(toAdj) > 0){
      segs[toAdj,2] = segs[toAdj,2] + offset
      segs[toAdj,3] = segs[toAdj,3] + offset
    }
  }

  ## do the lowres expansion if specified
  if(lowRes){
    lrpos = which((segs[,3]-segs[,2]) > lowResMin)
    segs[lrpos,] = makeVisible(segs[lrpos,],lowResMax)
  }

  ## draw the segments
  rect(segs[,2], 0, segs[,3], 1, col=getColor(segs[,5]),lty="blank")
}


  ## finally, do the drawing:

  ##plot normal
  if(showNorm){
    a2=segs[which((segs[,5] <= gainThresh) & (segs[,5] >=lossThresh)),]
    if(length(a2[,1])>0){
      drawSegs(a2,color="grey")
    }
  }
  ##plot gain
  a2=segs[which(segs[,5] > gainThresh),]
  if(length(a2[,1])>0){
    drawSegs(a2,color=gainColor)
  }
  ##plot loss
  a2=segs[which(segs[,5] < lossThresh),]
  if(length(a2[,1])>0){
    drawSegs(a2,color=lossColor)
  }

  ## draw chromosome labels
  abline(v = 0, col="gray25")
  for(i in 1:(length(offsets)-1)){
    abline(v = offsets[i+1], col="gray25")
  }

  #add sample titles
  if(title==FALSE){
    rect(0, 0, offsets[2], 1,col="grey50")
    text(0, 0.5, sampleName,pos=4,col="white",cex=1.2)
  } else {
    text(0, 0.5, "Chromosome",pos=4,col="white",cex=1.)
  }
  
  

  ## add top annotations, if specfied
  if(!(is.null(annotationsTop))){
    addAnnos(annTopRegions, segs, top=TRUE, offset=TRUE)
  }
  ## add btm annotations, if specfied
  if(!(is.null(annotationsBottom))){
    addAnnos(annBtmRegions, segs, top=FALSE, offset=TRUE)
  }
}

require(gtools)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 5) {
  usage()
  quit()
}

#MAIN()
microarray_cnv = args[1]
wgs_cnv = args[2]
exome_cnv = args[3]
plot_file = args[4]
sample_name = args[5]

#plotting params
pdf(plot_file, height=10, width=12);

##plot all samples, changed to sort by chromosome and by filename
files = c(microarray_cnv, wgs_cnv, exome_cnv) 
files = mixedsort(files)
files = c("header", files)
mat <- matrix(c(1, 2, 3, 4, 5), 5)
layout(mat, widths = c(1, 1, 1, 1), heights = c(0.5, 0.5, 0.5, 0.5))
oma = c(4,2,3,2)
par(oma=oma, mar=c(0, 1, 0, 1))
plotEm(files, sample_name)
addLegend()

#adds title, labels to X and Y axis
dev.off()
