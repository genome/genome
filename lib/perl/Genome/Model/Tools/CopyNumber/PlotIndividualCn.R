library(sqldf)

##---------------------------------------------------------
## read in the entrypoints file
##
readEntrypoints <- function(file){
  p=read.table(file,sep="\t",quote="",comment.char="#",
    colClasses=c("character","numeric","numeric"))
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



##---------------------------------------------------------
## main function - plot the segments
##

plotSegments <- function(chr="ALL", filename, entrypoints, ymax=NULL, ymin=NULL,
                         highlightedSegs=NULL, lowRes=FALSE, lowResMin=NULL,
                         lowResMax=NULL, showNorm=FALSE, baseline=2,
                         gainThresh=2.5, lossThresh=1.5, annotationsTop=NULL,
                         annotationsBottom = NULL, plotTitle="", gaps=NULL, tumorNormalRatio=1,
                         gainColor="indianred4", lossColor="midnightblue", correctedWindows=NULL,
                         xmin=NULL, xmax=NULL, label_size=0.6, pdfOpen=FALSE,
                         pdfClose=FALSE, pdfFile=NULL, pdfWidth=12, pdfHeight=9,
                         tumorWindows=NULL, normalWindows=NULL,coverageTracks=FALSE){

  if(pdfOpen){
    if(!is.null(pdfFile)){
      pdf(file=pdfFile, width=pdfWidth, height=pdfHeight)
    }
  }


  if(coverageTracks==TRUE){
    layout(matrix(c(1,2,3,4),4,1,byrow=TRUE), widths=c(1), heights=c(2,20,5,5))
    par(oma=c(1,1,1,1), mar=c(1.5,4,1,3))
    plot(-110,-110, xaxt="n", yaxt="n", ylab="", xlab="",xlim=c(0,100), ylim=c(0,100), bty="n")
    text(50,50,plotTitle, cex=1.5)
  }




  xlim = NULL

  ## read in the segments
  segs=read.table(filename,comment.char="#")

  ## read in the entrypoints
  entrypoints=addOffsets(readEntrypoints(entrypoints))
  names(entrypoints) = c("chr","length","ploidy","offset")

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
  for(i in 1:length(chrnames)){
                                        #raise an error if entrypoints and chrs don't match
    if(length(which(entrypoints$chr==chrnames[i])) < 1){
      cat("\nERROR - no entrypoint found for chromosome ",chrnames[i]," found in segs file\n")
      cat("maybe you meant to use the male entrypoints?\n")
      q(save="no")
    }
  }

  ##get this chromosome's entrypoints and segments
  entry=entrypoints[which(entrypoints$chr==chr),]
  segs = segs[which(segs$V1==chr),]

  ## if we haven't set a ymax, set it to be just
  ## higher than the maximum peak
  if(is.null(ymax)){
    ymax=max(segs[,5])*1.1
  }
  if(is.null(ymin)){
    ymin=min(segs[,5])*1.1
  }

  ## if there wasn't an xlim value passed in, use the whole chromosome
  ## otherwise, find the sub-region
  if(is.null(xmax)){
    xlim=c(1,entry$length)
  }else{
    if(!(is.null(xmin))){
      xlim=c(xmin,xmax)
    } else {
      xlim=c(1,xmax)
    }
    segs = segs[which(((segs$V3 >= xlim[1]) & (segs$V3 <= xlim[2])) | ((segs$V2 >= xlim[1]) & (segs$V2 <= xlim[2]))),]
  }

  

  ##draw the plot region
  plot(0,0,xlim=xlim,ylim=c(ymin,ymax),pch=".",ylab="Copy Number", xlab="position (Mb)",cex.lab=1, cex.axis=1)

  #color gap regions
  if(!is.null(gaps)){
    drawGaps(gaps, chr)
  }

  ##add segment highlights
  if(!is.null(highlightedSegs)){
    drawHighlights(highlightedSegs, chr)
  }
  
  ##add the title
  if(!(coverageTracks)){
    if(!(is.null(plotTitle))){
      title(main=plotTitle)
    }
  }

  ## draw baselines
  abline(h=baseline,col="grey25")

  ##grid lines
  xaxp=par()$xaxp
  abline(v=seq(xaxp[1], xaxp[2], (xaxp[2]-xaxp[1])/xaxp[3]), col="grey50", lty=3)
  yaxp=par()$yaxp
  abline(h=seq(yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3]), col="grey50", lty=3)

  ## --- do the plotting ---

  ##plot bins
  drawCorrectedWindows(correctedWindows)


  ##plot normal segs
  if(showNorm){
    a2=segs[which((segs[,5] <= gainThresh) & (segs[,5] >=lossThresh)),]
    if(length(a2[,1])>0){
      drawSegs(a2[,2], a2[,3], a2[,5], color="grey", type="l", lowRes, lowResMin, lowResMax, baseline)
    }
  }
  ##plot gain segs
  a2=segs[which(segs[,5] > gainThresh),]
  if(length(a2[,1])>0){
    drawSegs(a2[,2], a2[,3], a2[,5], color=gainColor, type="l", lowRes, lowResMin, lowResMax, baseline)
  }
  ##plot loss segs
  a2=segs[which(segs[,5] < lossThresh),]
  if(length(a2[,1])>0){
    drawSegs(a2[,2], a2[,3], a2[,5], color=lossColor, type="l", lowRes, lowResMin, lowResMax, baseline)
  }

  ## add top annotations, if specfied
  if(!(is.null(annotationsTop))){
    addAnnos(annTopRegions, segs, top=TRUE, chr=chr, leftEdge=5000000)
  }
  ## add btm annotations, if specfied
  if(!(is.null(annotationsBottom))){
    addAnnos(annBtmRegions, segs, top=FALSE, chr=chr, leftEdge=5000000)
  }


  

  ##draw coverage tracks
  if(coverageTracks==TRUE){
    ##exclude crazy outliers when doing height    
    z=c(tumorWindows[,2],normalWindows[,2])
    height = sort(z)[round(length(z)*0.99)]

    
    ##spacing
    par(mar=c(1,4,2,3),new=TRUE)
    plot(-110,-110, xaxt="n", yaxt="n", ylab="", xlab="",xlim=c(0,100), ylim=c(0,100), bty="n")

    ##normal
    plot(0,0,xlim=xlim,ylim=c(0,height), pch=".", ylab="Normal Covg", xlab="", xaxt="n", cex.lab=0.8, cex.axis=0.8)
    ##color gap regions
    if(!is.null(gaps)){
      drawGaps(gaps, chr)
    }
    ##grid
    abline(v=seq(xaxp[1], xaxp[2], (xaxp[2]-xaxp[1])/xaxp[3]), col="grey50", lty=3)
    yaxp=par()$yaxp
    abline(h=seq(yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3]), col="grey50", lty=3)
    #coverage
    size = normalWindows[2,1]-normalWindows[1,1]
    drawSegs(normalWindows[,1], normalWindows[,1]+size, normalWindows[,2], color="darkgreen", type="r", baseline=0)

    ##spacing
    par(mar=c(2,4,1,3),new=TRUE)
    plot(-110,-110, xaxt="n", yaxt="n", ylab="", xlab="",xlim=c(0,100), ylim=c(0,100), bty="n")
    
    #tumor 
    plot(0,0,xlim=xlim,ylim=c(0,height), pch=".", ylab="Tumor Covg", xlab="", xaxt="n", cex.lab=0.8, cex.axis=0.8)
    ##color gap regions
    if(!is.null(gaps)){
      drawGaps(gaps, chr)
    }
    ##grid
    abline(v=seq(xaxp[1], xaxp[2], (xaxp[2]-xaxp[1])/xaxp[3]), col="grey50", lty=3)
    yaxp=par()$yaxp
    abline(h=seq(yaxp[1], yaxp[2], (yaxp[2]-yaxp[1])/yaxp[3]), col="grey50", lty=3)
    ##coverage
    size = tumorWindows[2,1]-tumorWindows[1,1]
    drawSegs(tumorWindows[,1], tumorWindows[,1]+size, tumorWindows[,2]/tumorNormalRatio, color="sienna4", type="r", baseline=0)
  }


  if(pdfClose){
    dev = dev.off();
  }

}


##-----------------------------------------
## grab the window counts from the potentially very large
## bam window file (100bp bins = 30M windows * multiple libraries, etc)
getWindows <- function(chr, xmin=NULL, xmax=NULL, header=TRUE, tableName=NULL, windowFile=NULL){

  
  ##subset out just the windows we need for this plot
  grabReads <- function(chr, xmin, xmax, header, tableName, windowFile=NULL){

    #nonsense for scope
    f = c();
    assign(tableName,c())

    #read the file
    if(!(is.null(windowFile))){
      assign(tableName, file(windowFile))
      
      if(is.null(xmin) & is.null(xmax)){
        print(paste("select * from ",tableName," where Chr=\"",chr,"\"",sep=""), file.format = list(header = header, row.names = F,sep="\t"))
        winds = sqldf(paste("select * from ",tableName," where Chr=\"",chr,"\"",sep=""), file.format = list(header = header, row.names = F,sep="\t"))
      } else if(is.null(xmin)){
        winds = sqldf(paste("select * from ",tableName," where Chr=\"",chr,"\" AND Start < ",xmax,sep=""), file.format = list(header = header, row.names = F,sep="\t"))
      } else if(is.null(xmax)){
        winds = sqldf(paste("select * from ",tableName," where Chr=\"",chr,"\" AND Start > ",xmin,sep=""), file.format = list(header = header, row.names = F,sep="\t"))
      } else {
        winds = sqldf(paste("select * from ",tableName," where Chr=\"",chr,"\" AND Start < ",xmax," AND Start > ",xmin,sep=""), file.format = list(header = header, row.names = F,sep="\t"))
      }
    } else { #table exists, so just query
      if(is.null(xmin) & is.null(xmax)){
        winds = sqldf(paste("select * from main.",tableName," where Chr=\"",chr,"\"",sep=""))
      } else if(is.null(xmin)){
        winds = sqldf(paste("select * from main.",tableName," where Chr=\"",chr,"\" AND Start < ",xmax,sep=""))
      } else if(is.null(xmax)){
        winds = sqldf(paste("select * from main.",tableName," where Chr=\"",chr,"\" AND Start > ",xmin,sep=""))
      } else {
        winds = sqldf(paste("select * from main.",tableName," where Chr=\"",chr,"\" AND Start < ",xmax," AND Start > ",xmin,sep=""))
      }      
    }    
  }
  
  closeAllConnections()
  winds = c();
  ##sql table already exists
  if(is.null(windowFile)){
    winds = grabReads(chr, xmin, xmax, header, tableName=tableName)
  ##need to get the file first
  } else {                                          
    winds = grabReads(chr, xmin, xmax, header, tableName=tableName, windowFile=windowFile)
  }

  
  if(length(winds[,1]) < 0){
    return(NULL)
  } else {  ##condense columns
    if(length(winds) > 3){
      return(cbind(winds[,2],rowSums(winds[,3:length(winds)])))
    } else {
      return(cbind(winds[,2],winds[,3]))
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

  print(annos)

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

    print(paste(mid2,(ypos+annos[i,5]),annos[i,4],sep=","))
    text(mid2,(ypos+annos[i,5]),annos[i,4],cex=0.5,font=3)
    lines(c(mid,mid2),c(ptop,(ypos+annos[i,5])*.90))
  }
}


##--------------------------------------------------
## draw points for the paired windows
##
drawWinds <- function(winds){
  points(winds$V1,winds$V2)
}


##--------------------------------------------------
## function to actually draw the segments
##
drawSegs <- function(sts, sps, val, color="black", type="r", lowRes=FALSE, lowResMin, lowResMax, baseline){

  ## do the lowres expansion if specified
  if(lowRes){
    lrpos = which((segs[,3]-segs[,2]) > lowResMin)
    segs[lrpos,] = makeVisible(segs[lrpos,],lowResMax)
  }
  if(type=="r"){
    ## draw the segments
    rect(sts, baseline, sps, val, col=color, lty="blank")
  } else if(type=="l"){
    segments(sts, val, sps, val, col=color, lwd=2.5)
  } else {
    print("type parameter for drawSegs must be \"r\" or \"l\"")
    stop()
  }
}


##----------------------------------------------------
## draw the points for the t/n ratio 
##
drawCorrectedWindows <- function(cwind){
  points(cwind[,1], cwind[,2], cex=1.3, pch=21, col=rgb(0,0,0,0.5), bg=rgb(0,0,0,0.5))
}


##----------------------------------------------------
## draw gaps
##
drawGaps <- function(gaps, chr){
  gaps=gaps[gaps$V1==chr,]
  drawSegs(gaps[,2], gaps[,3], baseline=-1e6, val=1e6, color="grey80", type="r")
}

##----------------------------------------------------
## draw gaps
##
drawHighlights <- function(hlsegs, chr){
  for(i in hlsegs){    
    chrpos = strsplit(i,":")[[1]]
    pos = strsplit(chrpos[2],"-")[[1]]
    if(chrpos[1]==chr){
      abline(v=pos[1],col="darkgoldenrod1")
      abline(v=pos[2],col="darkgoldenrod1")
    }
  }
}


##----------------------------------------------------
## calculate tumor normal ratio
##
getTumorNormalRatio <- function(tumorTableName, normalTableName){
  #get the column names
  t = sqldf(paste("select * FROM main.",tumorTableName," where Chr=22", sep=""))[1,]
  n = sqldf(paste("select * FROM main.",normalTableName," where Chr=22", sep=""))[1,]

  tsum = 0
  for(i in names(t)[3:length(t)]){
    tsum = tsum + as.numeric(sqldf(paste("select SUM(\"",i,"\") FROM main.",tumorTableName,sep="")))
  }

  nsum = 0
  for(i in names(n)[3:length(n)]){
    nsum = nsum + as.numeric(sqldf(paste("select SUM(\"",i,"\") FROM main.",normalTableName,sep="")))
  }

  print(paste("Normal reads:",nsum))
  print(paste("Tumor reads:",tsum))
  print(paste("Tumor/Normal depth ratio is",tsum/nsum))
  return(tsum/nsum)
}


