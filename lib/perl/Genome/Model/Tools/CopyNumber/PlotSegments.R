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
                         highlights=NULL, lowRes=FALSE, lowResMin=NULL,
                         lowResMax=NULL, showNorm=FALSE, baseline=2,
                         gainThresh=2.5, lossThresh=1.5, annotationsTop=NULL,
                         annotationsBottom = NULL, plotTitle="",
                         gainColor="red", lossColor="blue", ylabel="",
                         xmin=NULL, xmax=NULL, label_size=0.6){

  ## add options for plotting just a smaller region - TODO
  xlim = NULL

  ## read in the segments
  segs=read.table(filename,comment.char="#")

  ## read in the entrypoints
  entrypoints=addOffsets(readEntrypoints(entrypoints))
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
  for(i in 1:length(chrnames)){
    #raise an error if entrypoints and chrs don't match
    if(length(which(entrypoints$chr==chrnames[i])) < 1){
      cat("\nERROR - no entrypoint found for chromosome ",chrnames[i]," found in segs file\n")
      cat("maybe you meant to use the male entrypoints?\n")
      q(save="no")
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
      #print(annos)
      #print(annos[i,5])
      print(paste(mid2,(ypos+annos[i,5]),annos[i,4],sep=","))
      text(mid2,(ypos+annos[i,5]),annos[i,4],cex=0.5,font=3)
      lines(c(mid,mid2),c(ptop,(ypos+annos[i,5])*.90))
    }
  }

  
  
################################################
  ## plot all chromosomes
  if(chr=="ALL"){

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
    plot(0, 0, xlim=xlim, ylim=c(ymin,ymax), pch=".",
         ylab=ylabel, xlab="", xaxt="n", cex.lab=1, cex.axis=label_size)

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
    abline(h=baseline,col="grey25")

    offsets = as.numeric(entrypoints[,4])
    offsets = append(offsets,sum(entrypoints$length))

    
    ## draw highlight regions, if specified
    if(!(is.null(highlights))){      
      chrNames = names(table(hlRegions[,1]))
      for(i in 1:length(chrNames)){
        #don't plot things that we aren't considering
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
      rect(segs[,2], baseline, segs[,3], segs[,5], col=color,lty="blank")
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
    abline(v=0,col="gray75")
    for(i in 1:(length(offsets)-1)){
      abline(v=offsets[i+1],col="gray75")
      text((offsets[i]+offsets[i+1])/2, ymax*0.9, labels= gsub("chr","",entrypoints[i,1]), cex=label_size)
    }


    ## add top annotations, if specfied
    if(!(is.null(annotationsTop))){
      addAnnos(annTopRegions, segs, top=TRUE, offset=TRUE)
    }
    ## add btm annotations, if specfied
    if(!(is.null(annotationsBottom))){
      addAnnos(annBtmRegions, segs, top=FALSE, offset=TRUE)
    }

    
    ## ## add top annotations, if specfied
    ## if(!(is.null(annotationsTop))){
    ##   ypos = ymax*0.8
    ##   for(i in 1:length(annTopRegions[,1])){
    ##     offset=as.numeric(entrypoints[which(entrypoints$chr==annTopRegions[i,1]),4])
    ##     st=as.numeric(annTopRegions[i,2])+offset
    ##     sp=as.numeric(annTopRegions[i,3])+offset
    ##     mid=((sp-st)/2)+st
    ##     ##get the height of the peak at this position (if it exists)
    ##     ptop = 0
    ##     peakNum=which(segs[,1] == as.numeric(annTopRegions[i,1]) &
    ##     (segs[,2] <= mid-offset) & (segs[,3] >= mid-offset))

    ##     if(length(peakNum > 0)){
    ##       ptop = max(segs[peakNum,5])+((ymax-as.numeric(baseline))*0.10)
    ##     }

    ##     ##adjust at edge of plot
    ##     mid2 = mid
    ##     if(mid <50000000){
    ##       mid2 = 50000000
    ##     }
    ##     text(mid2,ypos,annTopRegions[i,4],cex=0.5,font=3)
    ##     lines(c(mid,mid2),c(ptop,ypos*.95))
    ##   }
    ## }

    ## ## add bottom annotations, if specfied
    ## if(!(is.null(annotationsBottom))){
    ##   ypos = ymin*0.85
    ##   for(i in 1:length(annBtmRegions[,1])){
    ##     offset=as.numeric(entrypoints[which(entrypoints$chr==annBtmRegions[i,1]),4])
    ##     st=as.numeric(annBtmRegions[i,2])+offset
    ##     sp=as.numeric(annBtmRegions[i,3])+offset
    ##     mid=(sp-st)+st
    ##     ##get the height of the peak at this position (if it exists)
    ##     ptop = 0
    ##     peakNum=which((segs[,1] == as.numeric(annBtmRegions[i,1])) &
    ##       (segs[,2] <= mid-offset) & (segs[,3] >= mid-offset))

    ##     if(length(peakNum > 0)){
    ##       ptop = min(segs[peakNum,5])+((ymin-as.numeric(baseline))*0.10)
    ##     }

    ##     ##adjust at edge of plot
    ##     mid2 = mid
    ##     if(mid <55000000){
    ##       mid2 = 55000000
    ##     }
    ##     text(mid2,ypos,annBtmRegions[i,4],cex=0.5,font=3)
    ##     lines(c(mid,mid2),c(ptop,ypos*.90))
    ##   }
    ## }



############################################################
    ## --------single chromosome-----------------
  } else { #chr != "ALL"

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
    plot(0,0,xlim=xlim,ylim=c(ymin,ymax),pch=".",ylab=ylabel, xlab="position (Mb)",xaxt="n",cex.lab=0.8, cex.axis=0.8)

    ##add the title
    if(!(is.null(plotTitle))){
      title(main=plotTitle)
    }

    ## draw baselines
    abline(h=baseline,col="grey25")

    ## set the x-axis labels
    axis(1, at=seq(0,entry$length,5e6), cex.axis=0.8)

    ## draw highlight regions if specified
    if(!(is.null(highlights))){
      hlRegions = hlRegions[which(hlRegions[,1] == chr),]
      if(length(hlRegions[,1]) > 0){
        for(i in 1:length(hlRegions[,1])){
          if(lowRes){
            lrpos = which((hlRegions[,3]-hlRegions[,2]) > lowResMin)
            hlRegions[lrpos,] = makeVisible(hlRegions[lrpos,],lowResMax)
          }
          
          rect(hlRegions[,2], ymin*2, hlRegions[,3], ymax*2, col="gold",lwd=0,lty="blank")
        }
      }
    }

    ## function to actually draw the segments
    drawSegs <- function(segs,color="black"){

      ## do the lowres expansion if specified
      if(lowRes){
        lrpos = which((segs[,3]-segs[,2]) > lowResMin)
        segs[lrpos,] = makeVisible(segs[lrpos,],lowResMax)
      }

      ## draw the segments
      rect(segs[,2], baseline, segs[,3], segs[,5], col=color,lty="blank")
    }


    ## do the plotting

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
    
    ## add top annotations, if specfied
    if(!(is.null(annotationsTop))){
      addAnnos(annTopRegions, segs, top=TRUE, chr=chr, leftEdge=5000000)
    }
    ## add btm annotations, if specfied
    if(!(is.null(annotationsBottom))){
      addAnnos(annBtmRegions, segs, top=FALSE, chr=chr, leftEdge=5000000)
    }
  }
}
