#This file contains code to graph a small region of interest from BamToCna data

graph_cna_region=function(cna_file="",output_file="",chromosome="",region_start="",region_end="",roi_start="",roi_end="",title="",xlab="",ylab="",color_roi=T) {

    png(output_file,width=5,height=6,res=300,units="in",type="cairo");
    x<-read.table(cna_file,header=TRUE);    #this assumes a normal bamtocna output file
    x_region <- x[x$CHR==chromosome & x$POS >= region_start & x$POS <= region_end,];
    x_region_roi <- x_region[x_region$POS >= roi_start & x_region$POS <= roi_end,];
    x_region_nonroi <- x_region[x_region$POS < roi_start | x_region$POS > roi_end,];
    roi_color="blue";
    if(color_roi) {
        roi_color="red";
    }
    plot(x_region_nonroi$POS/1000000,x_region_nonroi$DIFF+2,ylim=c(0,6),col="blue",cex=0.5,pch=19,ylab=ylab,xlab=xlab,main=title);
    points(x_region_roi$POS/1000000,x_region_roi$DIFF+2,ylim=c(0,6),col=roi_color,cex=0.5,pch=19);
    lines(c(roi_start,roi_end)/1000000,c(median(x_region_roi$DIFF+2),median(x_region_roi$DIFF+2)),col="green",lwd=4);
    dev.off();
}


#code below plots normalized coverage from a ROI 
plot_tumor_normal_read_depth <- function(coverage_file,plot_title="Sequence Coverage",xcoord_view=NULL,output_file=NULL,transcript.info=NULL) {

  library("ggplot2");
  mbp <-1e6;  #convert the raw genomic coordinates to mbp
  Ylim_percentage <- 0.7; #set the percent of y axis to label as ymax
  #xcoord_view=NULL;
  #xcoord_view=c(8400000,9100000);


  #read in coverage data
  cov.data=read.table(coverage_file,header=F,sep="\t");
  names(cov.data) <- c('chr','pos1','pos2','normal','tumor','sample');
  cov.data$diff <- cov.data$tumor-cov.data$normal;
  cov.data$ratio <- cov.data$tumor/cov.data$normal;
  cov.data$mean_pos<- (cov.data$pos1+cov.data$pos2)/2;
  with(cov.data,data.frame(chr,mean_pos,diff,sample,ratio)) -> cov.diff.data;
  cov.data <- subset(cov.data,select=-c(pos1,pos2,diff));
  melt(cov.data,measure.vars=c('normal','tumor')) -> cov.data;

  #set the ylim so that the maxium value is 70% of the max(Y axis limit)
  ylim_max <- ceiling( max(cov.data$value)/Ylim_percentage); #set max Y scale 

  subset(cov.data,variable=='normal')-> normal.df;
  normal.df$variable <- factor(normal.df$variable);
  subset(cov.data,variable=='tumor')-> tumor.df;
  tumor.df$variable <- factor(tumor.df$variable);
  #transform genomic coordinates
  normal.df$mean_pos <- normal.df$mean_pos/mbp;
  tumor.df$mean_pos <- tumor.df$mean_pos/mbp;
  cov.diff.data$mean_pos <- cov.diff.data$mean_pos/mbp;

  if(!is.null(xcoord_view)) {
    xlim_min <- xcoord_view[1]/mbp;
    xlim_max <- xcoord_view[2]/mbp;
  }else {
    #common x-axis tick mark for all graphs
    xlim_min <- round((min(cov.diff.data$mean_pos)),2);
    xlim_max <- round((max(cov.diff.data$mean_pos)),2);
  }
  x_tick_mark <- round((seq(xlim_min,xlim_max,len=4)),2);

  p_normal <- get_normalized_reads_plot(normal.df,paste(plot_title,"(normal)",sep=" "),c(xlim_min,xlim_max),plot_color=colors()[51],transcript_file=transcript.info);
  p_tumor  <- get_normalized_reads_plot(tumor.df,paste(plot_title,"(tumor)",sep=" "),c(xlim_min,xlim_max),plot_color=colors()[556],transcript_file=transcript.info);
  #p_diff   <- get_tumor_normal_diff_plot(cov.diff.data,plot_title,c(xlim_min,xlim_max));
  p_ratio   <- get_tumor_normal_ratio_plot(cov.diff.data,plot_title,c(xlim_min,xlim_max),transcript_file=transcript.info);

  #if outputfile is defined, print plots to outputfile
  #else, return ggplot_obj to caller
  if(!is.null(output_file)) {
    pdf(file=output_file,width=12,height=8);
    grid.newpage();
    pushViewport(viewport(layout = grid.layout(3,1)))
    #pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    print(p_normal,vp=viewport(layout.pos.row=1,layout.pos.col=1));
    print(p_tumor, vp=viewport(layout.pos.row=2,layout.pos.col=1));
    print(p_ratio,  vp=viewport(layout.pos.row=3,layout.pos.col=1));
    dev.off();
  }
  else {
    return(list(a=p_normal,b=p_tumor,c=p_diff));
  }

}

get_normalized_reads_plot<- function(data_input,plot_title='normalized_reads',xlim,plot_color=colors()[170],ymax_percentage=0.7,transcript_file=NULL) {

  #determine the xlab breakpoints
  x_tick_mark <- round((seq(xlim[1],xlim[2],len=4)),2);
  #set the ylim so that the maxium value is 60% of the Y axis
  ylim_max <- ceiling( max(data_input$value)/ymax_percentage); #set max Y scale 

  ylim_min <- 0;  #default ylim_min is zero
  p1 <- ggplot(data_input);
  if(!is.null(transcript_file)) { #if transcript file is defined, draw it
    transcript_info <- process_transcript_struct(transcript_file,ylim_max);
    ylim_min <- transcript_info$exon$ymax[1]; #override the default value of zero 
    p1 <- p1 + geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='black',data=transcript_info$exon);  #draw exons
    p1 <- p1 + geom_segment(aes(x=xmin,xend=xmax,y=ymin,yend=ymin),data=transcript_info$intron,colour='black'); #draw introns
    p1 <- p1 + geom_segment(aes(x=xmin,xend=xmax,y=ymin,yend=ymin),data=transcript_info$arrow,colour='black',arrow=arrow(length=unit(0.25,'cm'))); #draw arrow line
    p1 <- p1 + scale_fill_manual('', c('utr_exon'=colors()[30],'cds_exon'= colors()[552]),breaks=c('utr_exon','cds_exon'),labels=c('UTR','CDS'),legend=FALSE);
  }
  p1 <- p1 + geom_area(aes(x=mean_pos,y=value),fill=plot_color);
  p1 <- p1 + scale_y_continuous(name="Read Depth",breaks=c(0,40,80),limits=c(ylim_min,ylim_max));
  p1 <- p1 + scale_x_continuous(name="MB",breaks=x_tick_mark);
  p1 <- p1 + coord_cartesian(xlim=xlim); #control the area to be displayed
  #p1 <- p1 + scale_fill_manual('', c('utr_exon'=colors()[30],'cds_exon'= colors()[552]),breaks=c('utr_exon','cds_exon'),labels=c('UTR','CDS'),legend=FALSE);
  p1 <- p1 + opts(title=plot_title);
  p1 <- p1 + opts(plot.title = theme_text(size=14, lineheight=.8, face="bold"),axis.text.y=theme_text(colour='black'),axis.text.x=theme_text(colour='black'),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank(),panel.background=theme_rect(fill=colors()[141]));

  return(p1);
  
}

get_tumor_normal_diff_plot<- function (data_input,plot_title='Tumor-Normal',xlim) {

  #determine the xlab breakpoints
  x_tick_mark <- round((seq(xlim[1],xlim[2],len=4)),2);

  pobj <- ggplot(data_input);
  #p1 <- p1 + geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=label),data=data.exon);  #draw exons
  #p1 <- p1 + geom_segment(aes(x=xmin,xend=xmax,y=ymin,yend=ymin),data=data.intron,colour='black'); #draw introns
  pobj <- pobj + geom_area(aes(x=mean_pos,y=diff),fill=colors()[28]);
  ylimits = round((range(data_input$diff)));
  pobj <- pobj + scale_y_continuous(name="Tumor-Normal",breaks=c(ylimits[1],0,ylimits[2]),limits=c(ylimits[1],ylimits[2]));
  pobj <- pobj + scale_x_continuous(name="MB",breaks=x_tick_mark);
  pobj <- pobj + coord_cartesian(xlim=xlim); #control the area to be displayed
  #p2 <- p2 + scale_fill_manual('', c('utr_exon'=colors()[30],'cds_exon'= colors()[50]),breaks=c('utr_exon','cds_exon'),labels=c('UTR','CDS'),legend=FALSE);
  pobj <- pobj + opts(title=plot_title);
  pobj <- pobj + opts(plot.title = theme_text(size=14, lineheight=.8, face="bold"),axis.text.y=theme_text(colour='black'),axis.text.x=theme_text(colour='black'),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank(),panel.background=theme_rect(fill=colors()[141]));

  return(pobj);
  
}

get_tumor_normal_ratio_plot<-function(data_input,plot_title="Tumor/Normal Coverage",xlim,transcript_file=NULL) {

  #determine the xlab breakpoints
  x_tick_mark <- round((seq(xlim[1],xlim[2],len=4)),2);

  ylim_max <- 4;
  ylim_min <- 0;
  pobj <- ggplot(data_input);
  if(!is.null(transcript_file)) { #if transcript file is defined, draw it
    transcript_info <- process_transcript_struct(transcript_file,ylim_max);
    ylim_min <- transcript_info$exon$ymax[1]; #override the default value of zero 
    pobj <- pobj + geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='black',data=transcript_info$exon);  #draw exons
    pobj <- pobj + geom_segment(aes(x=xmin,xend=xmax,y=ymin,yend=ymin),data=transcript_info$intron,colour='black'); #draw introns
    pobj <- pobj + geom_segment(aes(x=xmin,xend=xmax,y=ymin,yend=ymin),data=transcript_info$arrow,colour='black',arrow=arrow(length=unit(0.25,'cm'))); #draw arrow line
    pobj <- pobj + scale_fill_manual('', c('utr_exon'=colors()[30],'cds_exon'= colors()[552]),breaks=c('utr_exon','cds_exon'),labels=c('UTR','CDS'),legend=FALSE);
  }
  pobj <- pobj + geom_line(aes(x=mean_pos,y=ratio),colour=colors()[28]);
  #ylimits = round((range(data_input$ratio)));
  pobj <- pobj + geom_hline(yintercept=1,size=0.5,colour='purple');
  pobj <- pobj + geom_hline(yintercept=0,size=0.5,colour='black');
  pobj <- pobj + scale_y_continuous(name="Tumor/Normal",breaks=c(0,1,2,4),limits=c(ylim_min,ylim_max));
  #pobj <- pobj + scale_y_continuous(name="Tumor/Normal",breaks=c(0,0.5,1,ceiling(ylimits[2])),limits=c(0,ceiling(ylimits[2])));
  pobj <- pobj + scale_x_continuous(name="MB",breaks=x_tick_mark);
  pobj <- pobj + coord_cartesian(xlim=xlim); #control the area to be displayed

  pobj <- pobj + opts(title=plot_title);
  pobj <- pobj + opts(plot.title = theme_text(size=14, lineheight=.8, face="bold"),axis.text.y=theme_text(colour='black'),axis.text.x=theme_text(colour='black'),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank(),panel.background=theme_rect(fill=colors()[141]));

  return(pobj);
  
  
}


process_transcript_struct <- function(transcript_file,y_max,arrow_size=0.01) {

  mbp <-1e6;
  #read in transcript structure data
  data.in= read.table(transcript_file,sep="\t",header=T);
  #data.in = read.table("test.alk.gene.txt",sep="\t",header=T);
  data.in$xmin = data.in$xmin/mbp;  #transform the coordinates to megabases
  data.in$xmax = data.in$xmax/mbp;
  strand <- data.in$strand[1];

  data.intron = subset(data.in,grepl("intron",label));
  data.exon = subset(data.in,grepl("exon",label));
  data.exon$label <- factor(data.exon$label);
  data.intron$avg <- (data.intron$xmin+data.intron$xmax)/2;
  num_exon = nrow(data.exon);
  num_intron = nrow(data.intron);
  
  #y_min <- -y_max/10;
  #introns  <- y_min/2;
  #exon_max <- 1/4*y_min;
  #exon_min <- 3/4*y_min;
  ##########################
  introns <- -y_max/10;
  exon_min <- 1/2*introns;
  exon_max <- exon_min*3;
  y_min <- exon_max;
  ##########################

  data.exon$ymin <- rep(exon_min,num_exon);
  data.exon$ymax <- rep(exon_max,num_exon);
  data.intron$ymin <- rep(introns,num_intron);
  #determine the direction of the transcript
  if(strand == '+') { #+strand
    endOfGene   <- max(data.in$xmax);
    beginOfGene <- min(data.in$xmin);
    line.arrow  <- data.frame(xmin=endOfGene,xmax=endOfGene+(arrow_size*endOfGene),ymin=introns);
    #p <- p + geom_segment(aes(x=xmin,xend=avg,y=ymin,yend=ymin),data=data.intron,colour='black',arrow=arrow(length=unit(0.25,'cm')));
  }else { #-strand
    endOfGene   <- min(data.in$xmin);
    beginOfGene <- max(data.in$xmax);
    line.arrow  <- data.frame(xmin=endOfGene,xmax=endOfGene-(arrow_size*endOfGene),ymin=introns);
    #p <- p + geom_segment(aes(x=xmax,xend=avg,y=ymin,yend=ymin),data=data.intron,colour='black',arrow=arrow(length=unit(0.25,'cm')));
  }

  return(list(exon=data.exon,intron=data.intron,arrow=line.arrow));

  
}


