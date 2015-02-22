

plot_tumor_normal_read_depth <- function(coverage_file,plot_title="Sequence Coverage",xcoord_view=NULL,output_file=NULL,transcript.info=NULL,ylim_max=NULL,arrow=NULL,highlight=NULL) {
    library("ggplot2");
    library("plyr");
    library("grid");
    library("reshape2");
    #library("scales");
    
    mbp <-1e6;  #convert the raw genomic coordinates to mbp
    Ylim_percentage <- 1; #set the percent of y axis to label as ymax
    #xcoord_view=NULL;
    #xcoord_view=c(8400000,9100000);


    #read in coverage data
    cov.data=read.table(coverage_file,header=F,sep="\t");
    names(cov.data) <- c('chr','pos1','pos2','normal','tumor','sample');
    cov.data <- subset(cov.data,normal!=0);
    cov.data$diff <- cov.data$tumor-cov.data$normal;
    cov.data$ratio <- cov.data$tumor/cov.data$normal;
    cov.data$mean_pos<- (cov.data$pos1+cov.data$pos2)/2;
    with(cov.data,data.frame(chr,mean_pos,diff,sample,ratio)) -> cov.diff.data;

    cov.data <- subset(cov.data,select=-c(pos1,pos2,diff));
    melt(cov.data,measure.vars=c('normal','tumor')) -> cov.data;

    #set the ylim so that the maxium value is 70% of the max(Y axis limit)
    #ylim_max <- ceiling( max(cov.data$value)/Ylim_percentage); #set max Y scale 

    subset(cov.data,variable=='normal')-> normal.df;
    normal.df$variable <- factor(normal.df$variable);
    ylim_max_normal <- ceiling( max(normal.df$value)/Ylim_percentage);

    subset(cov.data,variable=='tumor')-> tumor.df;
    tumor.df$variable <- factor(tumor.df$variable);
    ylim_max_tumor <- ceiling( max(tumor.df$value)/Ylim_percentage);
    #pick the bigger of the 2 scales.
    #set the max scale on Y axis for depth
    #default:  ylim_max/Ylim_percentage
    #ylim_max = 40; <--- hardcode the value if you want to set it
    if(is.null(ylim_max)) {
        ylim_max <- max(c(ylim_max_normal,ylim_max_tumor));
    }else {
        #ylim_max <- 44
    }
    #transform genomic coordinates
    normal.df$mean_pos <- normal.df$mean_pos/mbp;
    tumor.df$mean_pos <- tumor.df$mean_pos/mbp;
    cov.diff.data$mean_pos <- cov.diff.data$mean_pos/mbp;
    highlight <- highlight/mbp;
    
    if(!is.null(xcoord_view)) {
        xlim_min <- xcoord_view[1]/mbp;
        xlim_max <- xcoord_view[2]/mbp;
    }else {
        #common x-axis tick mark for all graphs
        xlim_min <- round((min(cov.diff.data$mean_pos)),8);
        xlim_max <- round((max(cov.diff.data$mean_pos)),8);
    }

    p_normal <- get_normalized_reads_plot(normal.df,paste(plot_title,"(Normal)",sep=" "),c(xlim_min,xlim_max),ymax_value=ylim_max,plot_color=colors()[51],transcript_file=transcript.info,arrow_size=arrow,highlight_segment=highlight);
    p_tumor  <- get_normalized_reads_plot(tumor.df,paste(plot_title,"(Tumor)",sep=" "),c(xlim_min,xlim_max),ymax_value=ylim_max,plot_color=colors()[556],transcript_file=transcript.info,arrow_size=arrow,highlight_segment=highlight);
    #p_ratio   <- get_tumor_normal_ratio_plot(cov.diff.data,plot_title,c(xlim_min,xlim_max),transcript_file=transcript.info,arrow_size=arrow);
    p_diff <- get_tumor_normal_diff_plot(cov.diff.data,paste(plot_title,"Differential",sep=" "),c(xlim_min,xlim_max),highlight_segment=highlight);

    #if outputfile is defined, print plots to outputfile
    #else, return ggplot_obj to caller
    if(!is.null(output_file)) {
        pdf(file=output_file,width=12,height=8);
        grid.newpage();
        pushViewport(viewport(layout = grid.layout(3,1)))
        #pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
        print(p_normal,vp=viewport(layout.pos.row=1,layout.pos.col=1));
        print(p_tumor, vp=viewport(layout.pos.row=2,layout.pos.col=1));
        print(p_diff, vp=viewport(layout.pos.row=3,layout.pos.col=1));

        # UNCOMMENT THIS IF STATEMENT BELOW TO NOT PRINT A SECOND PAGE IN THE PDF WITH TUMOR ONLY
        #if (0) {
            grid.newpage();
            pushViewport(viewport(layout = grid.layout(3,1)))
            #pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
            print(p_tumor,vp=viewport(layout.pos.row=1,layout.pos.col=1));
            #print(p_tumor, vp=viewport(layout.pos.row=2,layout.pos.col=1));
            #print(p_diff, vp=viewport(layout.pos.row=3,layout.pos.col=1));
            #print(p_tumor);
        #}


            dev.off();
        }
        else {
            return(list(a=p_normal,b=p_tumor,c=p_diff));
        }

    }

    get_normalized_reads_plot<- function(data_input,plot_title='normalized_reads',xlim,plot_color=colors()[170],ymax_value=NULL,ymax_percentage=0.7,transcript_file=NULL,arrow_size=NULL,highlight_segment=NULL) {

        #determine the xlab breakpoints
        #x_tick_mark <- round((seq(xlim[1],xlim[2],len=6)),2);
        x_tick_mark <- axisTicks(c(xlim[1],xlim[2]),log=FALSE,nint=6);
        #set the ylim so that the maxium Y value is either 70% of the max Y axis or overriden by user specified ymax
        if(is.null(ymax_value)) {
            ylim_max <- ceiling( max(data_input$value)/ymax_percentage); #set max Y scale 
        }else {
            ylim_max <- ymax_value;
            y_tick_mark <- axisTicks(c(0,ylim_max),log=FALSE,nint=3);
            #ylim_max <- 40;
        }

        if(!is.null(highlight_segment)) {
            highlight_segment<-sort(highlight_segment);  #sort in ascending 
            highlight_seg_df<- data.frame(highlight_segment[1],highlight_segment[2],0,ylim_max);
            names(highlight_seg_df) <- c('xmin','xmax','ymin','ymax');
        }
        #ylim_max <- ceiling( max(data_input$value)/ymax_percentage); #set max Y scale 

        ylim_min <- -2;  #default ylim_min is zero
        p1 <- ggplot(data_input);
        if(!is.null(transcript_file)) { #if transcript file is defined, draw it
            transcript_info <- process_transcript_struct(transcript_file,ylim_max,arrow_size);
            ylim_min <- transcript_info$exon$ymax[1]; #override the default value of zero 
            p1 <- p1 + geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=label),data=transcript_info$exon);  #draw exons
            p1 <- p1 + geom_segment(aes(x=xmin,xend=xmax,y=ymin,yend=ymin),data=transcript_info$intron,colour='black'); #draw introns
            p1 <- p1 + geom_segment(aes(x=xmin,xend=xmax,y=ymin,yend=ymin),data=transcript_info$arrow,colour=colors()[261],arrow=arrow(length=unit(0.25,'cm'))); #draw arrow line
            p1 <- p1 + scale_fill_manual('', c('utr_exon'=colors()[556],'cds_exon'= colors()[30]),breaks=c('utr_exon','cds_exon'),labels=c('UTR','CDS'),legend=FALSE);
        }
        p1 <- p1 + geom_area(aes(x=mean_pos,y=value),fill=plot_color);
        if(!is.null(highlight_segment)) {
          p1 <- p1 + geom_rect(aes(xmin=xmin,xmax=xmax, ymin=ymin, ymax=ymax),fill=colors()[11],colour='gray',alpha=0.3,data=highlight_seg_df);
        }
        p1 <- p1 + geom_abline(slope=0,intercept=0,colour='black');
        #p1 <- p1 + scale_y_continuous(name="Read Depth",breaks=c(0,20,40,ylim_max));
        p1 <- p1 + scale_y_continuous(name="Read Depth",breaks=y_tick_mark,limits=c(ylim_min,ylim_max));
        p1 <- p1 + scale_x_continuous(name="MB",breaks=x_tick_mark);
        p1 <- p1 + coord_cartesian(xlim=xlim,ylim=c(ylim_min,ylim_max)); #control the area to be displayed
        #p1 <- p1 + scale_fill_manual('', c('utr_exon'=colors()[30],'cds_exon'= colors()[552]),breaks=c('utr_exon','cds_exon'),labels=c('UTR','CDS'),legend=FALSE);
        p1 <- p1 + opts(title=plot_title);
        p1 <- p1 + opts(plot.title = theme_text(size=14, lineheight=.8, face="bold"),axis.text.y=theme_text(colour='black'),axis.text.x=theme_text(colour='black'),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank(),panel.background=theme_rect(fill=colors()[141]));

        return(p1);

    }

    get_tumor_normal_diff_plot<- function (data_input,plot_title='Tumor-Normal',xlim,ymax=20,highlight_segment=NULL) {

        #determine the xlab breakpoints
        x_tick_mark <- axisTicks(c(xlim[1],xlim[2]),log=FALSE,nint=6);
        y_tick_mark <- axisTicks(c(-ymax,ymax),log=FALSE,nint=4);

        if(!is.null(highlight_segment)) {
            highlight_segment<-sort(highlight_segment);  #sort in ascending 
            highlight_seg_df<- data.frame(highlight_segment[1],highlight_segment[2],-ymax,ymax);
            names(highlight_seg_df) <- c('xmin','xmax','ymin','ymax');
        }

        
        pobj <- ggplot(data_input);
        #p1 <- p1 + geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=label),data=data.exon);  #draw exons
        #p1 <- p1 + geom_segment(aes(x=xmin,xend=xmax,y=ymin,yend=ymin),data=data.intron,colour='black'); #draw introns
        pobj <- pobj + geom_area(aes(x=mean_pos,y=diff),fill=colors()[28]);
        if(!is.null(highlight_segment)) {
          pobj <- pobj + geom_rect(aes(xmin=xmin,xmax=xmax, ymin=ymin, ymax=ymax),fill=colors()[11],colour='gray',alpha=0.3,data=highlight_seg_df);
        }

        # NDees changed the use of ylimits - removing it so that the positive and negative diff would be scaled similar to the other images. There is still some issue whereby the max and min on the differential plot still seems to sometimes be a different number than the max value of the tumor and normal plot. My guess is that it occurred to me when I was using normalization. Perhaps the differences between the normalization and actual heights for tumor coverage can explain the phenomenon.
        #ylimits = round((range(data_input$diff)));
        #pobj <- pobj + scale_y_continuous(name="Tumor-Normal",breaks=c(0),limits=c(ylimits[1],ylimits[2]));

        pobj <- pobj + scale_y_continuous(name="Tumor-Normal",breaks=y_tick_mark,limits=c(-ymax,ymax));
        pobj <- pobj + scale_x_continuous(name="MB",breaks=x_tick_mark);
        pobj <- pobj + coord_cartesian(xlim=xlim,ylim=c(-ymax,ymax)); #control the area to be displayed
        #p2 <- p2 + scale_fill_manual('', c('utr_exon'=colors()[30],'cds_exon'= colors()[50]),breaks=c('utr_exon','cds_exon'),labels=c('UTR','CDS'),legend=FALSE);
        pobj <- pobj + opts(title=plot_title);
        pobj <- pobj + opts(plot.title = theme_text(size=14, lineheight=.8, face="bold"),axis.text.y=theme_text(colour='black'),axis.text.x=theme_text(colour='black'),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank(),panel.background=theme_rect(fill=colors()[141]));

        return(pobj);

    }

    get_tumor_normal_ratio_plot<-function(data_input,plot_title="Tumor/Normal Coverage",xlim,transcript_file=NULL,arrow_size=NULL) {

        #determine the xlab breakpoints
        x_tick_mark <- round((seq(xlim[1],xlim[2],len=5)),2);

        ylim_max <- 2.5;
        ylim_min <- 0;
        pobj <- ggplot(data_input);
        if(!is.null(transcript_file)) { #if transcript file is defined, draw it
            transcript_info <- process_transcript_struct(transcript_file,ylim_max,arrow_size);
            ylim_min <- transcript_info$exon$ymax[1]; #override the default value of zero 
            pobj <- pobj + geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='black',data=transcript_info$exon);  #draw exons
            pobj <- pobj + geom_segment(aes(x=xmin,xend=xmax,y=ymin,yend=ymin),data=transcript_info$intron,colour='black'); #draw introns
            pobj <- pobj + geom_segment(aes(x=xmin,xend=xmax,y=ymin,yend=ymin),data=transcript_info$arrow,colour='black',arrow=arrow(length=unit(0.25,'cm'))); #draw arrow line
            pobj <- pobj + scale_fill_manual('', c('utr_exon'=colors()[30],'cds_exon'= colors()[552]),breaks=c('utr_exon','cds_exon'),labels=c('UTR','CDS'),legend=FALSE);
        }
        #ylimits = round((range(data_input$ratio)));
        pobj <- pobj + geom_hline(yintercept=1,size=0.5,colour='grey72',linetype=2);
        pobj <- pobj + geom_hline(yintercept=0,size=0.5,colour='black');
        pobj <- pobj + geom_point(aes(x=mean_pos,y=ratio),colour=colors()[28],size=1.5);
        pobj <- pobj + scale_y_continuous(name="Tumor/Normal",breaks=c(0,1,2,4),limits=c(ylim_min,ylim_max));
        #pobj <- pobj + scale_y_continuous(name="Tumor/Normal",breaks=c(0,0.5,1,ceiling(ylimits[2])),limits=c(0,ceiling(ylimits[2])));
        pobj <- pobj + scale_x_continuous(name="MB",breaks=x_tick_mark);
        pobj <- pobj + coord_cartesian(xlim=xlim); #control the area to be displayed

        pobj <- pobj + opts(title=plot_title);
        pobj <- pobj + opts(plot.title = theme_text(size=14, lineheight=.8, face="bold"),axis.text.y=theme_text(colour='black'),axis.text.x=theme_text(colour='black'),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank(),panel.background=theme_rect(fill=colors()[141]));

        return(pobj);


    }



    get_tumor_normal_ratio_plot2<-function(data_input,plot_title="Tumor/Normal Coverage",xlim,transcript_file=NULL) {

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
        ylimits = round((range(data_input$ratio)));
        pobj <- pobj + geom_hline(yintercept=1,size=0.5,colour='purple');
        pobj <- pobj + geom_hline(yintercept=0,size=0.5,colour='black');
        #pobj <- pobj + scale_y_continuous(name="Tumor/Normal",breaks=c(0,1,2,4),limits=c(0,4));
        pobj <- pobj + scale_y_continuous(name="Tumor/Normal",breaks=c(0,0.5,1,ceiling(ylimits[2])),limits=c(-5,ceiling(ylimits[2])));
        pobj <- pobj + scale_x_continuous(name="MB",breaks=x_tick_mark);
        pobj <- pobj + coord_cartesian(xlim=xlim); #control the area to be displayed

        pobj <- pobj + opts(title=plot_title);
        pobj <- pobj + opts(plot.title = theme_text(size=14, lineheight=.8, face="bold"),axis.text.y=theme_text(colour='black'),axis.text.x=theme_text(colour='black'),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank(),panel.background=theme_rect(fill=colors()[141]));

        return(pobj);


    }


    process_transcript_struct <- function(transcript_file,y_max,arrow_size=0.00010) {

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



    plot_gene_coverage <- function (gene_coordinate_file,coverage_file,title="Depth for Gene",depth_limit=4,xcoord_view=NULL,xcoord_buffer=50000) {

        mbp <-1e6;
        #print(mbp);
        #read in transcript structure data
        #data.in = read.table("test.atrx.gene.txt",sep="\t",header=T);
        data.in = read.table(gene_coordinate_file,sep="\t",header=T);
        data.in$xmin = data.in$xmin/mbp;  #transform the coordinates to megabases
        data.in$xmax = data.in$xmax/mbp;
        strand <- data.in$strand[1];


        #read in coverage data
        #cov.data <- read.table("ATRX_intron_exon.depth",sep="\t",header=F);
        cov.data=read.table(coverage_file,header=F,sep="\t");
        cov.data$ratio <- cov.data$V6/cov.data$V4;
        cov.data$V2 <- cov.data$V2/mbp;
        cov.data2 <- subset(cov.data,V3 > 29);  #require at least 30 reads for normal sample

        #separate intron and exons
        data.intron = subset(data.in,grepl("intron",label));
        data.exon = subset(data.in,grepl("exon",label));
        data.exon$label <- factor(data.exon$label);
        data.intron$avg <- (data.intron$xmin+data.intron$xmax)/2;
        num_exon = nrow(data.exon);
        num_intron = nrow(data.intron);
        ylim_max <- depth_limit;
        ylim_min <- -ylim_max/10;
        introns  <- ylim_min/2;
        exon_min <- 3/4*ylim_min;
        exon_max <- 1/4*ylim_min;
        data.exon$ymin <- rep(exon_min,num_exon);
        data.exon$ymax <- rep(exon_max,num_exon);
        data.intron$ymin <- rep(introns,num_intron);
        #determine the direction of the transcript
        if(strand == '+') { #+strand
            endOfGene <- max(data.in$xmax);
            beginOfGene <- min(data.in$xmin);
            xlim_min <- beginOfGene-xcoord_buffer/mbp;
            xlim_max <- endOfGene+xcoord_buffer/mbp;
            line.arrow <- data.frame(xmin=endOfGene,xmax=endOfGene+(xcoord_buffer/2)/mbp,ymin=introns);
            #p <- p + geom_segment(aes(x=xmin,xend=avg,y=ymin,yend=ymin),data=data.intron,colour='black',arrow=arrow(length=unit(0.25,'cm')));
        }else { #-strand
            endOfGene <- min(data.in$xmin);
            beginOfGene <- max(data.in$xmax);
            xlim_min <- endOfGene-xcoord_buffer/mbp;
            xlim_max <- beginOfGene+xcoord_buffer/mbp;
            line.arrow <- data.frame(xmin=endOfGene,xmax=endOfGene-(xcoord_buffer/2)/mbp,ymin=introns);
            #p <- p + geom_segment(aes(x=xmax,xend=avg,y=ymin,yend=ymin),data=data.intron,colour='black',arrow=arrow(length=unit(0.25,'cm')));
        }

        if(!is.null(xcoord_view)) {
            xlim_min <- xcoord_view[1]/mbp;
            xlim_max <- xcoord_view[2]/mbp;
        }
        else {
            #p <- p + coord_cartesian(xlim=xlimit,ylim=ylimit); #default area of gene to display (whole gene)
        }

        #plot_obj <- make_area_plot_1sample(data.in,data.exon,data.intron,line.arrow,coverage_data=cov.data,xlimit=c(xlim_min,xlim_max),ylimit=c(ylim_min,ylim_max+.01*ylim_max),plot_title=title,);
        p <- ggplot(cov.data2);
        p <- p + geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=label),data=data.exon);  #draw exons
        p <- p + geom_segment(aes(x=xmin,xend=xmax,y=ymin,yend=ymin),data=data.intron,colour='black'); #draw introns
        p <- p + geom_segment(aes(x=xmin,xend=xmax,y=ymin,yend=ymin),data=line.arrow,colour='black',arrow=arrow(length=unit(0.25,'cm'))); #draw arrow line

        p <- p + geom_point(aes(x=V2,y=ratio),colour=colors()[170],size=0.5) + facet_grid(V7 ~ ., scales="free_y"); #plot coverage
        #p <- p + geom_area(aes(x=V2,y=ratio),data=cov.data,fill='grey29'); #plot coverage
        p <- p + geom_hline(yintercept=0.5,size=0.5,colour='blue4');
        p <- p + geom_hline(yintercept=1.5,size=0.5,colour='blue4');
        p <- p + scale_y_continuous(name="Tumor/Normal",breaks=c(0,1,2,4,6));
        p <- p + scale_x_continuous(name="MB");
        p <- p + coord_cartesian(xlim=c(xlim_min,xlim_max),ylim=c(ylim_min,ylim_max)); #control the area to be displayed
        p <- p + scale_fill_manual('', c('utr_exon'=colors()[30],'cds_exon'= colors()[35]),breaks=c('utr_exon','cds_exon'),labels=c('UTR','CDS'),legend=FALSE);
        p <- p + coord_cartesian(xlim=c(xlim_min,xlim_max),ylim=c(ylim_min,ylim_max));

        p <- p + opts(title=title, plot.title = theme_text(size=14, lineheight=.8, face="bold"),axis.text.y=theme_text(colour='black'),axis.text.x=theme_text(colour='black'),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank(),panel.background=theme_rect(fill=colors()[141])); 


        return(p);

    }

