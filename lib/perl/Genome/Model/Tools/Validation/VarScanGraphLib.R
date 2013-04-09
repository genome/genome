#The following are functions for graphing density plots for varscan
#All of this is fragile to changes in varscan file format

#This graphs the density curves of the difference in tumor and normal variant frequencies. If type is all (the default) then the frequency is graphed for all sites in x
#If the type is category then 3 curves are plotted color coded by their VarScan called class (reference Germline Somatic only)
#If a filename is provided then a pdf is written. The .pdf extension is not assigned by default. By default, the function uses the screen to plot.
varscan.graph_frequency_difference_density = function(x,genome="",type = "all", file="") {
    if(file != "") {
        pdf(file=file,width=10,height=7.5);
    }
    if(type=="all") {
        plot(density(x$V20),xlab="Allele Frequency Difference (Tumor - Normal)",ylab="Kernel Density", main=paste(genome," Allele Frequency Difference",sep=""),xlim=c(0,100));
    }
    if(type=="category") {
        y=subset(x, x$V13 == 'Reference');
        a=subset(x, x$V13 == 'Germline');
        z=subset(x, x$V13 == 'Somatic');
        b=subset(x, x$V13 != 'Reference' & x$V13 != 'Germline' & x$V13 != 'Somatic');
        y_density = density(y$V20);
        a_density = density(a$V20);
        z_density = density(z$V20);
        b_density = density(b$V20);

        plot(y_density,xlim=c(0,100),xlab="Allele Frequency Difference (Tumor - Normal)",ylab="Kernel Density Estimation", main=paste(genome," Allele Frequency Difference",sep=""));
        lines(a_density,col="blue");
        lines(z_density,col="red");
        lines(b_density,col="green");
        legend(x="topright",c(paste("Reference (bw=",y_density$bw,")"),paste("Germline (bw=",a_density$bw,")"),paste("Somatic (bw=",z_density$bw,")"),paste("Other (bw=",b_density$bw,")")),col=c("black","blue","red","green"),lwd=1);
    }
    if(file != "") {
        dev.off();
    }
}

#This graphs a scatter of the variant frequencies in tumor and normal. The tumor is the y-axis and the normal the x-axis.
#If type is all then the frequency is graphed for all sites in x and not color coded
#If the type is category (the default) then points are color coded by their VarScan called class (Reference, Germline, Somatic only)
#If a filename is provided then a pdf is written. The .pdf extension is not assigned by default. By default, the function uses the screen to plot.
varscan.graph_frequency_scatter = function(x,genome="",type = "category", file="",xtissue="Normal",ytissue="Tumor",width=10, height=7.5) {
    ylabel = paste(ytissue,"Variant Allele Frequency");
    xlabel = paste(xtissue,"Variant Allele Frequency");
    if(file != "") {
        pdf(file=file,width=width,height=height);
    }
    if(type=="all") {
        plot.default(x=x$V7,y=x$V11,xlab=xlabel,ylab=ylabel, cex.axis=1.25, cex.lab=1.25,  family="sans", font.lab=2, main=paste(genome," Allele Frequency",sep=""), type="p",pch=19,cex=0.4, col="black",xlim=c(0,100),ylim=c(0,100))   
    }
    if(type=="category") {
        par(c(1,2));
        y=subset(x, x$V13 == 'Reference');
        a=subset(x, x$V13 == 'Germline');
        z=subset(x, x$V13 == 'Somatic');
        b=subset(x, x$V13 != 'Reference' & x$V13 != 'Germline' & x$V13 != 'Somatic');
        plot.default(x=y$V7,y=y$V11,xlab=xlabel,ylab=ylabel, main=paste(genome," Allele Frequency",sep=""), type="p",pch=19,cex=0.4, col="black", cex.axis=1.25, cex.lab=1.25,  family="sans", font.lab=2, xlim=c(0,100),ylim=c(0,100));
        points(x=a$V7,y=a$V11, type="p",pch=19,cex=0.4, col="blue",xlim=c(0,100),ylim=c(0,100));
        points(x=z$V7,y=z$V11, type="p",pch=19,cex=0.4, col="red",xlim=c(0,100),ylim=c(0,100));
        points(x=b$V7,y=b$V11, type="p",pch=19,cex=0.4, col="green",xlim=c(0,100),ylim=c(0,100));
        legend(x="bottomright",c("Reference","Germline","Somatic","Other"),col=c("black","blue","red","green"),pch=19); #bottom right will rarely have any points
    }
    if(type=="somatic") {
        par(c(1,2));
        z=subset(x, x$V13 == 'Somatic');
        plot.default(x=z$V7,y=z$V11,xlab=xlabel,ylab=ylabel, main=paste(genome," Allele Frequency",sep=""), type="p",pch=19,cex=0.4, col="red",xlim=c(0,100),ylim=c(0,100),cex.axis=1.25, cex.lab=1.25,  family="sans", font.lab=2);
        legend(x="bottomright",c("Somatic"),col=c("red"),pch=19); #bottom right will rarely have any points
    }
    if(file != "") {
        dev.off();
    }
}

#This graphs a transformed scatter of the categories
varscan.graph_contamination_scatter = function(x, genome="", file="") {
    if(file != "") {
        pdf(file=file,width=10,height=7.5);
    }
    plot.default(x=x$V11,y=x$V20/x$V11,xlab="Tumor Variant Allele Frequency",ylab="(Tumor Freq - Normal Freq)/Tumor Freq",main=paste(genome, " Contamination Plot",sep=""), type="p", pch=19,cex=0.4, col="black");
    if(file != "") {
        dev.off();
    }
}

#This reads the varscan snp output and normalizes the percentages to numbers in R.
#It also adds a column to the end, the difference in variant allele frequency.
#If header is True then the first line is skipped. This means column names are not propagated into R.
#Min tumor depth and min normal depth can be used to filter to higher depth data only.

varscan.load_snp_output = function(filename="varScan.output.snp",header=T,sep="\t", min_tumor_depth=0, min_normal_depth=0) {
    skip=NULL;
    if(header) {
        skip = 1;
    }
    else {
        skip = 0;
    }
    x <- read.table(file=filename,skip=skip,sep=sep);

    #apply read depth filter
    x = x[ ((x$V5+x$V6) >= min_normal_depth) & ((x$V9+x$V10) >= min_tumor_depth), ];

    #convert percents to numbers
    x$V7 = as.numeric(sub('%','',(x$V7)));
    x$V11 = as.numeric(sub('%','',(x$V11)));
    x = cbind(x,V20=x$V11-x$V7);   #add on allele frequency difference

    invisible(x); #return the loaded file
}
    
    
