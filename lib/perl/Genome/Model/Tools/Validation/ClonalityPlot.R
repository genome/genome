drawPlot <- function(z1, cn, xchr, additional_plot_points=0, cncircle=0, num_clusters=0, output_filename=NULL){

    # define plot colors
    ptcolor = NULL;
    circlecolor = NULL;
    if (cncircle == 0 || cncircle == 2) { ptcolor = "#67B32E99"; circlecolor = "#67B32E"; }
    if (cncircle == 1) { ptcolor = "#1C366099"; circlecolor = "#1C3660"; }
    if (cncircle == 3) { ptcolor = "#F4981999"; circlecolor = "#F49819"; }
    if (cncircle == 4) { ptcolor = "#E5242099"; circlecolor = "#E52420"; }

    # default plot space
    plot.default(x=(z1$V11), y=(z1$V9+z1$V10), log="y", type="p", pch=19, cex=0.4, col="#00000000", xlim=c(-1,101), ylim=c(5,absmaxx), axes=FALSE, ann=FALSE, xaxs="i", yaxs="i");
    # plot data points (autosomes)
    points(y=(cn$V9+cn$V10), x=(cn$V11), type="p", pch=16, cex=0.75, col=ptcolor);
    # plot data points (X chr)
    points(y=(xchr$V9+xchr$V10),x=(xchr$V11),type="p",pch=2,cex=0.8,col=ptcolor);

    # add in highlight of points selected for by script input
    if(length(additional_plot_points) >1) {
        points(x=additional_plot_points$V2,y=additional_plot_points$V3,type="p",pch=7,cex=0.8,col="#555555FF");
    }

    # define the axis
    axis(side=2,las=1,tck=0,lwd=0,cex.axis=0.6,hadj=0.5);
    for (i in 2:length(axTicks(2)-1)) {
        lines(c(-1,101),c(axTicks(2)[i],axTicks(2)[i]),col="#00000022");
    }

    # plot the background color
    rect(-1, 5, 101, axTicks(2)[length(axTicks(2))]*1.05, col = "#00000011",border=NA);

    ################# Plot clustered points (optional) #######################

    if (num_clusters > 0) {

        # function to plot clusters
        plot_clusters <- function(data,clusters) {
            point_colors <- rainbow(max(clusters),alpha=0.2);
            text_colors <- rainbow(max(clusters),alpha=1.0);

            for(k in 1:max(clusters)) {
                vector_size = rep(7.1,length(subset(data,clusters==k)));
                if (length(subset(data,clusters==k)) > 0) { #sometimes it doesn't assign any points to a cluster...don't understand why
                    points(subset(data,clusters==k),vector_size,col=point_colors[k],cex=0.75,type="p",pch=16);
                    #text(x=mean(subset(data,clusters==k)),y=10,col=text_colors[k],labels=paste("C",k),cex=0.4);
                    text(x=median(subset(data,clusters==k)),y=14.1,col=text_colors[k],labels=round(mean(subset(data,clusters==k)),1),cex=0.7);
                    abline(v = (min(subset(data,clusters==k))-.4),col=text_colors[k],cex=0.5,lwd=0.4) # add vertical line at min of x
                    abline(v = (max(subset(data,clusters==k))+.4),col=text_colors[k],cex=0.5,lwd=0.4) # add vertical line at max of x
                    print(min(subset(data,clusters==k)));
                    print(max(subset(data,clusters==k)));
                }
            } 
        }

        # cluster tumor VAF data using mixtools
        tumor_VAF_data = cn$V11;
        library(mixtools);
        mix_results <- normalmixEM(tumor_VAF_data, k = num_clusters, maxit = 10000, maxrestarts=20);
        posteriors <- mix_results$posterior;
        clusters = NULL;
        for (n in 1:(dim(posteriors)[1])) { clusters[n]=as.numeric(which(posteriors[n,]==max(posteriors[n,]),arr.ind=T)); }
        plot_clusters(tumor_VAF_data,clusters);

        # print output file displaying data points clustered and their associated cluster
        if (!is.null(output_filename)) {
            output = cbind(cn,clusters);
            write.table(output, file=output_filename, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE);
        }

    }

    ########### End Plot clustered points ###################################################

    # add cn circle
    if(cncircle != 0){
        points(x=c(97),y=c(absmaxx*2/5),type="p",pch=19,cex=3,col=circlecolor);
        text(c(97),y=c(absmaxx*2/5), labels=c(cncircle), cex=1, col="#FFFFFFFF")
    }

    # y axis label
    mtext("Tumor Coverage",side=2,cex=0.5,padj=-4.2);
}
