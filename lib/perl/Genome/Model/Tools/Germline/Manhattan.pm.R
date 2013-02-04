# This code is based on code from Stephen Turner here: http://gettinggeneticsdone.blogspot.com/2010/01/gwas-manhattan-plots-and-qq-plots-using.html
# http://StephenTurner.us/
# http://GettingGeneticsDone.blogspot.com/
# See license at http://gettinggeneticsdone.blogspot.com/p/copyright.html

#will need to preprocess as so
#d<-read.table(filename,header=T)
#strs=strsplit(as.character(d$x),"_")
#d=transform(d,CHR=sapply(strs,"[[",1),BP=as.numeric(sapply(strs,"[[",2)))
#d$CHR = sub("^X(.+)",replacement="\\1",d$CHR,fixed=FALSE)
#d$CHR = sub("^X",replacement="23",d$CHR,fixed=FALSE)
#d$CHR = as.numeric(sub("^Y",replacement="24",d$CHR,fixed=FALSE))
#names(d)[8]<-"P"

#the above will give a CHR and BP field that is ok. Probably would be easier and more productive to make those fields the output of the actual code...

# manhattan plot using ggplot2
manhattan = function(dataframe, title=NULL, max.y="max", suggestiveline=0, genomewideline=-log10(5e-8), size.x.labels=9, size.y.labels=9, annotate=F, SNPlist=NULL) {
    library(ggplot2)
    if (annotate & is.null(SNPlist)) stop("You requested annotation but provided no SNPlist!")
	d=dataframe
	
	#limit to only chrs 1-23?
	d=d[d$CHR %in% 1:23, ]
	if ("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d) ) {
		d=na.omit(d)
		d=d[d$P>0 & d$P<=1, ]
		d$logp = -log10(d$P)
		d$pos=NA
		ticks=NULL
		lastbase=0
		#new 2010-05-10
		numchroms=length(unique(d$CHR))
		if (numchroms==1) {
			d$pos=d$BP
		} else {
		
			for (i in unique(d$CHR)) {
				if (i==1) {
					d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
				}	else {
					lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
					d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
				}
				ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
			}
			ticklim=c(min(d$pos),max(d$pos))

		}

		mycols=rep(c("gray10","gray60"),max(d$CHR))
		if (max.y=="max") maxy=ceiling(max(d$logp)) else maxy=max.y
		if (maxy<8) maxy=8
		if (annotate) d.annotate=d[as.numeric(substr(d$SNP,3,100)) %in% SNPlist, ]
		if (numchroms==1) {
			plot=qplot(pos,logp,data=d,ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"))
		}	else {
			plot=qplot(pos,logp,data=d, ylab=expression(-log[10](italic(p))) , colour=factor(CHR))
			plot=plot+scale_x_continuous(name="Chromosome", breaks=ticks, labels=(unique(d$CHR)))
			plot=plot+scale_y_continuous(limits=c(0,maxy), breaks=1:maxy, labels=1:maxy)
			plot=plot+scale_colour_manual(value=mycols)
		}
		if (annotate) 	plot=plot + geom_point(data=d.annotate, colour=I("green3")) 
		plot=plot + opts(legend.position = "none") 
		plot=plot + opts(title=title)
		plot=plot+opts(
			panel.background=theme_blank(), 
			panel.grid.minor=theme_blank(),
			axis.text.x=theme_text(size=size.x.labels, colour="grey50",angle=90), 
			axis.text.y=theme_text(size=size.y.labels, colour="grey50"), 
			axis.ticks=theme_segment(colour=NA)
		)
		if (suggestiveline) plot=plot+geom_hline(yintercept=suggestiveline,colour="blue", alpha=I(1/3))
		if (genomewideline) plot=plot+geom_hline(yintercept=genomewideline,colour="red")
		plot
	}	else {
		stop("Make sure your data frame contains columns CHR, BP, and P")
	}
}
