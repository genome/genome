# Fetch command line arguments
args = commandArgs();
input_matrix = as.character(args[4]);
genes_to_plot = as.character(args[5]);
output_pdf = as.character(args[6]);
preserveGeneOrder = as.numeric(as.character(args[7]));

sort.data.frame <- function( x, by ) {
    if(by[[1]] != "~")
        stop("Argument 'by' must be a one-sided formula.")

    ## Make the formula into character and remove spaces
    formc <- as.character(by[2])
    formc <- gsub(" ", "", formc)
    ## If the first character is not + or -, add +
    if(!is.element(substring(formc, 1, 1), c("+", "-")))
        formc <- paste("+", formc, sep = "")

    ## Extract the variables from the formula
    vars <- unlist(strsplit(formc, "[\\+\\-]"))
    vars <- vars[vars != ""] # Remove any extra "" terms

    ## Build a list of arguments to pass to "order" function
    calllist <- list()
    pos <- 1 # Position of + or -
    for(i in 1:length(vars)){
        varsign <- substring(formc, pos, pos)
        pos <- pos + 1 + nchar(vars[i])
        if(is.factor(x[, vars[i]])){
        if(varsign == "-") {
            calllist[[i]] <- -rank(x[, vars[i]])
        } else {
            calllist[[i]] <- rank(x[, vars[i]])
        }
        } else {
        if(varsign == "-") {
            calllist[[i]] <- -x[, vars[i]]
        } else {
            calllist[[i]] <- x[,vars[i]]
        }
        }
    }
    return( x[do.call( "order", calllist ), ] )
}

plotMutRel <- function( infile, genes, outfile, preserveGeneOrder=FALSE) {

    ##------------------
    a = read.table(infile,row.names=1,header=T)
    gene_list = unlist(strsplit(genes, split=","))
    print(genes)
    df=a[,gene_list]

    numSamp=length(df[,1])
    numGenes=length(df)

    if(numGenes < 1){
        return("Error: genes to plot not found in matrix")
    }

    ##adjustments to plot and text sizes for different numbers of samples
    samptext=0.3
    genetext=0.75
    pdfwidth=numSamp/10
    offset=2

    if(numSamp < 50){
        samptext=0.3
        #genetext=0.85
        pdfwidth=numSamp/5
        offset=1.25
    }
    if(numSamp < 35){
        samptext=0.4
    }
    if(pdfwidth < 3) { pdfwidth = 3 }
#    if(numSamp >=100) { genetext=1.1 }
#    if(numSamp >=200) { genetext=1.1 }

    pdfheight=3+(0.25*numGenes-1)
    ##-------

    

    #sort the data using the number of mutations in each, desc
    sortdf = data.frame(g=gene_list[1],s=sum(df[,gene_list[1]]))
    for(i in 2:length(gene_list)){
        sortdf = rbind(sortdf, data.frame(g=gene_list[i],s=sum(df[,gene_list[i]])))
    }
    if(!(preserveGeneOrder)){
      gene_list = as.vector(sort.data.frame(sortdf,~-s)$g)
    }

    print(preserveGeneOrder)
    print(gene_list)

    for(i in rev(gene_list)){
        df = sort.data.frame(df,c("~",paste("-",i,sep="")))
    }

    #output pdf here
    pdf(outfile,height=pdfheight,width=pdfwidth)
    par(xpd=T)
    plot(-100,-100, xlim=c(0,numSamp), ylim=c(-numGenes,0), axes=F, xlab="", ylab="")

    hspace=0.10
    vspace=0.05

    ##plot grey rects
    for(i in -(1:numGenes)){
        rect( (1:numSamp)+hspace,
            rep(i,numSamp)+(1-vspace),
            (1:numSamp)+(1-hspace),
            rep(i,numSamp)+vspace,
            col="grey80",
            border=F );
    }

    ##plot color rects
    for(i in -(1:numGenes)){
        pos=which(df[,gene_list[-i]]==1)
        rect(pos+hspace,
            rep(i,length(pos))+(1-vspace),
            pos+(1-hspace),
            rep(i,length(pos))+vspace,
            col="darkgreen",
            border=F )
    }

    ##gene labels
    for(i in -(1:numGenes)){
        text(1,i+0.5,gene_list[-i],pos=2, cex=genetext)
    }

    ##sample labels
    for(i in 1:numSamp){
        text(i+offset,-numGenes,row.names(df)[i],srt=90,pos=2,cex=samptext)
    }

    dev.off()
}

plotMutRel( input_matrix, genes_to_plot, output_pdf, preserveGeneOrder );
