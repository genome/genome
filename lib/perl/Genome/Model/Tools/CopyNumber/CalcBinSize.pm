package Genome::Model::Tools::CopyNumber::CalcBinSize;

use strict;
use warnings;
use Genome;
use Cwd;
use FileHandle;

#this is a wrapper to precess illumina file and make copy number graph, swt test. 

class Genome::Model::Tools::CopyNumber::CalcBinSize {
    is => 'Command',
    has => [
	read_count => {
	    is => 'Integer',
	    is_optional => 0,	
	    doc => 'The number of reads in the sample',
	},
	output_dir => {
	    is => 'String',
	    is_optional => 1,
	    default => getcwd(),
	    doc => 'Directory to use for (small) temporary file',
	},
	entrypoints => {
	    is => 'String',
	    is_optional => 0,
	    doc => 'path to an entrypoints file containing chr, length and ploidy. See (/gscmnt/sata921/info/medseq/cmiller/annotations/entrypoints.hg18.[male|female])',
	},
	mapability => {
	    is => 'String',
	    is_optional => 0,
	    doc => 'path to an entrypoints file containing the mapability of the genome for the read length that you are using. See (~cmiller/annotations/mapability.hg18.[readlength]bpReads.dat)',
	},
	matched_normal => {
	    is => 'Boolean',
	    is_optional => 1,
	    default => 1,
	    doc => ' Flag indicating sample is standalone (not a tumor/normal pair)',
	},
	p_value => {
	    is => 'Float',
	    is_optional => 1,
	    default => 0.01,
	    doc => ' The probability that any given window is misclassified',
	},
	gain_fraction => {
	    is => 'Float',
	    is_optional => 1,
	    default => 0.05,
	    doc => 'The fraction of the genome that is copy number amplified. A good estimation will help the tool choose a better bin size, but isn\'t strictly necessary.',
	},
	loss_fraction => {

	    is => 'Float',
	    is_optional => 1,
	    default => 0.05,
	    doc => 'The fraction -of the genome that is copy number deleted. A good estimation will help the tool choose a better bin size, but isn\'t strictly necessary.',
	},
	overdispersion => {
	    is => 'String',
	    is_optional => 1,
	    default => 3,
	    doc => 'The amount of overdispersion observed in the data. The default is sensible for Illumina reads',
	},
	verbose => {
	    is => 'Boolean',
	    is_optional => 1,
	    default => 0,
	    doc => 'Include extra output, including some plots (written to Rplots.pdf)',
	},
	]
};

sub help_brief {
    "calculate a good window size for estimating copy number from WGS data"
}

sub help_detail {
    "This script takes the number of reads you've got as input, models the reads using a negative binomial distribution, then calculates a bin size. This bin size will enable good separability between the diploid and triploid peaks without misclassifying a greater fraction of windows than specified by the input p-value. 

It will also return info that gives you an idea of how many consecutive altered windows should be required during the segmentation step (CNAseg.pl: option -n) to maximize the resolution while minimizing the chance of obtaining a false positive alteration call.
"
}

sub execute {
    my $self = shift;
    my $read_count = $self->read_count;
    my $entrypoints = $self->entrypoints;
    my $mapability = $self->mapability;
    my $p_value = $self->p_value;
    my $gain_fraction = $self->gain_fraction;
    my $loss_fraction = $self->loss_fraction;
    my $overdispersion = $self->overdispersion;
    my $verbose = $self->verbose;
    my $matched_normal = $self->matched_normal;
    my $output_dir = $self->output_dir;


    #create a temp file for the R command to run
    my ($tfh,$tmpFile) = Genome::Sys->create_temp_file;
    unless($tfh) {
        $self->error_message("Unable to create temporary file $!");
        die;
    }
    my $outFh =  open (my $tmpFileH, ">>$tmpFile") || die "Can't open output file.\n";

    print $tmpFileH '

##-------------------------------------------------
## read in entrypoints file, return a data frame
##
readEntrypoints <- function(afile){
  p=read.table(afile,sep="\t",quote="",colClasses=c("character","numeric","numeric"))
   names(p)=c("chr","length","ploidy")  
  return(p)  
}


##-------------------------------------------------
## calculate the appropriate bin size,
## gain/loss thresholds, number of chromosomes, etc
##
getBinSize <-function(numReads,entrypoints,expectedLoss,expectedGain,fdr,overDispersion,minSize,matchedNormal){

  ## get genome size from entrypoints file, adjust by mappability estimates
  genomeSize = sum(entrypoints$length)
  mapPerc=sum(entrypoints$mapPerc*entrypoints$length)/sum(entrypoints$length)
  effectiveGenomeSize = genomeSize * mapPerc

  if(verbose){
    cat(numReads," total reads\n")
    cat("genome size:",genomeSize,"\n")
    cat("genome mapability percentage:",mapPerc,"\n")  
    cat("effectiveGenomeSize:",effectiveGenomeSize,"\n")
  }


  ## median value has to be adjusted if we have chromosomes with single ploidy
  ## or expected copy number alterations
  ploidyPerc = ploidyPercentages(effectiveGenomeSize, entrypoints, expectedLoss, expectedGain)

  if(verbose){
    cat("expect ",
        ploidyPerc$haploidPerc*100,"% haploid,",
        ploidyPerc$diploidPerc*100,"% diploid,",
        ploidyPerc$triploidPerc*100,"% triploid\n")
  }
    
  ## calculate window size based on triploid peak, since it
  ## produces larger (more conservative) windows
  ploidy = 3
  medAdj = 1
#  if(verbose){
#    if("medAdjustment" %in% names(params)){
#      medAdj = params[["medAdjustment"]]
#    }
#  }
  
  pTrip <- calcWindParams(numReads=numReads,
                          fdr=fdr,
                          genomeSize=effectiveGenomeSize,
                          oDisp=overDispersion,
                          ploidy=ploidy,
                          minSize=minSize,
                          ploidyPerc=ploidyPerc,
                          medAdj=medAdj,
                          matchedNormal=matchedNormal)


  binSize = pTrip$binSize
  binSize = round(binSize/100)*100
  
  return(data.frame(binSize=binSize, hapPerc=ploidyPerc$haploidPerc, dipPerc=ploidyPerc$diploidPerc, tripPerc=ploidyPerc$triploidPerc))
}


##--------------------------------------------------
## calculate adjusted genome size, based on the fact that we
## may have haploid chromosomes and/or expected CN alterations 
ploidyPercentages <- function(effectiveGenomeSize,ents,expectedLoss,expectedGain){
  ##first, get the coverage that come from diploid chrs
  diploidPerc = sum((ents$length*ents$mapPerc)[which(ents$ploidy==2)])/effectiveGenomeSize
  diploidPerc = diploidPerc - expectedLoss 
  diploidPerc = diploidPerc - expectedGain 

  ##coverage from haploid chrs
  haploidPerc = sum((ents$length*ents$mapPerc)[which(ents$ploidy==1)])/effectiveGenomeSize
  haploidPerc = haploidPerc + expectedLoss

  return(data.frame(haploidPerc=haploidPerc,
                    diploidPerc=diploidPerc,
                    triploidPerc=expectedGain))
}




##--------------------------------------------------
## calculate FDR rate and optimal dividing line
##
fdrRate <- function(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj, matchedNormal){
  numWinds <- genomeSize/windSize

  med <- (numReads/(((genomeSize*ploidyPerc$haploidPerc/2) +
                     (genomeSize*ploidyPerc$triploidPerc*3/2) +
                     (genomeSize*ploidyPerc$diploidPerc)) / windSize))*medAdj
  
  medAlt <- med*(ploidy/2)
  


  mislabeledWinds <- function(thresh, amed, bmed, oDisp,
                              aNum, bNum, matchedNormal){

    if(matchedNormal){
      lowreads = (rnbinom(n=aNum,size=(amed/oDisp-1),mu=amed)/
                  rnbinom(n=aNum,size=(amed/oDisp-1),mu=amed))
      hireads = (rnbinom(n=bNum,size=(bmed/oDisp-1),mu=bmed)/
                 rnbinom(n=bNum,size=(amed/oDisp-1),mu=amed))
    } else {
      lowreads = (rnbinom(n=aNum,size=(amed/oDisp-1),mu=amed)/med)
      hireads = (rnbinom(n=bNum,size=(bmed/oDisp-1),mu=bmed)/med)
    }

    ## probability of getting a miscall is simulated
    ## number of normal bins that are called amplified
    low=length(which(lowreads > thresh/2))
    ## number of amp bins called normal
    high=length(which(hireads < thresh/2))


    
    ##plotting
    if(verbose){
      hist(2*lowreads,xlim=c(0,5),breaks=seq(-5,20,0.05),col=rgb(0, 1, 0, 0.5),xlab="tumor/normal ratio",main="binned read distribution")
      hist(2*hireads,xlim=c(0,5),breaks=seq(-5,20,0.05),col=rgb(0, 0, 1, 0.5),add=T)
      abline(v=thresh)
    }
    return(low+high)
  }

  
  fdr <- 0
  if(ploidy < 2){

    aNum=numWinds*ploidyPerc$haploidPerc
    bNum=numWinds*ploidyPerc$diploidPerc
    fdr <- mislabeledWinds(1.5,medAlt, med, oDisp, aNum, bNum, matchedNormal)/(aNum+bNum)
  } else {
    aNum=numWinds*ploidyPerc$diploidPerc
    bNum=numWinds*ploidyPerc$triploidPerc
    fdr <- mislabeledWinds(2.5, med, medAlt, oDisp, aNum, bNum, matchedNormal)/(aNum+bNum)
  }
  return(fdr)
                                        #  return(data.frame(fdr=divFdr$fdr, div=divFdr$div, med=med))
}




##--------------------------------------------------
## find the optimal dividing line between the 
## two specified peaks. Assumes a poisson 
## distribution with overdispersion
##
dividePeaks <- function(amed,bmed,aNum,bNum,oDisp){  
  thresholds =(amed+1):(bmed-1)

  mislabeledWinds <- function(thresh){
    low=(1-pnbinom(thresh,size=(amed/oDisp-1),mu=amed))*aNum
    high=(pnbinom(thresh,size=(bmed/oDisp-1),mu=bmed))*bNum
    return(low+high)
  }

  vals=sapply(thresholds,mislabeledWinds)
  lows = which(vals == min(vals))

  ##choose the low point closest
  ##to the halfway point
  diffFromMed = abs((lows+amed)-(amed+bmed/2))
  pos = which(diffFromMed == min(diffFromMed))
  if(length(pos) > 1){
    pos = pos[1]
  }
  return(data.frame(div=lows[pos]+amed,fdr=mislabeledWinds(lows[pos]+amed)/(aNum+bNum)))
}



#-------------------------------------------------
# Calculates the window size that conforms to the
# given FDR rate and the threshold that best
# separates the peaks of ploidy
#
calcWindParams <- function(numReads,fdr,genomeSize,oDisp, ploidy, minSize, ploidyPerc, medAdj, matchedNormal, startDiv=100){#nullThis, nullBoth, startDiv=100){
  ## answer has to be this close to the FDR rate (on the lower side)
  tolerance = 0.005
  #starting point for search
  windSize = genomeSize/startDiv
  divider = NULL
  
  ## first, halve size until we get above FDR threshold
  found <- FALSE
  pfdr <- fdrRate(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj,matchedNormal)

  if(pfdr > fdr){
    stop("not enough reads to achieve the desired FDR rate")
  }
  
  while((found == FALSE) & (windSize > minSize)){
    windSize <- round(windSize / 2)
    pfdr <- fdrRate(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj,matchedNormal)
    
    if(pfdr > fdr){
      found = TRUE
    }
  }

  if(windSize > minSize){
    ## zero in on a size that is within the desired parameters
    found = FALSE    
    adj = windSize/2

    while((found == FALSE) & (windSize > minSize)){    
      if(pfdr > fdr){
        windSize <- round(windSize + adj)
        pfdr <- fdrRate(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj,matchedNormal)
        
      }else if(pfdr < (fdr - tolerance)){
        windSize <- round(windSize - adj)
        pfdr <- fdrRate(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj,matchedNormal)
        
      } else{
        found = TRUE
      }
      adj <- adj/2
    }
  }
    
  ##if the window size is below the minimum size, have to
  ##recalculate params
  if(windSize < minSize){
    windSize = minSize
  } else {
    ## round to multiple of minSize
    windSize = floor(windSize/minSize)*minSize
  }
  
  pfdr <- fdrRate(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj,matchedNormal)
  div=((ploidy/2)+1)
  return(data.frame(binSize=windSize,div=div))
}



addMapability <- function(mapFile,entrypoints){
  entrypoints = cbind(entrypoints,mapPerc=rep(1,length(entrypoints$chr)))
  tmp = scan(mapFile,what="",quiet=TRUE)
  for(i in seq(1,length(tmp),2)){
    mapp=as.numeric(tmp[i+1])
    if(tmp[i] %in% entrypoints$chr){
      entrypoints[which(entrypoints$chr==tmp[i]),]$mapPerc = mapp
    }
  }
  return(entrypoints)
}

############################################################

# perl script inserts code here


';



    print $tmpFileH "cat(\"input: $read_count reads\\n\")\n";
    print $tmpFileH "cat(\"matched normal: $matched_normal\\n\\n\")\n";
    print $tmpFileH "verbose <<- $verbose\n";
    print $tmpFileH "numReads = $read_count\n";
    print $tmpFileH "entrypointsFile = \"$entrypoints\"\n";
    print $tmpFileH "entrypoints = readEntrypoints(entrypointsFile)\n";
    print $tmpFileH 'entrypoints <- addMapability("' . $mapability . '",entrypoints)' ."\n";
    print $tmpFileH "expectedLoss = $loss_fraction\n";
    print $tmpFileH "expectedGain = $gain_fraction\n";
    print $tmpFileH "fdr = $p_value\n";
    print $tmpFileH "overDispersion = $overdispersion\n";
    print $tmpFileH "minSize = 100\n";
    print $tmpFileH "matchedNormal = $matched_normal\n";
    print $tmpFileH "genomeSize = sum(entrypoints\$length)\n";
    print $tmpFileH "binInfo=getBinSize(numReads,entrypoints,expectedLoss,expectedGain,fdr,overDispersion,minSize,matchedNormal)\n";
    print $tmpFileH "cat(\"binSize: \",binInfo\$binSize,\"\\n\")\n";
    print $tmpFileH "cat(\"expected false calls with N consecutive windows:\\n\")\n";
    print $tmpFileH "for(i in seq(2,5)){\n";
    print $tmpFileH "cat(i,\": \",genomeSize/binInfo\$binSize*(fdr^i),\"\\n\")\n";
    print $tmpFileH "}\n";
    $tmpFileH->close;

    system("Rscript --vanilla $tmpFile");

    return 1;
}
1;
