package Genome::Model::Tools::CopyNumber::SimulateReads;

##############################################################################
#
#
#	AUTHOR:		Chris Miller (cmiller@genome.wustl.edu)
#
#	CREATED:	05/05/2011 by CAM.
#	MODIFIED:
#
#	NOTES:
#
##############################################################################

use strict;
use Genome;
use IO::File;
use Statistics::R;
use File::Basename;
use warnings;
require Genome::Sys;
use FileHandle;

class Genome::Model::Tools::CopyNumber::SimulateReads {
    is => 'Command',
    has => [


        chr => {
            is => 'String',
            is_optional => 0,
            doc => 'the chromosome to simulate',
        },

        nreads => {
            is => 'Integer',
            is_optional => 1,
            doc => 'the number of reads to generate',

        },

        coverage => {
            is => 'Float',
            is_optional => 1,
            doc => 'generate enough reads to cover the chromosome this deeply',

        },

        readlength => {
            is => 'Integer',
            is_optional => 0,
            doc => 'the read length',
        },

        overdispersion => {
            is => 'Float',
            is_optional => 1,
            default => 3,
            doc => 'the amount of overdispersion to simulate. Overdispersion is defined as the variance-to-mean ratio of the negative binomal distribution that describes the genomic windows ',
        },

        window_size => {
            is => 'Integer',
            is_optional => 1,
            doc => 'the window size to use for generating the distribution. The default should be fine for almost all applications',
            default => 10000,
        },


        output_file => {
            is => 'String',
            is_optional => 0,
            doc => 'the sam file that will contain the simulated reads',
        },

        entrypoints_file => {
            is => 'String',
            is_optional => 1,
            default => "/gscmnt/sata921/info/medseq/cmiller/annotations/entrypoints.hg18.male",
            doc => 'entrypoints for the genome build to simulate (also determines gender if doing WGS simulation).',
        },

        alt_size => {
            is => 'Integer',
            is_optional => 1,
            doc => 'Size in bp of alteration to add',
        },

        alt_cn => {
            is => 'Integer',
            is_optional => 1,
            doc => 'the copy number of the alteration to add',
            default => 3,
        },

        alt_file =>{
            is => 'String',
            is_optional => 1,
            doc => 'output the position of the alteration for validating the output',
        },

        ]
};

sub help_brief {
    "generate a sam file full of simulated reads for benchmarking CN algorithms"
}

sub help_detail {
    "generate a sam file full of simulated reads for benchmarking CN algorithms"
}


################################################################################

sub printSamLines{
    my ($index, $chr, $readlength, $outfile, @reads) = @_;
        foreach my $pos (@reads){
            #half on pos strand, half on neg
            my $flag = 0;
            my $strand = "+";
            if (int(rand(2))){
                $flag = 16;
                $strand = "-";
            }
            print $outfile join("\t", $index, $flag, $chr, $pos, "255", $readlength ."M", $strand, 0, 0, "*", "*") . "\n";
            $index++;
    }
    return($index);
}

#------------------------------------

sub getChrLength{
    my ($entrypoints_file, $chr) = @_;
    my $inFh = IO::File->new($entrypoints_file) || die "can't open file\n";
    
    while( my $line = $inFh->getline )
    {
        chomp($line);
        my @fields = split("\t",$line);
        if ($fields[0] == $chr){
            return($fields[1]);
        }        
    }
    die("chr $chr doesn't exist in entrypoints file")
}



#------------------------------------
sub getDist{
    my ($regionSize, $nreads, $overdispersion, $window_size, $self) = @_;
    #one temp file for the output dist
    my ($tfh,$tmpfile) = Genome::Sys->create_temp_file;
    unless($tfh) {
        $self->error_message("Unable to create temporary file $!");
        die;
    }


    # and a temp file for the r commands
    my ($rfh,$rfile) = Genome::Sys->create_temp_file;
    unless($rfh) {
        $self->error_message("Unable to create temporary file $!");
        die;
    }


    open(OUTFILE,">$rfile") || die "can't open temp file for writing ($rfile)\n";


    print OUTFILE '

##--------------------------------------------------
## generate a distribution by modelling the overdispersed
## poisson with the negative binomial distribution
## d = var/mean
##
rpois.od<-function (n, lambda,d=1) {
  if (d==1)
    rpois(n, lambda)
  else
    rnbinom(n, size=(lambda/(d-1)), mu=lambda)
}


createDist <- function(regsize, nreads, overdispersion, wind.size, output.file){

  num.winds = ceiling(regsize/wind.size)
  median = nreads/num.winds
  ##we know that when binned, the reads follow a poisson-family
  ##distribution (possible with overdispersion), so generate
  ##such a distribution so we know how many reads should be in
  ##each region
  adist = rpois.od(num.winds, median, overdispersion)

  ## due to rounding, we sometimes end up with a smidge more
  ## or less reads than we specified. Add or remove reads 
  ## randomly to get the right number
  if(sum(adist) > nreads){
    sp=sum(adist)-nreads
    for(i in 1:sp){
      bin=round(runif(1,1,length(adist)));
      adist[bin]=adist[bin]-1;
    }
  } else if(nreads > sum(adist)){
    sp=nreads-sum(adist)
    for(i in 1:sp){
      bin=round(runif(1,1,length(adist)));
      adist[bin]=adist[bin]+1;
    }
  }

  write.table(adist,file=output.file, sep="\t",
                col.names=FALSE, quote=FALSE,
                row.names=FALSE)

}

###########################
';

    print OUTFILE "createDist(";
    print OUTFILE "regsize=" . $regionSize . ",";
    print OUTFILE "nreads=" . $nreads . ",";
    print OUTFILE "overdispersion=" . $overdispersion . ",";
    print OUTFILE "wind.size=" . $window_size . ",";
    print OUTFILE "output.file=\"" . $tmpfile . "\"";
    print OUTFILE ")\n";

    print OUTFILE "q()\n";
    close OUTFILE;

#    `cp $rfile /tmp/asdf.R`;

    #now run the R command
    my $cmd = "R --vanilla --slave \< $rfile";
    my $returnval = Genome::Sys->shellcmd(
        cmd => "$cmd",
        );
    unless($returnval) {
        $self->error_message("Failed to execute: Returned $returnval");
        die $self->error_message;
    }

    #now, read in the distribution from the R file
    my @dist;
    my $inFh = IO::File->new( $tmpfile ) || die "can't open file\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        push(@dist,$line);
    }
    return @dist;
}

#--------------------------------
#generate the specified number of reads in a given interval
sub genReads{
    my ($st, $sp, $num) = @_;    
    my $count = 0;
    my @positions;
    
    #randomly place reads within each window
    while($count < $num){
        push(@positions, (int(rand($sp-$st))+$st));
        $count++;
    }
    @positions = sort {$a <=> $b} @positions;
    return(@positions);
}


############################################################################

sub execute {
    my $self = shift;
    my $chr = $self->chr;
    my $nreads = $self->nreads;
    my $coverage = $self->coverage;
    my $readlength = $self->readlength;
    my $overdispersion = $self->overdispersion;
    my $window_size = $self->window_size;
    my $output_file = $self->output_file;
    my $entrypoints_file = $self->entrypoints_file;
    my $alt_size = $self->alt_size;
    my $alt_cn = $self->alt_cn;
    my $alt_file = $self->alt_file;



    if(!(defined($coverage)) && !(defined($nreads))){
        die "either coverage or nreads must be specified";
    }
   
    my $chrlength = getChrLength($entrypoints_file, $chr);
    if(defined($coverage)){
        $nreads = ($chrlength/$readlength)*$coverage;
    }
    #get normal windows for the whole chr
    my @chrwinds = getDist($chrlength, $nreads, $overdispersion, $window_size, $self);
     
    my @altwinds;
    my $alt_pos;
    my $alt_stop;
    if(defined($alt_size)){
        #get alt windows for the alt region
        my $altreads = int($nreads * ($alt_size/$chrlength));
        @altwinds = getDist($alt_size, $altreads*($alt_cn/2), $overdispersion, $window_size, $self);
        
        #randomly choose a position for the alt    
        $alt_pos = int(rand($chrlength-$alt_size));
        $alt_stop = $alt_pos + $alt_size - 1;

        print STDERR "alt at: $chr:$alt_pos-$alt_stop\n";
    }        

    #step through each window, generate the right number of reads 
    #output them in sam format

    open(my $outfile,">$output_file") || die "can't open output file for writing)\n";
    #step through each window
    my $index = 0;        
    my $altindex = 0;
    my $readindex = 0;

    while ($index < @chrwinds){
        my $st = $index*$window_size;
        my $sp = $index*$window_size + ($window_size-1);
        #if no alteration, just output it
        if(!(defined($alt_size))){
            my @reads = genReads($st, $sp, $chrwinds[$index]);
            $readindex = printSamLines($readindex, $chr, $readlength, $outfile, @reads);
        #else have to deal with alterations
        } else { 
            
            #if this window is completely covered by the alt
            if (($st >= $alt_pos) && ($sp < $alt_stop)){
                #replace the window with an alt window
                my @reads = genReads($st, $sp, $altwinds[$altindex]);
                $readindex = printSamLines($readindex, $chr, $readlength, $outfile, @reads);
                $chrwinds[$index] = $altwinds[$altindex];
                $altindex++;
            
            #else if it is the front window, partially covered
            } elsif (($st < $alt_pos) && ($sp > $alt_pos)){
                #get the number of bases in the window that are altered
                my $altbases = ($sp-$alt_pos);
                my $altperc = $altbases/$window_size;                

                #generate normal reads                
                my $num_normreads = sprintf("%.0f", ($chrwinds[$index]*(1-$altperc)));
                my @norm_reads = genReads($st, $alt_pos, $num_normreads);
                #generate altered reads
                my $num_altreads = sprintf("%.0f",$altwinds[@altwinds-1]*$altperc);
                my @alt_reads = genReads($alt_pos+1,$sp, $num_altreads);

                #print them all
                $readindex = printSamLines($readindex, $chr, $readlength, $outfile, (sort{$a <=> $b}(@norm_reads,@alt_reads)));

            #else if it is the rear window, partially covered
            } elsif (($st < $alt_stop) && ($sp > $alt_stop)){

                #get the number of bases in the window that are nonaltered
                my $normbases = ($sp-$alt_stop);
                my $normperc = $normbases/$window_size;                

                #generate altered reads
                my $num_altreads = sprintf("%.0f",$altwinds[@altwinds-1]*(1-$normperc));
                my @alt_reads = genReads($st, $alt_stop, $num_altreads);

                #generate normal reads
                my $num_normreads = sprintf("%.0f", ($chrwinds[$index]*($normperc)));
                my @norm_reads = genReads($alt_stop, $sp, $num_normreads);

                #print them all
                $readindex = printSamLines($readindex, $chr, $readlength, $outfile, (sort {$a <=> $b}(@alt_reads,@norm_reads)));
                
            #else we're in a normal window, just output the normal number of reads
            } else {
                my @reads = genReads($st, $sp, $chrwinds[$index]);
                $readindex = printSamLines($readindex, $chr, $readlength, $outfile, @reads);
            }

        }
        $index++
    }



    ##if specified, print out the alteration information
    #chr, st, sp, ploidy
    if(defined($alt_file)){
        open(ALTFILE,">$alt_file") || die "can't open file for writing ($alt_file)\n";
        print ALTFILE join("\t",($chr,$alt_pos,$alt_stop,$alt_cn)) . "\n";
        close(ALTFILE)
    }



    ##print out some info
    print STDERR "\n";
    print STDERR "  --- SAM file created at $output_file ---\n";
    print STDERR "\n";
    print STDERR "  to convert to BAM, run:\n";
    print STDERR "  samtools view -bt ~/sata921/NCBI-human-build36/$chr.fa.fai $output_file > " . $output_file . ".bam\n";
    print STDERR "\n";
    print STDERR "  to convert to BED, run:\n";
    print STDERR '  awk \'{print $2"\t"$3"\t"$3+' . $readlength . "}' < $output_file >" . $output_file . ".bed\n";
    print STDERR "\n";
    print STDERR "\n";
    close($outfile);   

    
    return 1;
}

