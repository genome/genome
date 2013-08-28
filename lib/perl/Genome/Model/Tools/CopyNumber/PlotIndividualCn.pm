package Genome::Model::Tools::CopyNumber::PlotIndividualCn;

use strict;
use Genome;
use IO::File;
use File::Basename;
use warnings;
require Genome::Sys;
use FileHandle;
use File::Spec;


class Genome::Model::Tools::CopyNumber::PlotIndividualCn{
    is => 'Command',
    has => [

        segment_file => {
	    is => 'String',
	    is_optional => 0,
	    doc => '5-column segmented CNA file - chr, start, stop, numWinds, copyNumber (usually named "alts.paired.dat")',
	},
	tumor_window_file => {
	    is => 'String',
	    is_optional => 1,
	    doc => 'tumor window file to get reads from (output of gmt copy-number bam-window)',
	},
	normal_window_file => {
	    is => 'String',
	    is_optional => 1,
	    doc => 'normal window file to get reads from (output of gmt copy-number bam-window)',
	},
        corrected_window_file => {
            is => 'String',
            is_optional => 1,
            doc =>'the corrected windows output by copyCat with the dumpBins option (usually named "rd.bins.dat")',
        },
        output_pdf => {
            is => 'String',
            is_optional => 0,
            doc =>'file to output plots into',
        },
        annotation_directory => {
            is => 'String',
            is_optional => 0,
            example_values => ['/gscmnt/gc6122/info/medseq/annotations/copyCat/'],
            doc =>'path to the cn annotation directory',
        },
        genome_build => {
            is => 'String',
            is_optional => 0,
            example_values => ['hg18', 'hg19'],
            doc =>'genome build that the data uses',
        },
        gaps_file => {
            is => 'String',
            is_optional => 1,
            example_values => ['/gscmnt/gc6122/info/medseq/annotations/gaps.hg18.ucsc.bed'],
            doc =>'path to a file containing genomic gaps (first three columns must be chr, start, stop)',
        },
        gain_threshold => {
            is => 'Number',
            is_optional => 1,
            default => 2.0,
            doc =>'copy number threshold for gains (default plot everything)',
        },
        loss_threshold => {
            is => 'Number',
            is_optional => 1,
            default => 2.0,
            doc =>'copy number threshold for losses (default plot everything)',
        },
        r_file => {
            is => 'String',
            is_optional => 1,
            doc =>'if provided, will output a file containing the R commands that are run',
        },
        plot_full_chrs => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc =>'output only a view of the entirety of each chromosome',
        }

        ]
};

sub help_brief {
    "Takes inputs and ouputs from CopyCat somatic cn caller and produces plots of each CN suitable for display or manual review"
}

sub help_detail {
    "Takes inputs and ouputs from CopyCat somatic cn caller and produces plots of each CN suitable for display or manual review"
}

#########################################################################

sub execute {
    my $self = shift;
    my $segment_file = $self->segment_file;

    my $rfile;
    my $newfile;
    if(defined($self->r_file)){
        $newfile = $self->r_file;
        open($rfile, ">$newfile") || die "can't open file\n";
    } else {
        ($rfile,$newfile) = Genome::Sys->create_temp_file;
        unless($rfile) {
            $self->error_message("Unable to create temporary file $!");
            die;
        }
    }

    my $dir_name = dirname(__FILE__);
    print $rfile 'source("' . $dir_name . "/PlotIndividualCn.R" . '")' . "\n";

    #read in common files
    if(defined($self->gaps_file)){
        print $rfile 'gaps=read.table("' . $self->gaps_file . '")' . "\n";
    }
    print $rfile 'tumorNormalRatio=NULL' . "\n";

    print $rfile 'db = sqldf() #open temp db' . "\n";

    my $count = 0;
    my $skipped = 0;

    #if plot full chrs, create a new segment file for each chromosome.
    if($self->plot_full_chrs){
        my ($seghandle,$segfile) = Genome::Sys->create_temp_file;
        unless($seghandle) {
            $self->error_message("Unable to create temporary file $!");
            die;
        }        

        my $entryfile = $self->annotation_directory . "/" . $self->genome_build . '/entrypoints.male';
        my $infile = IO::File->new( $entryfile ) || die "can't open segment file\n";
        while( my $line = $infile->getline )
        {
            chomp($line);
            my @F = split("\t",$line);
            print $seghandle join("\t",($F[0],"1",$F[1],"0","2")) . "\n";
        }
        close($infile);
        close($seghandle);
        $segment_file = $segfile;        
    }



    my $inFh = IO::File->new( $segment_file ) || die "can't open segment file\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        my @F = split("\t",$line);

        unless ($F[4] >= $self->gain_threshold || $F[4] <= $self->loss_threshold){
            $skipped++;
            next;
        }

        print $rfile 'chr = "' . $F[0] . "\"\n";
        print $rfile 'xmin = ' . $F[1] . "\n";
        print $rfile 'xmax = ' . $F[2] . "\n";
        print $rfile 'print("plotting alt: ' . "$F[0]:$F[1]-$F[2]" . '")' . "\n";

        #add the padding
        unless($self->plot_full_chrs){

            if(($F[2]-$F[1]) > 50000){
                print $rfile 'zxmin = xmin-(xmax-xmin)' . "\n";
                print $rfile 'xmax = xmax+(xmax-xmin)' . "\n";
                print $rfile 'xmin=zxmin' . "\n";
            } else {
                print $rfile 'xmin = xmin-50000' . "\n";
                print $rfile 'xmax = xmax+50000' . "\n";
            }
        }

        if(defined($self->tumor_window_file)){
            print $rfile 'tumwinds = getWindows(chr, xmin, xmax, header=T';
            if($count == 0){
                print $rfile ', windowFile="' . $self->tumor_window_file . '"';
            }
            print $rfile ', tableName="tumor")' . "\n";
        }
        if(defined($self->normal_window_file)){
            print $rfile 'nrmwinds = getWindows(chr, xmin, xmax, header=T';
            if($count == 0){
                print $rfile ', windowFile="' . $self->normal_window_file . '"';
            }
            print $rfile ', tableName="normal")' . "\n";
        }
        if(defined($self->corrected_window_file)){
            print $rfile 'cwinds = getWindows(chr, xmin, xmax, header=T';
            if($count == 0){
                print $rfile ', windowFile="' . $self->corrected_window_file . '"';
            }
            print $rfile ', tableName="corrected")' . "\n";
        }

        if(defined($self->tumor_window_file)){
            print $rfile "if(is.null(tumorNormalRatio)){\n";
            print $rfile "  tumorNormalRatio=getTumorNormalRatio(\"tumor\",\"normal\")\n";
            print $rfile "}\n";
        }

        #use original segment file here so they all get drawn
        print $rfile 'plotSegments(filename="' . $self->segment_file;
        print $rfile '", entrypoints="' . $self->annotation_directory . "/" . $self->genome_build . '/entrypoints.male"';
        print $rfile ", ymin=0, ymax=8, xmin=xmin, xmax=xmax, chr=chr";
        if(defined($self->tumor_window_file)){
            print $rfile ", tumorWindows=tumwinds";
        }
        if(defined($self->normal_window_file)){
            print $rfile ", normalWindows=nrmwinds";
        }
        if(defined($self->corrected_window_file)){
            print $rfile ", correctedWindows=cwinds";
        }
        if(defined($self->gaps_file)){
            print $rfile ", gaps=gaps";
        }

        print $rfile ", lossThresh=" . $self->loss_threshold;
        print $rfile ", gainThresh=" . $self->gain_threshold;


        if(defined($self->normal_window_file) || defined($self->tumor_window_file)){
            print $rfile ", coverageTracks=TRUE"
        }

        print $rfile ', plotTitle="' . "$F[0]:$F[1]-$F[2]" . '"';
        print $rfile ', highlightedSegs=c("' . "$F[0]:$F[1]-$F[2]" . '")';

        if($count == 0){
            print $rfile ', pdfOpen="TRUE"';
        }
        print $rfile ', pdfFile="' . $self->output_pdf . '"' . ")\n";
        $count++;
    }
    close($inFh);

    print $rfile "dev=dev.off()\n";
    print $rfile 'sqldf(); #close db' . "\n";
    close($rfile);

    if($skipped){
        print STDERR "skipped $skipped sites that did not exceed the gain/loss thresholds\n";
    }

    if($count == 0){
        print STDERR "no sites to plot, exiting\n";
        return 1;
    }

    my $cmd = "Rscript $newfile";
    my $return = Genome::Sys->shellcmd(
        cmd => "$cmd",
        );
    unless($return) {
        $self->error_message("Failed to execute: Returned $return");
        die $self->error_message;
    }
    return $return;
}
