package Genome::Model::Tools::CopyNumber::ConvertCnvhmmOutputToSaneFormat;

##############################################################################
#
#
#	AUTHOR:		Chris Miller (cmiller@genome.wustl.edu)
#
#	CREATED:	05/05/2011 by CAM.
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

class Genome::Model::Tools::CopyNumber::ConvertCnvhmmOutputToSaneFormat {
    is => 'Command',
    has => [
	cnvhmm_file => {
	    is => 'String',
	    is_optional => 1,
	    doc => 'CNVHMM output file to convert (also works for CNAHMM)',
	},
	output_file => {
	    is => 'String',
	    is_optional => 1,
	    doc => 'an output file (1-based format, 5 col like CBS)',
	},

    ]
};

sub help_brief {
    "convert the awful cnvhmm format to something usable (5 col: chr, start, stop #markers, cn)"
}

sub help_detail {
    "convert the awful cnvhmm format to something usable (5 col: chr, start, stop #markers, cn)"
}


#########################################################################

#-----------------------------------------------------
#convert cnvhmm output to a format we can use here
sub execute {
    my $self = shift;

    open(OUTFILE,">" . $self->output_file) || die "can't open file\n";
    
    if(defined($self->cnvhmm_file)){
        #read and convert the cnvhmm output
        my $inFh = IO::File->new( $self->cnvhmm_file ) || die "can't open file\n";
        my $inCoords = 0;
        while( my $line = $inFh->getline )
        {
            chomp($line);
            if ($line =~ /^#CHR/){
	    $inCoords = 1;
	    next;
            }
            if ($line =~ /^---/){
                $inCoords = 0;
                next;
            }
            
            if ($inCoords){
                my @fields = split("\t",$line);
                print OUTFILE join("\t",($fields[0],$fields[1],$fields[2],$fields[4],$fields[6])) . "\n";
            }
        }
        close(OUTFILE);
        $inFh->close;                
    } 

    return 1;
}

1;
