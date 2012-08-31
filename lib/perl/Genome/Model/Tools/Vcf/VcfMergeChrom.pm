package Genome::Model::Tools::Vcf::VcfMergeChrom;

use strict;
use warnings;
use Genome;
use File::stat;
use IO::File;
use File::Basename;
use Getopt::Long;
use FileHandle;

class Genome::Model::Tools::Vcf::VcfMergeChrom {
    is => 'Command',
    has => [
        output_file => {
            is => 'Text',
            is_output => 1,
            is_optional => 0,
            doc => "Output merged VCF",
        },
        vcf_files => {
            is => 'Text',
            is_optional => 0,
            doc => "comma-seperated list of VCF files containing mutations from the same sample but not the same positions -- filelines merged in order of this list (lines not re-sorted after merge)",
        },
	],
};


sub help_brief {                            # keep this to just a few words <---
    "Merge multiple VCFs - keep the quality scores from files in desc order"
}


sub help_synopsis {
<<'HELP';
Merge multiple VCFs - keep the FORMAT lines from files in desc order.
HELP
}

sub help_detail {                  # this is what the user will see with the longer version of help. <---
<<'HELP';
Merge multiple VCFs. For identical calls made by different algorithms, merge them, keeping the FORMAT/scores from the file that is listed first in the vcf_files string.
HELP
}

###############

sub execute {                               # replace with real execution logic.
    my $self = shift;

    my $vcf_files = $self->vcf_files;
    my $output_file = $self->output_file;

    my @vcffiles = split(",",$vcf_files);
    if (@vcffiles < 1){
        die ("requires multiple VCF files to be input (comma-sep)")
    }

    #output
    open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";

    #get headers by printing the first file
    my $inFh = IO::File->new( $vcffiles[0] ) || die "can't open file\n";
    while(my $line = $inFh->getline )
    {
	print OUTFILE "$line";
    }
    close($inFh);

    #add data from subsequent files if data does not exist in first file
    for(my $i=1; $i<@vcffiles; $i++){
        $inFh = IO::File->new( $vcffiles[$i] ) || die "can't open file\n";

        while(my $line = $inFh->getline ) {
		if ($line =~ /^#/){
			next;
		}
		print OUTFILE "$line";
	}
    }

    return 1;
}

