
package Genome::Model::Tools::Varscan::Readcounts;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Varscan::Somatic	Runs Varscan somatic pipeline on Normal/Tumor BAM files
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	12/09/2009 by D.K.
#	MODIFIED:	12/29/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Varscan::Readcounts {
	is => 'Genome::Model::Tools::Varscan',
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		bam_file	=> { is => 'Text', doc => "Path to BAM file", is_optional => 0 },
		samtools_path	=> { is => 'Text', doc => "Path to SAMtools executable", is_optional => 0, is_input => 1, default => "samtools" },
		variants_file	=> { is => 'Text', doc => "Path to variant positions file", is_optional => 0 },
		output_file	=> { is => 'Text', doc => "Path to output file" , is_optional => 0},
		reference        => { is => 'Text', doc => "Reference FASTA file for BAMs" , is_optional => 1, default_value => (Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa')},
		min_coverage	=> { is => 'Text', doc => "Minimum base coverage to report readcounts [8]" , is_optional => 1},
		min_base_qual	=> { is => 'Text', doc => "Minimum base quality to count a read [30]" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Run the Varscan readcounts tool"                 
}

sub help_synopsis {
    return <<EOS
Runs Varscan readcounts from BAM files
EXAMPLE:	gmt varscan readcounts --bam-file [sample.bam] --variants-file [variants.tsv] --output-file readcounts.txt ...
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 

EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;

	## Get required parameters ##
	my $bam_file = $self->bam_file;
	my $reference = $self->reference;
	my $variants_file = $self->variants_file;
	my $output_file = $self->output_file;
	my $min_base_qual = $self->min_base_qual;
	
	if(-e $bam_file)
	{
		## Prepare pileup commands ##
		
#		my $pileup = "samtools view -b -u -q 10 $bam_file | samtools pileup -f $reference -";
		my $pileup = $self->samtools_path . " mpileup -q 10 -f $reference $bam_file";
		my $cmd = $self->java_command_line("readcounts <\($pileup\) --variants-file $variants_file --output-file $output_file --min-base-qual $min_base_qual");
		print "RUN: $cmd\n";
		system($cmd);
	}
	else
	{
		die "Error: One of your BAM files doesn't exist!\n";
	}
	
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


1;

