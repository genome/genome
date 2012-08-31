
package Genome::Model::Tools::Varscan::Consensus;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::Varscan::Consensus {
	is => 'Genome::Model::Tools::Varscan',
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		bam_file	=> { is => 'Text', doc => "Path to BAM file", is_optional => 0 },
		positions_file	=> { is => 'Text', doc => "Path to variant positions file", is_optional => 0 },
		output_file	=> { is => 'Text', doc => "Path to output file" , is_optional => 0},
		min_coverage	=> { is => 'Text', doc => "Minimum base coverage to report readcounts [4]" , is_optional => 1},
		min_avg_qual	=> { is => 'Text', doc => "Minimum base quality to count a read [20]" , is_optional => 1},
		min_var_freq	=> { is => 'Text', doc => "Minimum variant allele frequency to call a variant [0.20]" , is_optional => 1},
		reference        => { is => 'Text', doc => "Reference FASTA file for BAMs" , is_optional => 1, default_value => (Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa')},		
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Run the Varscan pileup2cns tool"                 
}

sub help_synopsis {
    return <<EOS
Runs Varscan readcounts from BAM files
EXAMPLE:	gmt varscan consensus --bam-file [sample.bam] --variants-file [variants.tsv] --output-file readcounts.txt ...
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
	my $positions_file = $self->positions_file;
	my $output_file = $self->output_file;
	
	my $min_coverage = 4;
	my $min_avg_qual = 20;
	my $min_var_freq = 0.20;
	
	$min_coverage = $self->min_coverage if(defined($self->min_coverage));
	$min_avg_qual = $self->min_avg_qual if(defined($self->min_avg_qual));
	$min_var_freq = $self->min_var_freq if(defined($self->min_var_freq));

	if(-e $bam_file)
	{
		## Prepare pileup commands ##
		print "Building ROI pileup file...\n";
		my $cmd = "samtools view -b -u -q 10 $bam_file | samtools pileup -f $reference - | java -classpath ~dkoboldt/Software/Varscan net.sf.varscan.Varscan limit --positions-file $positions_file --output-file $output_file.pileup";
		system($cmd);
		print "Running Varscan pileup2cns...\n";
		$cmd = $self->java_command_line("pileup2cns $output_file.pileup --min-coverage $min_coverage --min-var-freq $min_var_freq --min-avg-qual $min_avg_qual >$output_file");
		system($cmd);
	}
	else
	{
		die "Error: One of your BAM files doesn't exist!\n";
	}
	
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


1;

