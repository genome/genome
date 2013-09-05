
package Genome::Model::Tools::Varscan::Consensus;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Genome::Model::Tools::Varscan::Consensus	Runs Varscan pileup2cns on the SAMtools mpileup output from one BAM file.
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	12/09/2009 by D.K.
#	MODIFIED:	08/14/2013 by D.K.
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
        bam_file => {
            is => 'Text',
            doc => "Path to BAM file",
            is_optional => 0,
        },
        output_file => {
            is => 'Text',
            doc => "Path to output file",
            is_optional => 0,
        },
        min_coverage => {
            is => 'Text',
            doc => "Minimum base coverage to report readcounts",
	    default => 3,
            is_optional => 1,
        },
        min_avg_qual => {
            is => 'Text',
            doc => "Minimum base quality to count a read",
	    default => 20,
            is_optional => 0,
        },
        min_var_freq => {
            is => 'Text',
            doc => "Minimum variant allele frequency to call a variant",
	    default => 0.20,
            is_optional => 0,
        },
        output_vcf => {
            is => 'Text',
            doc => "If set to 1, tells VarScan to output in VCF format (rather than native CNS)",
            is_optional => 1,
        },
        position_list_file => {
            is => 'Text',
            doc => "Optionally, provide a tab-delimited list of positions to be given to SAMtools with -l",
            is_optional => 1,
        },
        reference => {
            is => 'Text',
            doc => "Reference FASTA file for BAMs",
            is_optional => 0,
            example_values => ['/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa'],
        },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Run VarScan consensus calling for one BAM file"                 
}

sub help_synopsis {
    return <<EOS
Runs mpileup and then VarScan consensus calling (pileup2cns) on a single BAM file
EXAMPLE:	gmt varscan consensus --bam-file sample.bam --reference reference.fa --output sample.bam.varScan.cns
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
	my $output_file = $self->output_file;
	
	my $min_coverage = $self->min_coverage;
	my $min_avg_qual = $self->min_avg_qual;
	my $min_var_freq = $self->min_var_freq;

	if(-e $bam_file)
	{
		## Prepare pileup commands ##
		my $mpileup = $self->samtools_path . " mpileup -B -f $reference -q 10 $bam_file";
		
		if($self->position_list_file)
		{
			$mpileup = $self->samtools_path . " mpileup -B -f $reference -q 10 -l " . $self->position_list_file . " " . $bam_file;
		}
		
		
		
		my $cmd = "";
		
		if($self->output_vcf)
		{
			$cmd = $self->java_command_line("mpileup2cns <\($mpileup\) --min-coverage $min_coverage --min-var-freq $min_var_freq --min-avg-qual $min_avg_qual --output-vcf 1 >$output_file 2>/dev/null");			
		}
		else
		{
			$cmd = $self->java_command_line("pileup2cns <\($mpileup\) --min-coverage $min_coverage --min-var-freq $min_var_freq --min-avg-qual $min_avg_qual >$output_file 2>/dev/null");
		}

		system($cmd);
	}
	else
	{
		die "Error: One of your BAM files doesn't exist!\n";
	}
	
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


1;

