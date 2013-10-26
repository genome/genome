
package Genome::Model::Tools::Varscan::SomaticParallel;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::Varscan::SomaticParallel {
    is => 'Genome::Model::Tools::Varscan',

    has => [                                # specify the command's single-value properties (parameters) <--- 
        normal_bam => {
            is => 'Text',
            doc => "Path to Normal BAM file",
            is_optional => 0,
            is_input => 1,
        },
        tumor_bam => {
            is => 'Text',
            doc => "Path to Tumor BAM file",
            is_optional => 0,
            is_input => 1,
        },
        output => {
            is => 'Text',
            doc => "Path to Tumor BAM file",
            is_optional => 1,
            is_input => 1,
            is_output => 1,
        },
        output_snp => {
            is => 'Text',
            doc => "Basename for SNP output, eg. varscan_out/varscan.status.snp",
            is_optional => 1,
            is_input => 1,
            is_output => 1,
        },
        output_indel => {
            is => 'Text',
            doc => "Basename for indel output, eg. varscan_out/varscan.status.snp",
            is_optional => 1,
            is_input => 1,
            is_output => 1,
        },
        reference => {
            is => 'Text',
            doc => "Reference FASTA file for BAMs" ,
            is_optional => 0,
            example_values => [(Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa')],
        },
        chromosome => {
            is => 'Text',
            doc => "Specify a single chromosome (optional)",
            is_optional => 1,
            is_input => 1,
        },
        heap_space => {
            is => 'Text',
            doc => "Megabytes to reserve for java heap [1000]",
            is_optional => 1,
            is_input => 1,
        },
        skip_if_output_present => {
            is => 'Text',
            doc => "If set to 1, skip execution if output files exist",
            is_optional => 1,
            is_input => 1,
        },
        varscan_params => {
            is => 'Text',
            doc => "Parameters to pass to Varscan [--min-coverage 3 --min-var-freq 0.08 --p-value 0.10 --somatic-p-value 0.05 --strand-filter 1]",
            is_optional => 1,
            is_input => 1,
        },
    ],
    has_param => [
        lsf_resource => {
            default_value => 'rusage[mem=4000]'
        },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Run the Varscan somatic variant detection"                 
}

sub help_synopsis {
    return <<EOS
Runs Varscan from BAM files
EXAMPLE:	gmt varscan somatic-parallel --normal-bam [Normal.bam] --tumor-bam [Tumor.bam] --output varscan_out/Patient.status ...
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
	my $normal_bam = $self->normal_bam;
	my $tumor_bam = $self->tumor_bam;

	## Get output directive ##
	my $output = my $output_snp = my $output_indel = "";

	if($self->output)
	{
		$output = $self->output;
	}
	else
	{
		die "Please provide an output basename (--output) or output files for SNPs (--output-snp) and indels (--output-indels)\n";
	}


	my $reference = $self->reference;

	my $index_file = "$reference.fai";
	
	if(!(-e $index_file))
	{
		die "Index file for reference ($index_file) not found!\n";
	}

	## Get Varscan parameters ##

	my $varscan_params = "--min-coverage 3 --min-var-freq 0.08 --p-value 0.10 --somatic-p-value 0.05 --strand-filter 1"; #--min-coverage 8 --verbose 1
	$varscan_params = $self->varscan_params if($self->varscan_params);

	## Check skip if output present ##
	
	if($self->skip_if_output_present)
	{
		if(-e $output_snp)
		{
			my $snp_len = `cat $output_snp | wc -l`;
			chomp($snp_len);
			if($snp_len > 1)
			{
				return 1;
			}
		}
	}

	if(-e $normal_bam && -e $tumor_bam)
	{	
		my $input = new FileHandle ($index_file);
		my $lineCounter = 0;
	
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;
			
			my ($chrom) = split(/\t/, $line);

			if($chrom =~ 'NT_' || $chrom =~ /GL/) #skipping nonassembled contigs (GL is for build37)
			{
#				print "Skipping $chrom\n";								
			}
			else
			{
				if(!$self->chromosome  || $chrom eq $self->chromosome)
				{
					print "$chrom\t";
					$output_snp = $output . ".$chrom.snp";
					$output_indel = $output . ".$chrom.indel";
					print "$output_snp\t$output_indel\n";
			
	#                                my $normal_pileup = "samtools view -b -u -q 10 $normal_bam $chrom:1 | samtools pileup -f $reference -";
	 #                               my $tumor_pileup = "samtools view -b -u -q 10 $tumor_bam $chrom:1 | samtools pileup -f $reference -";
					my $normal_pileup = $self->samtools_path . " mpileup -f $reference -q 10 -r $chrom:1 $normal_bam";
					my $tumor_pileup = $self->samtools_path . " mpileup -f $reference -q 10 -r $chrom:1 $tumor_bam";
							
					my $cmd = $self->command_line(" somatic <\($normal_pileup\) <\($tumor_pileup\) --output-snp $output_snp --output-indel $output_indel $varscan_params");
	
					print "Running $cmd\n";                
					system("bsub -q long -J varscan -R\"select[mem>2000 && tmp>2000] rusage[mem=2000]\" $cmd");
	 #                               system($cmd);					
				}


			}

		}
	
		close($input);

	}
	else
	{
		die "Error: One of your BAM files doesn't exist!\n";
	}
	
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


1;

