
package Genome::Model::Tools::Varscan::SomaticParallelFilter;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Varscan::Somatic	Runs VarScan somatic pipeline on Normal/Tumor BAM files
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

class Genome::Model::Tools::Varscan::SomaticParallelFilter {
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
            doc => "Reference FASTA file for BAMs",
            is_optional => 0,
            example_values => [(Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa')],
        },
        chromosome => {
            is => 'Text',
            doc => "Specify a single chromosome (optional)",
            is_optional => 1,
            is_input => 1,
        },
        filter_loh => {
            is => 'Text',
            doc => "If set to 1, apply filter to LOH-HC calls using normal BAM",
            is_optional => 1,
            is_input => 1,
            default => 1,
        },
        filter_germline => {
            is => 'Text',
            doc => "If set to 1, apply filter to Germline-HC calls using tumor BAM",
            is_optional => 1,
            is_input => 1,
            default => 0,
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
            doc => "Parameters to pass to VarScan [--min-coverage 3 --min-var-freq 0.08 --p-value 0.10 --somatic-p-value 0.05 --strand-filter 1]" ,
            is_optional => 1,
            is_input => 1,
        },
    ],
    has_param => [
        lsf_resource => {
            default_value => 'select[rusage[mem=4000]'
        },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Applies filters to the VarScan calls"                 
}

sub help_synopsis {
    return <<EOS
Runs VarScan from BAM files
EXAMPLE:	gmt varscan somatic-parallel-filter --normal-bam [Normal.bam] --tumor-bam [Tumor.bam] --output varscan_out/Patient.status ...
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

	## Get VarScan parameters ##

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

			if($chrom =~ 'NT_' || $chrom =~ /GL/)
			{
#				print "Skipping $chrom\n";								
			}
			else
			{
				if(!$self->chromosome || $chrom eq $self->chromosome)
				{
					print "CHROMOSOME $chrom\n";
					$output_snp = $output . ".$chrom.snp";
					$output_indel = $output . ".$chrom.indel";
	
					## FORMAT SNVS ##
	
					if(-e $output_snp)
					{
						if($self->skip_if_output_present && -e "$output_snp.formatted")
						{
							## Skip because already present 
						}
						else
						{
							print "Formatting SNVs...\n";
							## Format the variants ##
							my $cmd_obj = Genome::Model::Tools::Capture::FormatSnvs->create(
							    variants_file => $output_snp,
							    output_file => "$output_snp.formatted",
							    append_line => 1,
							);
							
							$cmd_obj->execute;                                                                                        
						}
	
					}
					else
					{
						warn "Warning: Missing SNP file $output_snp\n";
					}
	
					## FORMAT INDELS ##
					
					if(-e $output_indel)
					{
						if($self->skip_if_output_present && -e "$output_snp.formatted")
						{
							## Skip because already present 
						}
						else
						{
							print "Formatting indels...\n";
							## Format the variants ##
							my $cmd_obj = Genome::Model::Tools::Capture::FormatIndels->create(
							    variants_file => $output_indel,
							    output_file => "$output_indel.formatted",
							);
							
							$cmd_obj->execute;
						}
					}
					else
					{
						warn "Warning: Missing indel file $output_indel\n";
						
					}
	
					## PROCESS SNVS ##
	
					if(-e "$output_snp.formatted")
					{
						if($self->skip_if_output_present && -e "$output_snp.formatted.Somatic")
						{
							## Skip because already present 
						}
						else
						{
							print "Processing SNVs...\n";
							my $cmd_obj = Genome::Model::Tools::Varscan::ProcessSomatic->create(
							    status_file => "$output_snp.formatted",
							);
							
							$cmd_obj->execute;                                                
						}                                        
					}
					
					## PROCESS INDELS ##
	
					if(-e "$output_indel.formatted")
					{
						if($self->skip_if_output_present && -e "$output_indel.formatted.Somatic")
						{
							## Skip because already present 
						}
						else
						{
							print "Processing Indels...\n";
	
							my $cmd_obj = Genome::Model::Tools::Varscan::ProcessSomatic->create(
							    status_file => "$output_indel.formatted",
							);
							
							$cmd_obj->execute;                                                
						}                                        
					}
	
	
					## FILTER SNVS ##
	
					# Somatic #                                
					my $variant_file = "$output_snp.formatted.Somatic.hc";
					my $bam_file = $tumor_bam;
					run_filter($self, $variant_file, $bam_file);
					# Germline #
					$variant_file = "$output_snp.formatted.Germline.hc";
					$bam_file = $tumor_bam;
					run_filter($self, $variant_file, $bam_file) if($self->filter_germline);
					# LOH using normal bam ##
					$variant_file = "$output_snp.formatted.LOH.hc";
					$bam_file = $normal_bam;
					run_filter($self, $variant_file, $bam_file) if($self->filter_loh);
	
	
					## FILTER INDELS ##
	
					# Somatic #                                
					$variant_file = "$output_indel.formatted.Somatic.hc";
					$bam_file = $tumor_bam;
					run_indel_filter($self, $variant_file, $bam_file);
					# Germline #
					$variant_file = "$output_indel.formatted.Germline.hc";
					$bam_file = $tumor_bam;
					run_indel_filter($self, $variant_file, $bam_file) if($self->filter_germline);
					# LOH using normal bam ##
					$variant_file = "$output_indel.formatted.LOH.hc";
					$bam_file = $normal_bam;
					run_indel_filter($self, $variant_file, $bam_file) if($self->filter_loh);
					
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




################################################################################################
# Run Filter - Launches the LSF job to filter the variants 
#
################################################################################################

sub run_filter
{                               # replace with real execution logic.
	my ($self, $variant_file, $bam_file) = @_;
        
        if(-e $variant_file)
        {
                if($self->skip_if_output_present && -e "$variant_file.fpfilter")
                {
                        ## Skip because output present ##
                }
                else
                {
                        my $cmd = "gmt somatic filter-false-positives --variant-file $variant_file --bam-file $bam_file --output-file $variant_file.fpfilter --filtered-file $variant_file.fpfilter.removed";
			if($self->reference)
			{
				$cmd .= " --reference " . $self->reference;
			}
                        system("bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} -oo $variant_file.err -J varscan -R\"select[mem>2000 && tmp>2000] rusage[mem=2000]\" $cmd");
                }
        }        
}


################################################################################################
# Run Filter - Launches the LSF job to filter the variants 
#
################################################################################################

sub run_indel_filter
{                               # replace with real execution logic.
	my ($self, $variant_file, $bam_file) = @_;
        
        if(-e $variant_file)
        {
                if($self->skip_if_output_present && -e "$variant_file.fpfilter")
                {
                        ## Skip because output present ##
                }
                else
                {
                        my $cmd = "gmt somatic filter-false-indels --variant-file $variant_file --bam-file $bam_file --output-file $variant_file.fpfilter --filtered-file $variant_file.fpfilter.removed";
                        system("bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} -oo $variant_file.err -J varscan -R\"select[mem>2000 && tmp>2000] rusage[mem=2000]\" $cmd");
                }
        }        
}

1;

