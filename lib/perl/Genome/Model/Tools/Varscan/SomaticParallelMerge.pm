
package Genome::Model::Tools::Varscan::SomaticParallelMerge;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::Varscan::SomaticParallelMerge {
	is => 'Genome::Model::Tools::Varscan',
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		normal_bam	=> { is => 'Text', doc => "Path to Normal BAM file", is_optional => 0, is_input => 1 },
		tumor_bam	=> { is => 'Text', doc => "Path to Tumor BAM file", is_optional => 0, is_input => 1 },
		output	=> { is => 'Text', doc => "Path to Tumor BAM file", is_optional => 1, is_input => 1, is_output => 1 },
		output_snp	=> { is => 'Text', doc => "Basename for SNP output, eg. varscan_out/varscan.status.snp" , is_optional => 1, is_input => 1, is_output => 1},
		output_indel	=> { is => 'Text', doc => "Basename for indel output, eg. varscan_out/varscan.status.snp" , is_optional => 1, is_input => 1, is_output => 1},
		reference        => { is => 'Text', doc => "Reference FASTA file for BAMs" , is_optional => 1, default_value => (Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa')},
		chromosome	=> { is => 'Text', doc => "Specify a single chromosome (optional)", is_optional => 1, is_input => 1},
		heap_space	=> { is => 'Text', doc => "Megabytes to reserve for java heap [1000]" , is_optional => 1, is_input => 1},
		skip_if_output_present	=> { is => 'Text', doc => "If set to 1, skip execution if output files exist", is_optional => 1, is_input => 1 },
		varscan_params	=> { is => 'Text', doc => "Parameters to pass to VarScan [--min-coverage 3 --min-var-freq 0.08 --p-value 0.10 --somatic-p-value 0.05 --strand-filter 1]" , is_optional => 1, is_input => 1},
	],	

	has_param => [
		lsf_resource => { default_value => 'select[model!=Opteron250 && type==LINUX64] rusage[mem=4000]'},
       ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Merges genome-wide filtered VarScan calls"                 
}

sub help_synopsis {
    return <<EOS
Runs VarScan from BAM files
EXAMPLE:	gmt varscan somatic-parallel-merge --normal-bam [Normal.bam] --tumor-bam [Tumor.bam] --output varscan_out/Patient.status ...
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

	my %stats = ();

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
		## Remove output files if they exist ##
		remove_output_files($output, "formatted.Somatic");
		remove_output_files($output, "formatted.LOH");
		remove_output_files($output, "formatted.Germline");

		remove_output_files($output, "formatted.Somatic.hc");
		remove_output_files($output, "formatted.LOH.hc");
		remove_output_files($output, "formatted.Germline.hc");

		remove_output_files($output, "formatted.Somatic.hc.fpfilter");
		remove_output_files($output, "formatted.LOH.hc.fpfilter");
		remove_output_files($output, "formatted.Germline.hc.fpfilter");		
		
		print "CHROM\tSOMATIC\t\tLOH\t\tGERMLINE\t\tHC_SOM\t\tHC_LOH\t\tHC_GERM\t\tFILT_SOM\t\tFILT_LOH\t\tFILT_GERM\n";
		print "chrom\tSNV\tIndel\tSNV\tIndel\tSNV\tIndel\tSNV\tIndel\tSNV\tIndel\tSNV\tIndel\tSNV\tIndel\tSNV\tIndel\tSNV\tIndel\n";

		my $input = new FileHandle ($index_file);
		my $lineCounter = 0;
	
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;
			
			my ($chrom) = split(/\t/, $line);

			if($chrom =~ 'NT_')
			{
#				print "Skipping $chrom\n";								
			}
			else
			{
				if(!$self->chromosome || $chrom eq $self->chromosome)
				{
					## Check for file truncation ##
					check_files($output, $chrom);
					my $chrom_output_line = $chrom;
	
					## Merge the output files ##
					my $num_snp = my $num_indel = 0;
	
					# All Calls, Formatted for Annotation ##
					($num_snp, $num_indel) = merge_output($chrom, $output, "formatted");
	
					# All Somatic ##
					($num_snp, $num_indel) = merge_output($chrom, $output, "formatted.Somatic");
					$chrom_output_line .= "\t$num_snp\t$num_indel";
					# All LOH ##
					($num_snp, $num_indel) = merge_output($chrom, $output, "formatted.LOH");
					$chrom_output_line .= "\t$num_snp\t$num_indel";
					# All Germline ##
					($num_snp, $num_indel) = merge_output($chrom, $output, "formatted.Germline");
					$chrom_output_line .= "\t$num_snp\t$num_indel";
	
	
					# HC Somatic ##
					($num_snp, $num_indel) = merge_output($chrom, $output, "formatted.Somatic.hc");
					$chrom_output_line .= "\t$num_snp\t$num_indel";
					# HC LOH ##
					($num_snp, $num_indel) = merge_output($chrom, $output, "formatted.LOH.hc");
					$chrom_output_line .= "\t$num_snp\t$num_indel";
					# HC Germline ##
					($num_snp, $num_indel) = merge_output($chrom, $output, "formatted.Germline.hc");
					$chrom_output_line .= "\t$num_snp\t$num_indel";
	
					
					# HC Somatic passing filter ##
					($num_snp, $num_indel) = merge_output($chrom, $output, "formatted.Somatic.hc.fpfilter");
					$chrom_output_line .= "\t$num_snp\t$num_indel";
					# HC LOH passing filter ##
					($num_snp, $num_indel) = merge_output($chrom, $output, "formatted.LOH.hc.fpfilter");
					$chrom_output_line .= "\t$num_snp\t$num_indel";
					# HC Germline passing filter ##
					($num_snp, $num_indel) = merge_output($chrom, $output, "formatted.Germline.hc.fpfilter");
					$chrom_output_line .= "\t$num_snp\t$num_indel";
	
	
					# HC Somatic failing filter ##
					($num_snp, $num_indel) = merge_output($chrom, $output, "formatted.Somatic.hc.fpfilter.removed");
					# HC LOH failing filter ##
					($num_snp, $num_indel) = merge_output($chrom, $output, "formatted.LOH.hc.fpfilter.removed");
					# HC Germline failing filter ##
					($num_snp, $num_indel) = merge_output($chrom, $output, "formatted.Germline.hc.fpfilter.removed");
					
					print "$chrom_output_line\n";					
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
# Remove output file - remove a file if it exists 
#
################################################################################################

sub check_files
{
	my ($output, $chrom) = @_;

	my $output_snp = $output . ".$chrom.snp";
	my $output_indel = $output . ".$chrom.indel";


	## Check for incomplete files ##

	my $num_snv = my $num_snv_formatted = 0;
	$num_snv = `cat $output_snp | wc -l`;
	$num_snv_formatted = `cat $output_snp.formatted | wc -l`;
	
	if(($num_snv - $num_snv_formatted) > 1)
	{
		## Truncated SNV file ##
		warn "$output_snp.formatted (lines=$num_snv_formatted) is smaller than $output_snp ($num_snv lines) and may be truncated!\n";
	}
	
	my $num_indel = my $num_indel_formatted = 0;
	$num_indel = `cat $output_indel | wc -l`;
	$num_indel_formatted = `cat $output_indel.formatted | wc -l`;
	
	if(($num_indel - $num_indel_formatted) > 1)
	{
		## Truncated SNV file ##
		warn "$output_indel.formatted (lines=$num_indel_formatted) is smaller than $output_indel ($num_indel lines) and may be truncated!\n";
	}	
}


################################################################################################
# Remove output file - remove a file if it exists 
#
################################################################################################

sub merge_output
{
	my ($chrom, $output, $suffix) = @_;

	my $chrom_file_snp = $output . ".$chrom.snp.$suffix" ;
	my $chrom_file_indel = $output . ".$chrom.indel.$suffix";

	my $output_file_snp = $output . ".snp.$suffix";
	my $output_file_indel = $output . ".indel.$suffix";

	my $num_snps = `cat $chrom_file_snp | wc -l`;
	chomp($num_snps);
	
	my $num_indels = `cat $chrom_file_indel | wc -l`;
	chomp($num_indels);

	## Append the files ##
	
	system("cat $chrom_file_snp >>$output_file_snp");
	system("cat $chrom_file_indel >>$output_file_indel");

	return($num_snps, $num_indels);
}



################################################################################################
# Remove output file - remove a file if it exists 
#
################################################################################################

sub remove_output_files
{                               # replace with real execution logic.
	my ($output, $suffix) = @_;
        
	my $output_file_snp = $output . ".snp.$suffix";
	my $output_file_indel = $output . ".indel.$suffix";

	system("rm -rf $output_file_snp") if (-e $output_file_snp);
	system("rm -rf $output_file_indel") if (-e $output_file_indel);	
}

1;

