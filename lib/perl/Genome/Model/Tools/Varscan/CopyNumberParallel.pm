
package Genome::Model::Tools::Varscan::CopyNumberParallel;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::Varscan::CopyNumberParallel {
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
        chromosome => {
            is => 'Text',
            doc => "Specify a single chromosome (optional)",
            is_optional => 1,
            is_input => 1,
        },
        output => {
            is => 'Text',
            doc => "Output file for copy number results",
            is_optional => 0,
            is_input => 1,
            is_output => 1,
        },
        target_regions => {
            is => 'Text',
            doc => "Optional target region(s) for limiting the BAM (e.g 5:1 or 6:11134-11158)",
            is_optional => 1,
            is_input => 1,
        },
        reference => {
            is => 'Text',
            doc => "Reference FASTA file for BAMs",
            is_optional => 0,
            example_values => ['/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa'],
        },
        heap_space => {
            is => 'Text',
            doc => "Megabytes to reserve for java heap [1000]" ,
            is_optional => 1,
            is_input => 1,
        },
        mapping_quality => {
            is => 'Text',
            doc => "Default minimum mapping quality" ,
            is_optional => 1,
            is_input => 1,
            default => 10,
        },
        skip_if_output_present => {
            is => 'Text',
            doc => "If set to 1, skip execution if output files exist",
            is_optional => 1,
            is_input => 1,
        },
        varscan_params => {
            is => 'Text',
            doc => "Parameters to pass to Varscan" ,
            is_optional => 1,
            is_input => 1,
            default => "--min-coverage 20 --min-segment-size 25 --max-segment-size 100",
        },
    ],
    has_param => [
        lsf_resource => {
            default_value => 'select[mem>4000 && tmp>1000] rusage[mem=4000]'
        },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Run the Varscan2-CopyNumber on a tumor-normal pair, parallel by chromosome"                 
}

sub help_synopsis {
    return <<EOS
This command runs the Varscan2-CopyNumber analysis on a tumor-normal pair in parallel by chromosome. 
EXAMPLE:	gmt varscan copy-number --normal-bam [Normal.bam] --tumor-bam [Tumor.bam] --output varscan_out/varScan.output.copynumber ...
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
	my $min_map_qual = $self->mapping_quality;
	my $reference = $self->reference;
	
	my $index_file = "$reference.fai";
	
	if(!(-e $index_file))
	{
		die "Index file for reference ($index_file) not found!\n";
	}

	## Get output directive ##
	my $output = "output";

	if($self->output)
	{
		$output = $self->output;
	}
	else
	{
		die "Please provide an output basename (--output) or output files\n";
	}


#	my $varscan_params = "--min-coverage 20 --min-segment-size 10 --max-segment-size 100"; #--min-coverage 8 --verbose 1
	my $varscan_params = $self->varscan_params if($self->varscan_params);

	## Check skip if output present ##
	
	if($self->skip_if_output_present)
	{

	}

	if(-e $normal_bam && -e $tumor_bam)
	{	
		## Get the flagstat ##
		print "Getting Flagstat...\n";	
		my %normal_flagstat = get_flagstat($normal_bam);
		my %tumor_flagstat = get_flagstat($tumor_bam);

		print "Computing Average Read Length...\n";
		my $normal_readlen = avg_read_len($normal_bam);
		my $tumor_readlen = avg_read_len($tumor_bam);
		
		## Determine the total unique GBP ##
		
		my $normal_unique_bp = ($normal_flagstat{'mapped'} - $normal_flagstat{'duplicates'}) * $normal_readlen;
		my $tumor_unique_bp = ($tumor_flagstat{'mapped'} - $tumor_flagstat{'duplicates'}) * $tumor_readlen;

		my $normal_tumor_ratio = $normal_unique_bp / $tumor_unique_bp;
		$normal_tumor_ratio = sprintf("%.4f", $normal_tumor_ratio);

		print "Normal: $normal_readlen ==> $normal_unique_bp\n";
		print "Tumor: $tumor_readlen ==> $tumor_unique_bp\n";
		print "Ratio: $normal_tumor_ratio\n";

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
				if(!$self->chromosome  || $chrom eq $self->chromosome)
				{
					if(substr($chrom, 0, 2) ne "GL")
					{
						print "$chrom\t";
						my $output_file = $output . ".$chrom";
						my $varscan_path = Genome::Model::Tools::Varscan->path_for_version($self->version);
						my $mpileup = $self->samtools_path . " mpileup -f $reference -B -q 10 -r $chrom:1 $normal_bam $tumor_bam";
	
	#					my $cmd = $self->command_line(" copynumber <\($mpileup\) $output_file --mpileup 1 $varscan_params");
						my $cmd = "bash -c \"java -jar $varscan_path copynumber <\($mpileup\) $output_file --mpileup 1 --data-ratio $normal_tumor_ratio $varscan_params\"";
		
						system("bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} -R\"select[mem>2000 && tmp>2000] rusage[mem=2000]\" -oo $output_file.log $cmd");
						
					}
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





###################################################
# get_flagstat - 
#
###################################################

sub get_flagstat
{
	my $bam_file = shift(@_);
	my $flagstat = "";
	
	if(-e "$bam_file.flagstat")
	{
		$flagstat = `cat $bam_file.flagstat`;	
	}
	else
	{
		$flagstat = `samtools flagstat $bam_file`;
	}


	## IF we got it, parse it ##

	if($flagstat)
	{
		my %cov_stats = ();	
		my @lines = split(/\n/, $flagstat);
		
		foreach my $line (@lines)
		{
			## Erase the plus zero ##
			
			$line =~ s/\s\+\s0//;
			
			(my $num_reads) = split(/\s+/, $line);
			my $category = $line;
			$category =~ s/$num_reads\s//;
			
			## Remove stuff with parentheses ##
			my $split_char = " \\(";
			($category) = split(/$split_char/, $category);
			
			$cov_stats{$category} = $num_reads if($category);			
		}
		
		return(%cov_stats);
	}
	
	return();
}






#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub avg_read_len
{
	my $FileName = shift(@_);
	return(100);
	my $len_sum = my $len_num = 0;
	my $read_seqs = `samtools view $FileName 2>/dev/null | head -10000 | cut -f 10`;

	my @lines = split(/\n/, $read_seqs);
	foreach my $line (@lines)
	{
		my $read_len = length($line);
		$len_sum += $read_len;
		$len_num++;
	}

	my $avg_readlen = $len_sum / $len_num if($len_num);
	$avg_readlen = sprintf("%.2f", $avg_readlen);
	return($avg_readlen);
}


1;

