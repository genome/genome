
package Genome::Model::Tools::Varscan::NormalizeCopySegments;     # rename this when you give the module file a different name <--

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

use Genome;

class Genome::Model::Tools::Varscan::NormalizeCopySegments {
    is => 'Genome::Model::Tools::Varscan',
    has => [
        basename=> {
            is => 'Text',
            doc => "Path and basename to varScan copynumber files e.g. dir/varScan.output",
            is_input => 1,
        },
        output=> {
            is => 'Text',
            doc => "Output file for copy number results",
            is_input => 1,
            is_output => 1,
        },
        chromosome => {
            is => 'Text',
            doc => "Optional chromosome to do one only",
            is_optional => 1,
            is_input => 1,
        },
        reference => {
            is => 'Text',
            doc => "Reference FASTA file for BAMs",
            is_optional => 0,
            example_values => [Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa'],
        },
        skip_if_output_present => {
            is => 'Text',
            doc => "If set to 1, skip execution if output files exist",
            is_optional => 1,
            is_input => 1,
        },
        min_depth => {
            is => 'Text',
            doc => "Minimum depth in both samples to include for segmentation",
            default => 20,
        },
        undo_sd => {
            is => 'Text',
            doc => "Remove change-points of less than this number of standard deviations",
            default => 4,
        },
    ],
    has_param => [
        lsf_resource => {
            default_value => 'select[tmp>1000] rusage[mem=4000]'
        },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Normalizes CBS segments based on GC content"                 
}

sub help_synopsis {
    return <<EOS
This command normalizes CBS segments based on GC content
EXAMPLE:	gmt varscan normalize-copy-segments --basename varscan_out/varScan.output.copynumber --output varscan_out/varScan.output.copynumber.cbs ...
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
	my $basename = $self->basename;
	my $output = $self->output;
	my $reference = $self->reference;
	my $index_file = "$reference.fai";
	
	if(!(-e $index_file))
	{
		die "Index file for reference ($index_file) not found!\n";
	}

	## Check skip if output present ##
	
	if($self->skip_if_output_present)
	{

	}

	if(-e $index_file)
	{	
		## Get the image basename ##
		my @tempArray = split(/\//, $output);
		my $numElements = @tempArray;
		my $image_basename = $tempArray[$numElements - 1];

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
						my $input_file = $basename . ".$chrom.infile.segments.p_value";
						
						if(-e $input_file)
						{
							my $output_basename = $output . ".$chrom";
							parse_regions($input_file, $output_basename, $chrom, $self);
						}
						
					}
				}


			}

		}
	
		close($input);
		
		close(HTML);

	}
	else
	{
		die "Error: One of your BAM files doesn't exist!\n";
	}

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




################################################################################################
# Execute - the main program logic
#
################################################################################################

sub parse_regions
{                               # replace with real execution logic.
	my $segments_file = shift(@_);
	my $output_basename = shift(@_);
	my $chrom = shift(@_);
	my $self = shift(@_);
	my $reference = $self->reference;
	
	my $result = "";

	my $chrom_filename = $output_basename . ".infile";
	my $normalized_file = "$output_basename.$chrom.tsv";

#	open(OUTFILE, ">$normalized_file") or die "Can't open outfile: $!\n";

	my $input = new FileHandle ($segments_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		if($lineCounter > 1)
		{
			my ($id, $name, $chrom, $chr_start, $chr_stop, $num_positions, $seg_mean, $bstat, $pval, $lcl, $ucl) = split(/\s+/, $line);

			## Retrieve the GC content of the sequence ##
			my $query_string = $chrom . ":" . $chr_start . "-" . $chr_stop;
			my $sequence = `samtools faidx $reference $query_string | grep -v \">\"`;
			chomp($sequence);
		}
	}
	
	close($input);

	print "$lineCounter regions\n";
#	close(OUTFILE);
}


1;

