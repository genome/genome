
package Genome::Model::Tools::Capture::SampleRegionCoverage;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SampleRegionCoverageForAnnotation - Merge glfSomatic/Varscan somatic calls in a file that can be converted to MAF format
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	10/23/2009 by D.K.
#	MODIFIED:	10/23/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Capture::SampleRegionCoverage {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		q20_files	=> { is => 'Text', doc => "Path to Pileup-Q20 files for ROI, comma-delimited", is_optional => 0, is_input => 1 },
		regions_file	=> { is => 'Text', doc => "Path to regions file in chrom start stop format", is_optional => 0, is_input => 1 },
		min_depth	=> { is => 'Text', doc => "Minimum depth of coverage to assess [20]", is_optional => 0, is_input => 1, default => 20 },
		output_file     => { is => 'Text', doc => "Output file to receive coverage report", is_optional => 0, is_input => 1, is_output => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Reports Q20 coverage for one or more samples across a set of target regions"                 
}

sub help_synopsis {
    return <<EOS
This command reports Q20 coverage for one or more samples across a set of target regions.
It's useful for collaborations like RP where we need to know how many exons in a linkage region were not covered
EXAMPLE:	gmt capture sample-region-coverage --q20-files sample1.pileup.roi.q20,sample2.pileup.roi.q20 --regions-file linkage-region.tsv --output-file Region-Coverage-Samples1-2.tsv
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
	my $regions_file = $self->regions_file;
	my $q20_files = $self->q20_files;
	my $output_file = $self->output_file;
	my $min_depth = $self->min_depth;
	
	
	my @q20_files = split(/\,/, $q20_files);
	
	my %master_coverage = ();
	
	my $num_samples = 0;
	foreach my $filename (@q20_files)
	{
		warn "Loading $filename...\n";
		$num_samples++;
		my %coverage = load_q20_file($filename);
		
		foreach my $key (keys %coverage)
		{
			$master_coverage{$filename}{$key} = $coverage{$key};
		}
	}

	
	print "Parsing $regions_file...\n";
	
	## Open outfile ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	print OUTFILE "chrom\tchr_start\tchr_stop";
	
	foreach my $filename (@q20_files)
	{
		print OUTFILE "\t$filename";
	}

	print OUTFILE "\n";
	
	## Parse the indels ##

	my $input = new FileHandle ($regions_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my ($chrom, $chr_start, $chr_stop) = split(/\t/, $line);
		
		my $region_bases = 0;
		my %region_samples_covered = ();
		$region_samples_covered{'0x'} = $region_samples_covered{'1x'} = $region_samples_covered{'10x'} = $region_samples_covered{'15x'} = $region_samples_covered{'20x'} = 0;
		
		## Get a base coverage for each file ##
		
		my %bases_covered = ();
		
		## Go through each filename ##
		my $sample_results = "";
		
		foreach my $filename (@q20_files)
		{
			my $num_bases = my $num_covered = 0;

			## Go through each position in the region ##
			for(my $position = $chr_start; $position <= $chr_stop; $position++)
			{			
				my $key = join("\t", $chrom, $position);
				my $coverage = $master_coverage{$filename}{$key};
				$num_bases++;
				$num_covered++ if($coverage && $coverage >= $min_depth);
			}
			
			## Calculate the pct bases covered at 20x ##
			
			my $pct_covered = 0.000;
			$pct_covered = ($num_covered / $num_bases) if($num_bases);
			$pct_covered = sprintf("%.3f", $pct_covered);
			
			## Append to sample results ##
			
			$sample_results .= "\t" if($sample_results);
			$sample_results .= "$pct_covered";
		}
		
		print OUTFILE join("\t", $chrom, $chr_start, $chr_stop, $sample_results) . "\n";

	}

	close($input);
	
	
	close(OUTFILE);
}


################################################################################################
# Load Q20 File - load Q20 coverage for a file 
#
################################################################################################

sub load_q20_file
{
	my $filename = shift(@_);
	my %coverage = ();
	
	## Parse the indels ##

	my $input = new FileHandle ($filename);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my ($chrom, $pos, $q20) = split(/\t/, $line);

		my $key = join("\t", $chrom, $pos);
		
		$coverage{$key} = $q20;
	}

	close($input);	

	return(%coverage);
}


1;

