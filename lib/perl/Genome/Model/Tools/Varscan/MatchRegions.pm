
package Genome::Model::Tools::Varscan::MatchRegions;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# RunVarscan - Run Varscan somatic on two BAM files.
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	12/09/2009 by D.K.
#	MODIFIED:	12/09/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

## SET DEFAULT PARAMS ##
my $min_loh_size = 10;
my $min_loh_snps = 3;


class Genome::Model::Tools::Varscan::MatchRegions {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		loh_regions	=> { is => 'Text', doc => "File of LOH calls", is_optional => 0 },
		gene_regions 	=> { is => 'Text', doc => "File of gene regions that could be affected", is_optional => 0 },
		output_file 	=> { is => 'Text', doc => "Output file", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Identifies gene regions that could be affected by LOH regions"                 
}

sub help_synopsis {
    return <<EOS
Identify gene regions that could be affected by LOH regions
EXAMPLE:	gmt capture match-regions ...
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
	my $loh_regions = $self->loh_regions;
	my $gene_regions = $self->gene_regions;
	my $output_file = $self->output_file;

	## Define some placeholder variables
	my $loh_chrom = my $loh_start = my $loh_stop = my $loh_snps = 0;

	my %gene_regions_by_chrom = load_gene_regions($gene_regions);

	## Open the output file ##
	if($output_file)
	{
		open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
		print OUTFILE "chrom\tchr_start\tchr_stop\tnum_snps\tregion_size\n";
	}
	
	## Parse the variant file ##

	my $input = new FileHandle ($loh_regions);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my @lineContents = split(/\t/, $line);
		
		if(($lineContents[0] eq "chrom" || $lineContents[0] eq "ref_name"))
		{

		}
		else
		{
			(my $chrom, my $chr_start, my $chr_stop, my $loh_snps, my $loh_size) = split(/\t/, $line);
			
			if($gene_regions_by_chrom{$chrom})
			{
				my @regions = split(/\n/, $gene_regions_by_chrom{$chrom});
				
				foreach my $region (@regions)
				{
					my $overlap_type = "";
					
					(my $region_start, my $region_stop, my $region_name) = split(/\t/, $region);
					
					## Case 1: Gene region contained within LOH region ##
					
					if($region_start >= $chr_start && $region_stop <= $chr_stop)
					{
						$overlap_type = "contains";	
					}
					
					## Case 2: LOH region contained within Gene region ##
					
					elsif($chr_start >= $region_start && $chr_stop <= $region_stop)
					{
						$overlap_type = "within";
					}
					
					## Case 3: LOH and gene regions overlap ##
					
					elsif(($region_start >= $chr_start && $region_start <= $chr_stop) || ($region_start >= $chr_stop && $region_stop <= $chr_stop))
					{
						$overlap_type = "overlaps";
					}
					
					if($overlap_type)
					{
						print "$chrom\t$chr_start\t$chr_stop\t$loh_snps\t$loh_size\t";
						print "$overlap_type\t";
						print "$region_start\t$region_stop\t$region_name\n";
					}
				}
			}
		}
	}
	
	close($input);
	
	if($output_file)
	{
		close(OUTFILE);
	}
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



################################################################################################
# Process results - filter variants by type and into high/low confidence
#
################################################################################################

sub load_gene_regions
{
	my $infile = shift(@_);	

	my %regions_by_chrom = ();

	my $input = new FileHandle ($infile);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $chrom, my $chr_start, my $chr_stop, my $region_name) = split(/\t/, $line);
		$regions_by_chrom{$chrom} .= "\n" if($regions_by_chrom{$chrom});
		$regions_by_chrom{$chrom} .= "$chr_start\t$chr_stop\t$region_name";
	}
	
	close($input);

	return(%regions_by_chrom)
}


1;

