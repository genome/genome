
package Genome::Model::Tools::Varscan::LohRegions;     # rename this when you give the module file a different name <--

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
my $num_loh_regions = 0;

class Genome::Model::Tools::Varscan::LohRegions {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		status_file	=> { is => 'Text', doc => "Path to file of Germline/Somatic/LOH calls", is_optional => 0 },
		output_file 	=> { is => 'Text', doc => "Output file for LOH regions", is_optional => 0 },
		min_loh_size 	=> { is => 'Text', doc => "Minimum size in bp for LOH region [10]", is_optional => 1 },
		min_loh_snps 	=> { is => 'Text', doc => "Minimum # of LOH SNPs for LOH region [3]", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Call regions of LOH using Varscan somatic output"                 
}

sub help_synopsis {
    return <<EOS
Call regions of LOH using Varscan somatic output
EXAMPLE:	gmt capture loh-regions ...
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
	my $status_file = $self->status_file;
	my $output_file = $self->output_file;

	$min_loh_size = $self->min_loh_size if($self->min_loh_size);
	$min_loh_snps = $self->min_loh_snps if($self->min_loh_snps);
	

	## Define some placeholder variables
	my $loh_chrom = my $loh_start = my $loh_stop = my $loh_snps = 0;

	## Open the output file ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	print OUTFILE "chrom\tchr_start\tchr_stop\tnum_snps\tregion_size\n";

	## Parse the variant file ##

	my $input = new FileHandle ($status_file);
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
			my $chrom = $lineContents[0];
			my $position = $lineContents[1];
			my $somatic_status = $lineContents[12];
			my $germline_p_value = $lineContents[13];
			my $somatic_p_value= $lineContents[14];
			if($lineContents[13] eq "Somatic" || $lineContents[13] eq "LOH" || $lineContents[13] eq "Germline")
			{
				$somatic_status = $lineContents[13];
				$germline_p_value = $lineContents[14];
				$somatic_p_value= $lineContents[15];
			}	
			if($chrom ne $loh_chrom)
			{
				report_loh_region($loh_chrom, $loh_start, $loh_stop, $loh_snps);
				
				$loh_chrom = $loh_start = $loh_stop = $loh_snps = 0;
			}
	
			if($somatic_status eq "LOH")
			{
				if(!$loh_snps)
				{
					$loh_chrom = $chrom;
					$loh_start = $position;
				}
				
				$loh_stop = $position;
				$loh_snps++;
			}
			elsif($somatic_status eq "Germline" && $germline_p_value <= 1.0E-06)
			{
				report_loh_region($loh_chrom, $loh_start, $loh_stop, $loh_snps);
				$loh_chrom = $loh_start = $loh_stop = $loh_snps = 0;
			}
		}
	}
	
	close($input);
		
	print "$num_loh_regions LOH regions identified\n";
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



################################################################################################
# Process results - filter variants by type and into high/low confidence
#
################################################################################################

sub report_loh_region
{
	(my $loh_chrom, my $loh_start, my $loh_stop, my $loh_snps) = @_;

	my $loh_size = $loh_stop - $loh_start + 1;
	
	if($loh_snps >= $min_loh_snps && $loh_size >= $min_loh_size)
	{
		print OUTFILE "$loh_chrom\t$loh_start\t$loh_stop\t$loh_snps\t$loh_size\n";	
		#print "$loh_chrom\t$loh_start\t$loh_stop\t$loh_snps\t$loh_size\n";
		$num_loh_regions++;
	}

}


1;

