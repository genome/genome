
package Genome::Model::Tools::Analysis::Sammy::AdjustCoverage;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# AdjustCoverage - Call somatic variants from normal/tumor BAM files
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	07/28/2009 by D.K.
#	MODIFIED:	07/28/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::Sammy::AdjustCoverage {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		coverage_file	=> { is => 'Text', doc => "Existing coverage file from Sammy3 coverage", is_optional => 0 },
		output_file	=> { is => 'Text', doc => "Output file for adjusted-coverage file", is_optional => 0 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Adjusts the coverage calculation for 454 efficiency pool analysis"                 
}

sub help_synopsis {
    return <<EOS
This command adjusts the coverage calculation for 454 efficiency pool analysis
EXAMPLE:	gmt analysis sammy
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

	my $coverage_file = $self->coverage_file;
	my $output_file = $self->output_file;

	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	print OUTFILE "chrom\tchr_start\tchr_stop\ttarget_name\tadjusted_depth\n";

	my $input = new FileHandle ($coverage_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		if($lineCounter == 1)
		{
			# Skip header 
		}
		else
		{
			(my $chrom, my $chr_start, my $chr_stop, my $target_name, my $target_size, my $covered_1x, my $covered_2x, my $covered_10x, my $covered_15x, my $covered_20x, my $avg_depth) = split(/\t/, $line);

			if($avg_depth < 0)
			{
				die "Average depth was less than zero at $line\n";
			}

			my $pct_1x = $covered_1x;
			$pct_1x =~ s/\%//;
			$pct_1x = $pct_1x / 100;
			
			my $adjusted_depth = $avg_depth * $pct_1x;
			$adjusted_depth = sprintf("%.1f", $adjusted_depth);
			print OUTFILE "$chrom\t$chr_start\t$chr_stop\t$target_name\t$adjusted_depth\n";
			
		}

	}
	
	close(OUTFILE);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}





sub commify
{
    local($_) = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}

1;

