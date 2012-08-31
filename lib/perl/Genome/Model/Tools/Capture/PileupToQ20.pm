
package Genome::Model::Tools::Capture::PileupToQ20;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# PileupToQ20 - Converts a pileup file to a three column file of chromosome, position, and bases with q>20.
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	02/10/2010 by D.K.
#	MODIFIED:	02/10/2010 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Capture::PileupToQ20 {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		pileup_file	=> { is => 'Text', doc => "Input file in pileup format", is_optional => 0, is_input => 1 },
		min_base_qual	=> { is => 'Text', doc => "Minimum base quality [20]", is_optional => 1, is_input => 1 },
		output_file     => { is => 'Text', doc => "Output file to receive per-base qual>min coverage", is_optional => 0, is_input => 1, is_output => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Calculates positional coverage with minimum base quality"                 
}

sub help_synopsis {
    return <<EOS
This command calculates positional coverage with minimum base quality
EXAMPLE:	gmt capture pileup-to-q20 --pileup-file [my.pileup] --output-file [my.pileup.q20]
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
	my $pileup_file = $self->pileup_file;
	my $output_file = $self->output_file;
	my $min_base_qual = 20;
	$min_base_qual = $self->min_base_qual if(defined($self->min_base_qual));
	
	## Open outfile ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	
	my %stats = ();
	$stats{'num_positions'} = $stats{'num_bases'} = $stats{'num_bases_qual'} = 0;
	
	## Parse the indels ##

	my $input = new FileHandle ($pileup_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my @lineContents = split(/\t/, $line);

		$stats{'num_positions'}++;

		my $chrom = $lineContents[0];
		my $position = $lineContents[1];
		my $ref_base = $lineContents[2];
		my $depth = $lineContents[3];
		my $qualities = $lineContents[5];

		## Go through each quality ##
		
		my @qualities = split(//, $qualities);
		my $num_quals = 0;
		my $qual_coverage = 0;
		
		foreach my $code (@qualities)
		{
			my $qual_score = ord($code) - 33;
			$num_quals++;

			$stats{'num_bases'}++;

			if($qual_score >= $min_base_qual)
			{
				$qual_coverage++;
				$stats{'num_bases_qual'}++;
			}
		}
		
		print OUTFILE "$chrom\t$position\t$qual_coverage\n";			
	}

	close($input);

	print "$stats{'num_positions'} positions in pileup file\n";	
	print "$stats{'num_bases'} total read bases\n";
	print "$stats{'num_bases_qual'} had base quality >= $min_base_qual\n";
	
	close(OUTFILE);
}


1;

