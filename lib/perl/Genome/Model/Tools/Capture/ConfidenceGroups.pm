
package Genome::Model::Tools::Capture::ConfidenceGroups;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ConfidenceGroups - Build Genome Models for Capture Datasets
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

## Declare global statistics hash ##

my %stats = ();


class Genome::Model::Tools::Capture::ConfidenceGroups {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		varscan_file	=> { is => 'Text', doc => "File of variants in Varscan format", is_optional => 0, is_input => 1 },
		glf_file	=> { is => 'Text', doc => "File of variants in glfSomatic format", is_optional => 0, is_input => 1 },
		variant_file	=> { is => 'Text', doc => "File of Tiered variants to classify", is_optional => 0, is_input => 1 },
		output_high	=> { is => 'Text', doc => "Output file for high confidence" , is_optional => 0, is_input => 1, is_output => 1},
		output_highest	=> { is => 'Text', doc => "Output file for highest confidence" , is_optional => 0, is_input => 1, is_output => 1},
	],
	
	has_param => [
		lsf_resource => { default_value => 'select[model!=Opteron250 && type==LINUX64 && mem>6000] rusage[mem=6000]'},
       ],	
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Classifies variant calls as high or highest confidence"                 
}

sub help_synopsis {
    return <<EOS
Classifies variant calls as high or "highest confidence based on detection by 2 algorithms
EXAMPLE:	gmt capture merge-variant-calls ...
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

	my $varscan_file = $self->varscan_file;
	my $glf_file = $self->glf_file;
	my $variant_file = $self->variant_file;
	my $output_high = $self->output_high;
	my $output_highest = $self->output_highest;
	
	if(!(-e $varscan_file && -e $glf_file))
	{
		die "One or more files didn't exist!\n";
	}

	$stats{'highest-conf'} = $stats{'high-conf'} = $stats{'varscan-only'} = $stats{'sniper-only'} = 0;

	my %varscan = load_positions($varscan_file);
	my %sniper = load_positions($glf_file);

	open(HIGH, ">$output_high") or die "Can't open outfile: $!\n";
	open(HIGHEST, ">$output_highest") or die "Can't open outfile: $!\n";

	## Parse the variants file ##

	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $chrom, my $position) = split(/\t/, $line);
		my $key = "$chrom\t$position";

		if($varscan{$key} && $sniper{$key})
		{
			$stats{'highest-conf'}++;
			print HIGHEST "$line\n";
		}
		elsif($sniper{$key})
		{
			$stats{'high-conf'}++;
			$stats{'sniper-only'}++;
			print HIGH "$line\n";
		}
		elsif($varscan{$key})
		{
			$stats{'high-conf'}++;
			$stats{'varscan-only'}++;
			print HIGH "$line\n";
		}
		else
		{
			warn "$key not found!\n";
		}
	
	}
	
	close($input);

	print "$lineCounter variants in $variant_file\n";
	print $stats{'highest-conf'} . " highest-confidence (both algorithms)\n";
	print $stats{'high-conf'} . " high confidence\n";
	print $stats{'varscan-only'} . " varscan-only\n";
	print $stats{'sniper-only'} . " sniper-only\n";

	close(HIGH);
	close(HIGHEST);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_positions
{
	my $filename = shift(@_);
	my %positions = ();	

	my $input = new FileHandle ($filename);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $chrom, my $position) = split(/\t/, $line);
		my $key = "$chrom\t$position";
		$positions{$key} = 1;
	
	}
	
	close($input);
	
	return(%positions);
}


	
1;

