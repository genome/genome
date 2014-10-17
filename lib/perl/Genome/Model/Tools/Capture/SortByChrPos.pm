
package Genome::Model::Tools::Capture::SortByChrPos;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SortByChrPos - Build Genome Models for Capture Datasets
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

class Genome::Model::Tools::Capture::SortByChrPos {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		input_file	=> { is => 'Text', doc => "File of unsorted variants", is_optional => 0, is_input => 1 },
		output_file	=> { is => 'Text', doc => "Output file to contain sorted results" , is_optional => 0, is_input => 1, is_output => 1},
	],
	
	has_param => [
		lsf_resource => { default_value => 'select[mem>6000] rusage[mem=6000]'},
       ],	
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Sort variants by chromosome and position, e.g. to facilitate annotation"                 
}

sub help_synopsis {
    return <<EOS
Sorts a file by chromosome and position
EXAMPLE:	gmt capture sort-by-chr-pos ...
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

	my $input_file = $self->input_file;
	my $output_file = $self->output_file;
	
	if(!(-e $input_file))
	{
		die "Input file $input_file not found!\n";
	}

	my $input = new FileHandle ($input_file);
	my @lines = ();
	my $lineCounter = 0;
	
	my $header = "";
	
	while (<$input>)
	{
		chomp;
		my $line = $_;

		if($lineCounter == 0 && $line && substr($line, 0, 5) eq "chrom")
		{
			$header = $line;
		}
		else
		{
			$lines[$lineCounter] = $line;
			$lineCounter++;					
		}

	}

	close($input);

	## Sort the entries ##
	warn "Sorting variants by chromosome and position...\n";
	
	@lines = sort byChrPos (@lines);	
	
	## Open outfile ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	print OUTFILE "$header\n" if($header);
	
#	foreach my $line (@lines)
	for(my $printCounter = 0; $printCounter < $lineCounter; $printCounter++)
	{
		my $line = $lines[$printCounter];
		print OUTFILE "$line\n";
	}
	
	close(OUTFILE);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



sub byChrPos
{
	(my $chrom_a, my $pos_a) = split(/\t/, $a);
	(my $chrom_b, my $pos_b) = split(/\t/, $b);
	
	## Remove lowercase chr ##
	
	$chrom_a =~ s/chr//;
	$chrom_b =~ s/chr//;
	
	## Replace characters with their numerical equivalents ##
	
	if($chrom_a eq 'X')
	{
		$chrom_a = 23;
	}
	elsif($chrom_a eq 'Y')
	{
		$chrom_a = 24;
	}
	elsif($chrom_a =~ 'M')
	{
		$chrom_a = 25;
	}
	elsif($chrom_a =~ 'NT')
	{
		$chrom_a =~ s/NT\_//;
	}


	if($chrom_b eq 'X')
	{
		$chrom_b = 23;
	}
	elsif($chrom_b eq 'Y')
	{
		$chrom_b = 24;
	}
	elsif($chrom_b =~ 'M')
	{
		$chrom_b = 25;
	}
	elsif($chrom_b =~ 'NT')
	{
		$chrom_b =~ s/NT\_//;
	}

	## Make sure there's nothing but numbers ##
	
	$chrom_a =~ s/[^0-9]//g;
	$chrom_b =~ s/[^0-9]//g;
	$pos_a =~ s/[^0-9]//g;
	$pos_b =~ s/[^0-9]//g;


	$chrom_a <=> $chrom_b
	or
	$pos_a <=> $pos_b;

	
}





1;

