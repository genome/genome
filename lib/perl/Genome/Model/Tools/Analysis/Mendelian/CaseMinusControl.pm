
package Genome::Model::Tools::Analysis::Mendelian::CaseMinusControl;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SearchRuns - Search the database for runs
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	04/01/2009 by D.K.
#	MODIFIED:	04/01/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::Mendelian::CaseMinusControl {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		case_files	=> { is => 'Text', doc => "Comma-separated path to case variant files", is_optional => 0, is_input => 1},
		control_files	=> { is => 'Text', doc => "Comma-separated path to control variant files", is_optional => 0, is_input => 1},
		output_file	=> { is => 'Text', doc => "Output file for passed results", is_optional => 1, is_input => 1, is_output => 1},
		removed_file	=> { is => 'Text', doc => "A file to contain variants that were removed", is_optional => 1, is_input => 1, is_output => 1},
		min_cases	=> { is => 'Text', doc => "The minimum number of cases to include a variant", is_input => 1, default => 1},		
		max_controls	=> { is => 'Text', doc => "The maximum number of controls to include a variant", is_input => 1, default => 0},		
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Runs a simplistic Mendelian analysis of variants: found in X cases but no more than Y controls"                 
}

sub help_synopsis {
    return <<EOS
This command runs a simplistic Mendelian analysis of variants: found in X cases but no more than Y controls
EXAMPLE:	gmt analysis mendelian case-minus-control --case-files affected1,affected2 --control-files control1
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

	my $output_file = $self->output_file;
	my $removed_file = $self->removed_file;
	
	my %stats = ();

	## Reset master hashes to save inforamtion ##
	
	my %variant_annotation = ();
	my %variant_cases = my %variant_controls = ();

	## Go through each case file ##
	
	my @case_files = split(/\,/, $self->case_files);
	$stats{'num_cases'} = @case_files;
	
	foreach my $case_file (@case_files)
	{
		## load the variant ##
		my %variants = load_variants($case_file);
		
		## Record variant annotation and count this case ##
		foreach my $variant (keys %variants)
		{
			$variant_annotation{$variant} = $variants{$variant};
			
			if($variant_cases{$variant})
			{
				$variant_cases{$variant}++;
			}
			else
			{
				$variant_cases{$variant} = 1;
			}
			
		}
	}


	## Go through each case file ##
	
	my @control_files = split(/\,/, $self->control_files);
	$stats{'num_controls'} = @control_files;
	
	foreach my $control_file (@control_files)
	{
		## load the variant ##
		my %variants = load_variants($control_file);
		
		## Record variant annotation and count this case ##
		foreach my $variant (keys %variants)
		{
			$variant_annotation{$variant} = $variants{$variant};
			
			if($variant_controls{$variant})
			{
				$variant_controls{$variant}++;
			}
			else
			{
				$variant_controls{$variant} = 1;
			}
			
		}
	}


	open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
	open(REMOVED, ">" . $self->removed_file) or die "Can't open removed file: $!\n";

	print OUTFILE "variant\tnum_cases\tnum_controls\n";
	print REMOVED "variant\tnum_cases\tnum_controls\treason\n";

	## Go through all variants ##
	
	foreach my $variant (sort byChrPos keys %variant_annotation)
	{
		$stats{'variants'}++;
		my $num_cases = my $num_controls = 0;
		
		$num_cases = $variant_cases{$variant} if($variant_cases{$variant});
		$num_controls = $variant_controls{$variant} if($variant_controls{$variant});
		
		if($num_cases >= $self->min_cases)
		{
			$stats{'variants_pass_cases'}++;
			if($num_controls <= $self->max_controls)
			{
				$stats{'variants_pass_controls'}++;
				print OUTFILE join("\t", $variant_annotation{$variant}, $num_cases, $num_controls) . "\n";
			}
			else
			{
#				$stats{'variants_too_many_controls'}++;
				print REMOVED join("\t", $variant_annotation{$variant}, $num_cases, $num_controls, "TooManyControls") . "\n";
			}
		}
		else
		{
#			$stats{'variants_too_few_cases'}++;
			print REMOVED join("\t", $variant_annotation{$variant}, $num_cases, $num_controls, "TooFewCases") . "\n";
		}
	}
	
	close(OUTFILE);
	close(REMOVED);

	
	
	## Print stats ##
	
	foreach my $key (sort keys %stats)
	{
		print "$stats{$key}\t$key\n";
	}


	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# Load Variants
################################################################################################

sub load_variants
{
	my $FileName = shift(@_);
	my %variants = ();
	## Print the variants ##

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);
		my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
		
		$variants{$key} = $line;
	}
	
	close($input);	

	return(%variants);
}


################################################################################################
# sorting subroutine by chromosome and then position
#
################################################################################################

sub byChrPos
{
	my ($chrom_a, $pos_a) = split(/\t/, $a);
	my ($chrom_b, $pos_b) = split(/\t/, $b);

	$chrom_a = 23 if($chrom_a eq 'X');
	$chrom_a = 24 if($chrom_a eq 'Y');
	$chrom_a = 25 if($chrom_a eq 'MT');
	
	$chrom_b = 23 if($chrom_b eq 'X');
	$chrom_b = 24 if($chrom_b eq 'Y');
	$chrom_b = 25 if($chrom_b eq 'MT');

	$chrom_a =~ s/[^0-9]//g;
	$chrom_b =~ s/[^0-9]//g;
	
	$chrom_a <=> $chrom_b
	or
	$pos_a <=> $pos_b;
}


1;


