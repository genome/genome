
package Genome::Model::Tools::Somatic::FilterVariantClusters;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# LimitToRoi - Converts a pileup file to a three column file of chromosome, position, and bases with q>20.
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

class Genome::Model::Tools::Somatic::FilterVariantClusters {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variant_file	=> { is => 'Text', doc => "Input file of variants or positions", is_optional => 0, is_input => 1 },
		window_size	=> { is => 'Text', doc => "The window size to count variants", is_optional => 0, is_input => 1, default => 40000 },
		max_variants	=> { is => 'Text', doc => "Maximum allowed variants in window size", is_optional => 0, is_input => 1, default => 4},
		output_file     => { is => 'Text', doc => "Output file to receive filter-passed variants", is_optional => 0, is_input => 1, is_output => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Filter clusters of variants within a specified window size"                 
}

sub help_synopsis {
    return <<EOS
This command filters clusters of variants within a specified window size
EXAMPLE:	gmt capture filter-variant-clusters --variant-file [my.variants] --output-file [my.variants.roi]
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
	my $variant_file = $self->variant_file;
	my $window_size = $self->window_size;
	my $max_variants = $self->max_variants;
	my $output_file = $self->output_file;
	
	my %variants = ();
	
	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;
	my $totalVariants = my $failedVariants = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $position) = split(/\t/, $line);
		my $key = join("\t", $chrom, $position);
		$variants{$key} = $line;
		$totalVariants++;
	}
	
	close($input);
	
	## Define windowing parameters ##
	
	my $window_chrom = my $window_start = my $window_stop = my $window_variants = my $window_variant_list = "";
	my @failedWindows = ();
	my $numFailedWindows = 0;
	my $numWindows = 0;

	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	open(FAILEDFILE, ">$output_file.removed") or die "Can't open outfile: $!\n";
	
	## Go through variants ##
	
	foreach my $variant_key (sort byChrPos keys %variants)
	{
		my ($chrom, $position) = split(/\t/, $variant_key);
		my $variant = $variants{$variant_key};

		if(!$window_chrom)
		{
			## Start new window ##
			$window_chrom = $chrom;
			$window_start = $window_stop = $position;
			$window_variants = 1;
			$window_variant_list = $variant;
		}
		elsif($chrom ne $window_chrom)
		{
			## Process Current Window ##
			my $result = process_window($self, $window_chrom, $window_start, $window_stop, $window_variants, $window_variant_list);
			$numWindows++;
			if($result)
			{
				$failedWindows[$numFailedWindows] = $result;
				$numFailedWindows++;
				$failedVariants += $window_variants;
				print FAILEDFILE "$window_variant_list\n";
			}
			else
			{
				print OUTFILE "$window_variant_list\n";
			}

			## Start new window ##
			$window_chrom = $chrom;
			$window_start = $window_stop = $position;
			$window_variants = 1;
			$window_variant_list = $variant;
		}
		elsif(($position - $window_start) > $window_size)
		{
			## Process Current Window ##
			my $result = process_window($self, $window_chrom, $window_start, $window_stop, $window_variants, $window_variant_list);
			$numWindows++;
			if($result)
			{
				$failedWindows[$numFailedWindows] = $result;
				$numFailedWindows++;
				$failedVariants += $window_variants;
				print FAILEDFILE "$window_variant_list\n";
			}
			else
			{
				print OUTFILE "$window_variant_list\n";
			}
			## Start new window ##
			$window_chrom = $chrom;
			$window_start = $window_stop = $position;
			$window_variants = 1;
			$window_variant_list = $variant;
		}
		else
		{
			## Add variant to window ##
			$window_variants++;
			$window_stop = $position;
			$window_variant_list .= "\n" . $variant;
		}
	}

	## Process last window ##
	my $result = process_window($self, $window_chrom, $window_start, $window_stop, $window_variants, $window_variant_list);
	$numWindows++;
	if($result)
	{
		$failedWindows[$numFailedWindows] = $result;
		$numFailedWindows++;
		$failedVariants += $window_variants;
		print FAILEDFILE "$window_variant_list\n";
	}
	else
	{
		print OUTFILE "$window_variant_list\n";
	}

	print "$totalVariants total variants\n";
	print "$numWindows windows assessed\n";
	print "$numFailedWindows failed cluster-filter\n";
	print "$failedVariants variants removed\n";
	
	close(OUTFILE);
	close(FAILEDFILE);
	
	return 1;
}


################################################################################################
# Process_window - Check window against specified rules
#
################################################################################################

sub process_window
{                               # replace with real execution logic.
	my ($self, $window_chrom, $window_start, $window_stop, $window_variants, $window_variant_list) = @_;
	if($window_variants > $self->max_variants)
	{
		print join("\t", $window_chrom, $window_start, $window_stop, $window_variants) . "\n";
		print $window_variant_list . "\n\n";
	}

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

