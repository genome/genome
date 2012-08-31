
package Genome::Model::Tools::Capture::MatchIndelsToDbsnp;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MatchIndelsToDbsnpForAnnotation - Merge glfSomatic/Varscan somatic calls in a file that can be converted to MAF format
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

class Genome::Model::Tools::Capture::MatchIndelsToDbsnp {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variant_file	=> { is => 'Text', doc => "File of indels in annotation format", is_optional => 0, is_input => 1 },
		dbsnp_file	=> { is => 'Text', doc => "UCSC file of dbSNP indels (pre-limited if possible)", is_optional => 0, is_input => 1 },		
		output_file     => { is => 'Text', doc => "Output file of matched indels", is_optional => 0, is_input => 1, is_output => 1 },
		not_file     => { is => 'Text', doc => "Output file of indels that did NOT match dbSNP", is_optional => 1, is_input => 1, is_output => 1 },

	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Match a set of small indels in annotation format to UCSC dbSNP indel file"                 
}

sub help_synopsis {
    return <<EOS
This command matches a set of small indels in annotation format to UCSC dbSNP indel file
EXAMPLE:	gmt capture match-indels-to-dbsnp --variant-file [file] --dbsnp-file [file] --output-file [file]
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
	my $dbsnp_file = $self->dbsnp_file;
	my $output_file = $self->output_file;
	
	my %dbsnp = load_dbsnp($dbsnp_file);
	
	## Open outfile ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	open(INFOFILE, ">$output_file.info") or die "Can't open outfile: $!\n";
	open(NOTFILE, ">" . $self->not_file) or die "Can't open outfile: $!\n" if($self->not_file);
	
	## Parse the indels ##

	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;

	my @formatted = ();
	my $formatCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);

		if(lc($chrom) =~ "chrom")
		{
			## Skip header ##
		}
		else
		{
			my $indel_type = "";
			my $indel_size = 0;
			if($ref eq '-' || $ref eq '0' || $ref eq '*' || length($var) > length($ref))
			{
				## INSERTION ##
				$indel_type = "INS";
				$indel_size = length($var);
			}
			else
			{
				$indel_type = "DEL";
				$indel_size = length($ref);
			}

			my $matches = "";
			
			for(my $position = $chr_start - $indel_size - 1; $position <= $chr_stop + $indel_size + 1; $position++)
			{
				my $key = join("\t", $chrom, $position);
				if($dbsnp{$key})
				{
					my ($dbsnp_chrom, $dbsnp_chr_start, $dbsnp_chr_stop, $dbsnp_ref, $dbsnp_var) = split(/\t/, $dbsnp{$key});
					my $dbsnp_type = "";
					my $dbsnp_size = 0;
					
					if($dbsnp_ref eq "-" || length($dbsnp_var) > length($dbsnp_ref))
					{
						$dbsnp_type = "INS";
						$dbsnp_size = length($var);
						$dbsnp_size -= length($ref) if($ref ne "-");
					}
					else
					{
						$dbsnp_type = "DEL";
						$dbsnp_size = length($ref);
						$dbsnp_size -= length($var) if($var ne "-");
					}
					
					if($dbsnp_type eq $indel_type)
					{
						my $size_diff = abs($indel_size - $dbsnp_size);
						my $size_diff_pct = $size_diff / $indel_size;
						if($size_diff <= 2 || $size_diff_pct <= 0.20)
						{
							$matches .= "\n" if($matches);
							$matches .= $dbsnp{$key};							
						}
					}

				}
			}
			
			if($matches)
			{
				print OUTFILE "$line\n";
				my $best_match = get_best_match($chr_start, $chr_stop, $ref, $var, $indel_type, $indel_size, $matches);
				print INFOFILE "$best_match\n";
			}
			elsif($self->not_file)
			{
				print NOTFILE "$line\n";
			}
		}
	}

	close($input);
	
	
	close(OUTFILE);
	close(INFOFILE);
}


#############################################################
# load dbsnp variants
#
#############################################################

sub get_best_match
{
	my ($chr_start, $chr_stop, $ref, $var, $indel_type, $indel_size, $matches) = @_;
	
	my @matches = split(/\n/, $matches);
	my $num_matches = @matches;
	return($matches) if($num_matches == 1);
	
	my @unsorted = ();
	my $num_unsorted = 0;
	foreach my $match (@matches)
	{
		my ($dbsnp_chrom, $dbsnp_chr_start, $dbsnp_chr_stop, $dbsnp_ref, $dbsnp_var) = split(/\t/, $match);
		my $dbsnp_type = "";
		my $dbsnp_size = 0;
		
		if($dbsnp_ref eq "-" || length($dbsnp_var) > length($dbsnp_ref))
		{
			$dbsnp_type = "INS";
			$dbsnp_size = length($var);
			$dbsnp_size -= length($ref) if($ref ne "-");
		}
		else
		{
			$dbsnp_type = "DEL";
			$dbsnp_size = length($ref);
			$dbsnp_size -= length($var) if($var ne "-");
		}
		
		if($dbsnp_type eq $indel_type)
		{
			my $size_diff = abs($indel_size - $dbsnp_size);
			my $position_diff = abs($chr_start - $dbsnp_chr_start);
			$unsorted[$num_unsorted] = join("\t", $size_diff, $position_diff, $match);
			$num_unsorted++;
		}	
	}
	
	my @sorted = sort bySizeDistance (@unsorted);
	
	my @bestContents = split(/\t/, $sorted[0]);
	my $best_match = "";
	my $numContents = @bestContents;
	for(my $colCounter = 2; $colCounter < $numContents; $colCounter++)
	{
		$best_match .= "\t" if($best_match);
		$best_match .= $bestContents[$colCounter];
	}

#	print join("\t", $chr_start, $chr_stop, $ref, $var, $indel_type, $indel_size) . "\n";
#	print "$matches\nChose: $best_match\n";

	return($best_match);
	
}



sub bySizeDistance
{
	my ($size_a, $dist_a) = split(/\t/, $a);
	my ($size_b, $dist_b) = split(/\t/, $b);
	$size_a <=> $size_b
	or
	$dist_a <=> $dist_b;
}

#############################################################
# load dbsnp variants
#
#############################################################

sub load_dbsnp
{
	my $FileName = shift(@_);
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	my %indels = ();

	print "Loading dbSNP...\n";

	my @formatted = ();
	my $formatCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my ($chrom, $chr_start, $chr_stop) = split(/\t/, $line);
		my $key = join("\t", $chrom, $chr_start);
		
		if($indels{$key})
		{
			$indels{$key} .= "\n" . $line;
		}
		else
		{
			$indels{$key} = $line;
		}
		
	}
	
	close($input);
	
	print "$lineCounter indels loaded\n";
	
	return(%indels);
}


1;

