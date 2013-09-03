package Genome::Model::Tools::Capture::CombineBase;

use strict;
use warnings;
use FileHandle;
use Genome;
use Genome::Model::Tools::Capture::Helpers qw/fix_chrom/;

class Genome::Model::Tools::Capture::CombineBase {
	is => 'Genome::Model::Tools::Capture',
};

################################################################################################
# Load variants - parse a variant file 
#
################################################################################################

sub load_variants
{                               # replace with real execution logic.
	my $FileName = shift(@_);	

	my %variants = ();

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	my @formatted = ();
	my $formatCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my @lineContents = split(/\t/, $line);

		if(!(lc($lineContents[0]) =~ "chrom" || lc($lineContents[0]) =~ "ref_name"))
		{
			my $chrom = $lineContents[0];
			$chrom = fix_chrom($chrom);
			my $chr_start = $lineContents[1];
			my $chr_stop = $lineContents[2];
			my $allele1= $lineContents[3];
			my $allele2 = $lineContents[4];
			my $numContents = @lineContents;
			
			$allele1 = "-" if($allele1 eq "0");
			$allele2 = "-" if($allele2 eq "0");
			
			my $rest_of_line = "";
			for(my $colCounter = 5; $colCounter < $numContents; $colCounter++)
			{
				$rest_of_line .= "\t" if($rest_of_line);
				$rest_of_line .= $lineContents[$colCounter];
			}

			if($chrom && $chr_start && $chr_stop)
			{
				my $key = join("\t", $chrom, $chr_start, $chr_stop);
				$variants{$key} = "$allele1\t$allele2";
				$variants{$key} .= "\t$rest_of_line" if($rest_of_line);
			}
		}
	}

	close($input);
	
	return(%variants);
	
}
