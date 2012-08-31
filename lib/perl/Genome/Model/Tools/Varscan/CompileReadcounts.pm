
package Genome::Model::Tools::Varscan::CompileReadcounts;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Varscan::Somatic	Runs Varscan somatic pipeline on Normal/Tumor BAM files
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	12/09/2009 by D.K.
#	MODIFIED:	12/29/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Varscan::CompileReadcounts {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		readcounts	=> { is => 'Text', doc => "A Varscan readcounts file or csv-separated list of one", is_optional => 0 },
		variants_file	=> { is => 'Text', doc => "Variants in annotation format", is_optional => 0 },
		output_file	=> { is => 'Text', doc => "Path to output file" , is_optional => 0},
		other_alleles	=> { is => 'Text', doc => "If set to 1, output non-expected alleles" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Compiles readcounts across multiple samples"                 
}

sub help_synopsis {
    return <<EOS
Runs Varscan readcounts from BAM files
EXAMPLE:	gmt varscan compile-readcounts --variants myVariants --readcounts [fof]
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

sub execute
{	
	my $self = shift;

	## Get required parameters ##
	my $variants_file = $self->variants_file;
	my $readcounts_file = $self->readcounts;
	my $reference = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa';
	my $output_file = $self->output_file;

	if($self->output_file)
	{
		open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
		print OUTFILE "variant";		
	}

	## Load the variants ##
	
	my %target_variants = load_variants($variants_file);


	my @readcounts_files = split(/\,/, $readcounts_file);

	my %results_by_filename = ();
	
	foreach my $filename (@readcounts_files)
	{
		print OUTFILE "\t$filename";
		my $input = new FileHandle ($filename);
		my $lineCounter = 0;
		
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;

			my ($chrom, $position, $ref, $depth, $quality_depth) = split(/\t/, $line);			

			my $position_key = "$chrom\t$position";

			my @lineContents = split(/\t/, $line);
			my $numContents = @lineContents;				
			
			## Proceed if targeted variant ##
			
			if($target_variants{$position_key} && $numContents > 4)
			{
				my $reads1 = my $reads2 = 0;
				my $var_freq = "0.00%";
				
				my @targetContents = split(/\t/, $target_variants{$position_key});
				my $variant_allele = $targetContents[4];
				

				my $ref_counts = $lineContents[5];
				my @refCounts = split(/\:/, $ref_counts);
				$reads1 = $refCounts[1];

				my $var_counts = "";
				my $other_alleles = "";

				for(my $colCounter = 6; $colCounter < $numContents; $colCounter++)
				{
					(my $allele) = split(/\:/, $lineContents[$colCounter]);
					if($allele eq $variant_allele)
					{
						$var_counts = $lineContents[$colCounter];
						my @varCounts = split(/\:/, $var_counts);
						$reads2 = $varCounts[1];
						$var_freq = sprintf("%.2f", ($reads2 / ($reads1 + $reads2) * 100)) . '%';
					}
					elsif($self->other_alleles)
					{
						my @varCounts = split(/\:/, $var_counts);
						$reads2 = $varCounts[1];
#						$var_freq = sprintf("%.2f", ($reads2 / ($reads1 + $reads2) * 100)) . '%';
						$other_alleles .= "|" if($other_alleles);
						$other_alleles .= "$allele" . ":" . "$reads2";
					}
				}

				my $result_key = "$position_key\t$filename";
				$results_by_filename{$result_key} = "$reads1\t$reads2\t$var_freq";
				$results_by_filename{$result_key} .= "\t$other_alleles" if($self->other_alleles);
			}


		}
		
		close($input);		
	}
	
	print OUTFILE "\n";
	
	## Once all filenames have been handled, print the variant results ##
	
	foreach my $position_key (sort byChrPos keys %target_variants)
	{
		print OUTFILE $target_variants{$position_key};
		
		foreach my $filename (@readcounts_files)
		{
			my $file_result = "-\t-\t-";
			my $result_key = "$position_key\t$filename";
			$file_result = $results_by_filename{$result_key} if($results_by_filename{$result_key});
			print OUTFILE "\t" . $file_result;
		}
		
		print OUTFILE "\n";
	}


	close(OUTFILE);

}



sub byChrPos
{
	my ($chrom_a, $pos_a) = split(/\t/, $a);
	my ($chrom_b, $pos_b) = split(/\t/, $b);
	
	$chrom_a = "23" if($chrom_a eq "X");
	$chrom_a = "24" if($chrom_a eq "Y");
	$chrom_a = "25" if(substr($chrom_a, 0, 1) eq "M");
	$chrom_a =~ s/[^0-9]//g;

	$chrom_b = "23" if($chrom_b eq "X");
	$chrom_b = "24" if($chrom_b eq "Y");
	$chrom_b = "25" if(substr($chrom_b, 0, 1) eq "M");
	$chrom_b =~ s/[^0-9]//g;

	$chrom_a <=> $chrom_b
	or
	$pos_a <=> $pos_b;
}



################################################################################################
# Load variants
#
################################################################################################

sub load_variants
{
	my $variants_file = shift(@_);
	my %variants = ();

	## Load readcounts ##
	
	my $input = new FileHandle ($variants_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		my ($chrom, $start, $stop, $ref, $var) = split(/\t/, $line);
		my $position_key = "$chrom\t$start";
		
		$variants{$position_key} = $line;
	}
	
	close($input);

	return(%variants);	
}


################################################################################################
# Load variants
#
################################################################################################

sub load_readcounts
{	

}








1;



