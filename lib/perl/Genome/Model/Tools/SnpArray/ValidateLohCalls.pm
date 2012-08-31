
package Genome::Model::Tools::SnpArray::ValidateLohCalls;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::SnpArray::ValidateLohCalls {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		tumor_genotype_file	=> { is => 'Text', doc => "Three-column file of genotype calls chrom, pos, genotype", is_optional => 0, is_input => 1 },
		normal_genotype_file	=> { is => 'Text', doc => "Three-column file of genotype calls chrom, pos, genotype", is_optional => 0, is_input => 1 },
		variant_file	=> { is => 'Text', doc => "Varscan LOH Calls in annotation format", is_optional => 0, is_input => 1 },
		sample_name	=> { is => 'Text', doc => "Name of the sample for output file", is_optional => 1, is_input => 1 },
		min_depth	=> { is => 'Text', doc => "Minimum depth to compare a het call [8]", is_optional => 1, is_input => 1, default => 8},		
		verbose	=> { is => 'Text', doc => "Turns on verbose output [0]", is_optional => 1, is_input => 1},
		flip_alleles 	=> { is => 'Text', doc => "If set to 1, try to avoid strand issues by flipping alleles to match", is_optional => 1, is_input => 1},
		output_file	=> { is => 'Text', doc => "Output file for QC result", is_optional => 1, is_input => 1}
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Validates Varscan LOH calls using SNP array data"                 
}

sub help_synopsis {
    return <<EOS
This command validates Varscan LOH calls using SNP array data
EXAMPLE:	gmt snp-array validate-germline-calls --genotype-file affy.genotypes --variant-file lane1.var
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
	my $tumor_genotype_file = $self->tumor_genotype_file;
	my $normal_genotype_file = $self->normal_genotype_file;
	my $min_depth = $self->min_depth;
	my $output_file = $self->output_file;

	my %stats = ();
	$stats{'num_loh_calls'} = $stats{'had_array_data'} = $stats{'met_min_depth'} = $stats{'array_was_germline'} = $stats{'array_was_LOH'} = $stats{'array_was_GOH'} = $stats{'array_was_mismatch'} = 0;


	print "Loading tumor genotypes...\n";
	my %tumor_genotypes = load_genotypes($tumor_genotype_file);

	print "Loading normal genotypes...\n";
	my %normal_genotypes = load_genotypes($normal_genotype_file);


	print "Parsing variant calls in $variant_file...\n" if($self->verbose);

	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;

	my $verbose_output = "";

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my @lineContents = split(/\t/, $line);
		my $chrom = $lineContents[0];
		my $position = $lineContents[1];
		my $ref_base = $lineContents[3];
		my $normal_reads1 = $lineContents[5];
		my $normal_reads2 = $lineContents[6];
		my $normal_call = $lineContents[8];

		my $tumor_reads1 = $lineContents[9];
		my $tumor_reads2 = $lineContents[10];
		my $tumor_call = $lineContents[12];

		my $normal_coverage = $normal_reads1 + $normal_reads2;
		my $tumor_coverage = $tumor_reads1 + $tumor_reads2;
		
		my $snp_key = $chrom . "\t" . $position;

		$stats{'num_loh_calls'}++;

		if($normal_genotypes{$snp_key} && $tumor_genotypes{$snp_key})
		{
			$stats{'had_array_data'}++;
			
			if($normal_coverage >= $min_depth && $tumor_coverage >= $min_depth)
			{
				$stats{'met_min_depth'}++;
				my $normal_gt = $normal_genotypes{$snp_key};
				my $tumor_gt = $tumor_genotypes{$snp_key};
				
				if($normal_gt eq $tumor_gt)
				{
					$stats{'array_was_germline'}++;
				}
				elsif(is_heterozygous($normal_gt) && is_homozygous($tumor_gt))
				{
					$stats{'array_was_LOH'}++;
				}
				elsif(is_homozygous($normal_gt) && is_heterozygous($tumor_gt))
				{
					$stats{'array_was_GOH'}++;
				}
				else
				{
					$stats{'array_was_mismatch'}++;
				}							
			}

		}

		
	}
	
	close($input);


	## Calculate the LOH concordance ##
	
	my $loh_concordance = $stats{'array_was_LOH'} / $stats{'met_min_depth'} * 100;
	$loh_concordance = sprintf("%.2f", $loh_concordance) . '%';


	print $stats{'num_loh_calls'} . " Varscan LOH calls\n";
	print $stats{'had_array_data'} . " had SNP array data\n";
	print $stats{'met_min_depth'} . " met min depth\n";
	print $stats{'array_was_germline'} . " array was Germline\n";
	print $stats{'array_was_LOH'} . " array was LOH\n";
	print $stats{'array_was_GOH'} . " array was GOH\n";
	print $stats{'array_was_mismatch'} . " array was mismatch\n";
	print "$loh_concordance LOH concordance\n";

	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
		print OUTFILE "Sample\tLOHcalls\tWithArrayData\tMetMinDepth\tArrayGermline\tArrayGOH\tArrayMismatch\tArrayLOH\tConcordance\n";
		print OUTFILE join("\t", $variant_file, $stats{'num_loh_calls'}, $stats{'had_array_data'}, $stats{'met_min_depth'}, $stats{'array_was_germline'}, $stats{'array_was_GOH'}, $stats{'array_was_mismatch'}, $stats{'array_was_LOH'}, $loh_concordance) . "\n";
		close(OUTFILE);
	}

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# Load Genotypes
#
################################################################################################

sub load_genotypes
{                               # replace with real execution logic.
	my $genotype_file = shift(@_);
	my %genotypes = ();
	
	my $input = new FileHandle ($genotype_file);
	my $lineCounter = 0;
	my $gtCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $chrom, my $position, my $genotype) = split(/\t/, $line);

		if($chrom && $position)
		{
			my $key = "$chrom\t$position";
			
			if($genotype && $genotype ne "--" && $genotype ne "NN")
			{
				$genotypes{$key} = $genotype;
				$gtCounter++;
			}			
		}

	}
	close($input);

	print "$gtCounter genotypes loaded\n";

#	print "$gtCounter genotypes loaded\n";
	
	return(%genotypes);                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

################################################################################################
# Sorting
#
################################################################################################


sub byBamOrder
{
	my ($chrom_a, $pos_a) = split(/\t/, $a);
	my ($chrom_b, $pos_b) = split(/\t/, $b);
	
	$chrom_a =~ s/X/9\.1/;
	$chrom_a =~ s/Y/9\.2/;
	$chrom_a =~ s/MT/25/;
	$chrom_a =~ s/M/25/;
	$chrom_a =~ s/NT/99/;
	$chrom_a =~ s/[^0-9\.]//g;

	$chrom_b =~ s/X/9\.1/;
	$chrom_b =~ s/Y/9\.2/;
	$chrom_b =~ s/MT/25/;
	$chrom_b =~ s/M/25/;
	$chrom_b =~ s/NT/99/;
	$chrom_b =~ s/[^0-9\.]//g;

	$chrom_a <=> $chrom_b
	or
	$pos_a <=> $pos_b;
}


################################################################################################
# Load Genotypes
#
################################################################################################

sub is_heterozygous
{
	my $gt = shift(@_);
	(my $a1, my $a2) = split(//, $gt);
	return(1) if($a1 ne $a2);
	return(0);
}



################################################################################################
# Load Genotypes
#
################################################################################################

sub is_homozygous
{
	my $gt = shift(@_);
	(my $a1, my $a2) = split(//, $gt);
	return(1) if($a1 eq $a2);
	return(0);
}



################################################################################################
# Load Genotypes
#
################################################################################################

sub flip_genotype
{
	my $gt = shift(@_);
	(my $a1, my $a2) = split(//, $gt);

	if($a1 eq "A")
	{
		$a1 = "T";
	}
	elsif($a1 eq "C")
	{
		$a1 = "G";
	}
	elsif($a1 eq "G")
	{
		$a1 = "C";
	}	
	elsif($a1 eq "T")
	{
		$a1 = "A";		
	}

	if($a2 eq "A")
	{
		$a2 = "T";
	}
	elsif($a2 eq "C")
	{
		$a2 = "G";
	}
	elsif($a2 eq "G")
	{
		$a2 = "C";
	}	
	elsif($a2 eq "T")
	{
		$a2 = "A";		
	}
	
	$gt = $a1 . $a2;
	$gt = sort_genotype($gt);
	return($gt);
}

################################################################################################
# Load Genotypes
#
################################################################################################

sub sort_genotype
{
	my $gt = shift(@_);
	(my $a1, my $a2) = split(//, $gt);

	my @unsorted = ($a1, $a2);
	my @sorted = sort @unsorted;
	$a1 = $sorted[0];
	$a2 = $sorted[1];
	return($a1 . $a2);
}



sub code_to_genotype
{
	my $code = shift(@_);
	
	return("AA") if($code eq "A");
	return("CC") if($code eq "C");
	return("GG") if($code eq "G");
	return("TT") if($code eq "T");

	return("AC") if($code eq "M");
	return("AG") if($code eq "R");
	return("AT") if($code eq "W");
	return("CG") if($code eq "S");
	return("CT") if($code eq "Y");
	return("GT") if($code eq "K");

#	warn "Unrecognized ambiguity code $code!\n";

	return("NN");	
}



sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;

