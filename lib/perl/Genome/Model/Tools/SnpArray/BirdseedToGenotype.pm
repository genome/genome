
package Genome::Model::Tools::SnpArray::BirdseedToGenotype;     # rename this when you give the module file a different name <--

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

my %stats = ();
my %snp_info = ();


class Genome::Model::Tools::SnpArray::BirdseedToGenotype {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		birdseed_file	=> { is => 'Text', doc => "Three-column file of genotype calls chrom, pos, genotype", is_optional => 0, is_input => 1 },
		output_file	=> { is => 'Text', doc => "Output file for QC result", is_optional => 1, is_input => 1},
		snp_info_file	=> { is => 'Text', doc => "SNP array information file downloaded from UCSC", is_optional => 0, is_input => 1, default=> "/gscmnt/xp4101/info/medseq/snp_array_data/snpArrayAffy6.txt" },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Convert birdseed files into chrom-pos-genotype files"                 
}

sub help_synopsis {
    return <<EOS
This command converts SNP array genotypes in birdseed format to a simpler chromosome-position-genotype format
EXAMPLE:	gmt snp-array birdseed-to-genotype --birdseed-file TCGAb49_SNP_1N_GenomeWideSNP_6_C07_653286.birdseed.data.txt --output-file converted_UCEC_Genome_Wide_SNP_6.Genotypes/H_LR-AP-A053-01A-21D-A013-09.genotype
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
	
	my $birdseed_file = $self->birdseed_file;
	my $output_file = $self->output_file;

	print "Loading the SNP information...\n";
	%snp_info = load_snp_info($self->snp_info_file);

	## Open the unsorted output file ##

	open(OUTFILE, ">$output_file.unsorted") or die "Can't open output file $output_file.unsorted : $!\n";
	print "Parsing birdseed file...\n";
	parse_birdseed_file($birdseed_file);

	close(OUTFILE);

	## Sort the output file ##
	
	my $cmd = "gmt capture sort-by-chr-pos --input $output_file.unsorted --output $output_file";
	system($cmd);

	## Remove unsorted file ##
	print "Removing unsorted file...\n";
	system("rm -rf $output_file.unsorted") if (-e "$output_file");

	print "Done\n";
}
	


#############################################################
# parse_birdseed_file - parses the birdseed genotype file
#
#############################################################

sub parse_birdseed_file
{
	my $FileName = shift(@_);
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		if($lineCounter > 2)
		{
			my ($snp_id, $call, $score) = split(/\t/, $line);
			
			if($snp_info{$snp_id})
			{
				my ($chrom, $position, $strand, $allele1, $allele2) = split(/\t/, $snp_info{$snp_id});
				
				my $genotype = "";
				
				if($call eq "0")
				{
					$genotype = $allele1 . $allele1;
				}
				elsif($call eq "1")
				{
					$genotype = $allele1 . $allele2;
				}
				elsif($call eq "2")
				{
					$genotype = $allele2 . $allele2;
				}
				else
				{
					$genotype = "NN";
				#	warn "Warning: unknown call value $call\n";	
				}
				
				if($genotype)
				{
					if($strand eq "-" && $genotype ne "NN")
					{
						$genotype = flip_genotype($genotype);
					}

					print OUTFILE join("\t", $chrom, $position, $genotype) . "\n";
				}

			}
		}
	}

	close($input);
}




#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub load_snp_info
{
	my $FileName = shift(@_);
	my %info = ();

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my @lineContents = split(/\t/, $line);
		my $chrom = $lineContents[0];
		my $position = $lineContents[2];
		my $snp_id = $lineContents[3];
		my $strand = $lineContents[5];
		my $alleles = $lineContents[6];
		(my $allele1, my $allele2) = split(/\//, $alleles);

		$chrom =~ s/chr//;

		$info{$snp_id} = "$chrom\t$position\t$strand\t$allele1\t$allele2";
	}

	close($input);

	return(%info);
}




#############################################################
# Flip genotype - reverse complement a genotype
#
#############################################################

sub flip_genotype
{
	my $genotype = shift(@_);
	
	my $allele1 = substr($genotype, 0, 1);
	my $allele2 = substr($genotype, 1, 1);

	$allele1 = flip_base($allele1);
	$allele2 = flip_base($allele2);
	
	$genotype = $allele1 . $allele2;
	
	return($genotype);
}


#############################################################
# Flip base - flip a base
#
#############################################################

sub flip_base
{
	my $base = shift(@_);
	
	return("A") if ($base eq "T");
	return("C") if ($base eq "G");
	return("G") if ($base eq "C");
	return("T") if ($base eq "A");
	return($base);
}




1;