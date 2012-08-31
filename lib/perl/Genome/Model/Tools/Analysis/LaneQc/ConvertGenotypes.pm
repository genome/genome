
package Genome::Model::Tools::Analysis::LaneQc::ConvertGenotypes;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::Analysis::LaneQc::ConvertGenotypes {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		birdseed_file	=> { is => 'Text', doc => "Raw genotype file in 3-column birdseed format (SNP, 0/1/2, conf)", is_optional => 0, is_input => 1 },
		snp_info_file	=> { is => 'Text', doc => "SNP array information file from UCSC", is_optional => 1, is_input => 1, default => "/gscmnt/sata180/info/medseq/biodb/shared/Affy6array/snpArrayAffy6.txt" },
		output_file	=> { is => 'Text', doc => "Output file for converted genotypes", is_optional => 1, is_input => 1, is_output => 1 },		
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Converts SNP array data from birdseed format to gentoype format"                 
}

sub help_synopsis {
    return <<EOS
This command converts SNP array data from birdseed format to gentoype format
EXAMPLE:	gmt analysis lane-qc convert-genotypes --birdseed-file GenomeWideSNP_6_A07_585392.birdseed.data.txt --output-file GenomeWideSNP_6_A07_585392.converted.txt
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
	my $snp_info_file = $self->snp_info_file;

	print "Loading SNP info from $snp_info_file...\n";
	my %snp_info = load_snp_info($snp_info_file);


	print "Parsing the birdseed genotype file...\n";

	my $input = new FileHandle ($birdseed_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;

		## Only parse genotype lines ##
		if($numContents == 3 && $lineContents[1] && ($lineContents[1] eq "0" || $lineContents[1] eq "1" || $lineContents[1] eq "2" || $lineContents[1] eq "-1"))
		{
			$stats{'num_snps'}++;
			
			my ($snp_id, $call, $score) = split(/\t/, $line);
			
			if($snp_info{$snp_id})
			{
				$stats{'snps_with_info'}++;
				
				my ($chrom, $position, $strand, $allele1, $allele2) = split(/\t/, $snp_info{$snp_id});
				
				my $genotype = "";
				
				if($call eq "0")
				{
					$genotype = $allele1 . $allele1;
					$stats{'num_homozygous'}++;
				}
				elsif($call eq "1")
				{
					$genotype = $allele1 . $allele2;
					$stats{'num_heterozygous'}++;
				}
				elsif($call eq "2")
				{
					$genotype = $allele2 . $allele2;
					$stats{'num_homozygous'}++;
				}
				else
				{
					$genotype = "NN";
					$stats{'num_blank'}++;
				#	warn "Warning: unknown call value $call\n";	
				}
				
				if($genotype)
				{
					if($strand eq "-" && $genotype ne "NN")
					{
						$genotype = flip_genotype($genotype);
					}

					print join("\t", $chrom, $position, $genotype) . "\n";
				}

			}
		}
	}
	close($input);


	print $stats{'num_snps'} . " positions genotyped\n";
	print $stats{'snps_with_info'} . " had allele information in $snp_info_file\n";
	print $stats{'num_blank'} . " were blank genotypes (NN)\n";
	print $stats{'num_heterozygous'} . " were heterozygous\n";
	print $stats{'num_homozygous'} . " were homozygous\n";


	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




#############################################################
# load_snp_info - Parse allele/strand info from SNP file
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
		my $chrom = $lineContents[1];
		my $position = $lineContents[3];
		my $snp_id = $lineContents[4];
		my $strand = $lineContents[6];
		my $alleles = $lineContents[7];
		(my $allele1, my $allele2) = split(/\//, $alleles);

		$chrom =~ s/chr//;

		$info{$snp_id} = "$chrom\t$position\t$strand\t$allele1\t$allele2";
	}

	close($input);

	return(%info);
}




#############################################################
# flip_genotype - convert a genotype to the other strand
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
# flip_base - Reverse-complement a base
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

