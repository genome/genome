
package Genome::Model::Tools::Analysis::Maf::FrequencyDistribution;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FrequencyDistribution - Assess the frequency distribution of variants in a MAF
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	03/01/2011 by D.K.
#	MODIFIED:	03/01/2011 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

my %stats = ();
my $max_proximity = 0;

class Genome::Model::Tools::Analysis::Maf::FrequencyDistribution {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file	=> { is => 'Text', doc => "Original MAF file" },
		output_file	=> { is => 'Text', doc => "Output file for recurrence report", is_optional => 1 },
		verbose		=> { is => 'Text', doc => "Print verbose output", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Assesses the frequency distribution of variants in a MAF file"                 
}

sub help_synopsis {
    return <<EOS
This command performs a proximity analysis on mutations in a MAF file
EXAMPLE:	gt analysis maf proximity --maf-file original.maf --output-file proximity-genes.tsv
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
	my $maf_file = $self->maf_file;

	if(!(-e $maf_file))
	{
		die "Error: MAF file not found!\n";
	}

	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
		print OUTFILE "chrom\tchr_start\tchr_stop\tref\tvar\tdbsnp_status\tvar_allele_freq\n";
	}


	my %variants = ();
	my %genotypes = ();
	


	## Column index for fields in MAF file ##
	
	my %column_index = ();
	my @columns = ();	
	
	my $input = new FileHandle ($maf_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my @lineContents = split(/\t/, $line);
	
		if($lineCounter <= 2 && $line =~ "Chrom")
		{
		
			my $numContents = @lineContents;
			
			for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
			{
				if($lineContents[$colCounter])
				{
					$column_index{$lineContents[$colCounter]} = $colCounter;
				}
			}
			
			foreach my $column (keys %column_index)
			{
				## Print out the columns as parsed ##
				#print "$column_index{$column}\t$column\n";
				$columns[$column_index{$column}] = $column;	## Save the column order ##
			}
		}
		elsif($lineCounter < 2)
		{

		}
		elsif($lineCounter > 2 && !@columns)
		{
			die "No Header in MAF file!\n";
		}
		elsif($lineCounter > 2 && @columns)
		{
			$stats{'mutations_in_file'}++;
			## Build a record for this line, assigning all values to respective fields ##
			
			my %record = ();

			foreach my $column_name (keys %column_index)
			{
				my $column_number = $column_index{$column_name};
				$record{$column_name} = $lineContents[$column_number];
			}			


			## Here's how to parse out information for this record ##
			
			my $sample = $record{'Tumor_Sample_Barcode'};

			## Parse the gene name ##
			my $chrom = $record{'Chromosome'};
			my $chr_start = $record{'Start_position'};
			my $chr_stop = $record{'End_position'};
			my $ref_allele = $record{'Reference_Allele'};
			my $var_allele = $record{'Tumor_Seq_Allele2'};
			$var_allele = $record{'Tumor_Seq_Allele1'} if($var_allele eq $ref_allele);
			my $var_type = $record{'Variant_Type'};
			my $dbsnp_rs = $record{'dbSNP_RS'};
			my $filter_status = $record{'Verification_Status'};

			## Only include variants that pass the strand filter ##

			if($filter_status ne "Strandfilter_Failed")
			{
				my $variant_key = join("\t", $chrom, $chr_start, $chr_stop, $ref_allele, $var_allele);				
				$variants{$variant_key} = join("\t", $var_type, $dbsnp_rs);

				my $genotype = $record{'Tumor_Seq_Allele1'} . $record{'Tumor_Seq_Allele2'};
				
				$genotypes{$variant_key} .= "\n" if($genotypes{$variant_key});
				$genotypes{$variant_key} .= join("\t", $sample, $genotype);
			}

		} ## Otherwise this line is non-recognizable in MAF.

	}

	close($input);	


	foreach my $variant_key (keys %variants)
	{
		$stats{'num_variants'}++;

		my ($var_type, $dbsnp_rs) = split(/\t/, $variants{$variant_key});
		
		my @genotypes = split(/\n/, $genotypes{$variant_key});
		my $num_genotypes = @genotypes;
		
		my $dbsnp_status = "novel";
		
		if($dbsnp_rs && $dbsnp_rs =~ 'rs')
		{
			$dbsnp_status = "known";
		}

		$stats{'num_' . $var_type}++;
		
		if($var_type eq "INS" || $var_type eq "DEL")
		{
			## Get size distribution ##	
		}
		elsif($var_type eq "SNP")
		{
			$stats{'num_' . $var_type . "_" . $dbsnp_status}++;
			$stats{'num_' . $var_type . "_" . $dbsnp_status . "_genotypes"} += $num_genotypes;			

			my $num_ref_alleles = my $num_var_alleles = "";

			foreach my $sample_genotype (@genotypes)
			{
				my ($sample, $genotype) = split(/\t/, $sample_genotype);
				
				if(substr($genotype, 0, 1) eq substr($genotype, 1, 1))
				{
					## Homoozygous variant ##
					$num_var_alleles += 2;
				}
				else
				{
					## Heterozygous ##
					$num_ref_alleles++;
					$num_var_alleles++;
				}
			}
			
			## Calculate MAF ##
			
			my $var_allele_freq = $num_var_alleles / (1940 * 2);
			$var_allele_freq = sprintf("%.4f", $var_allele_freq);
			
			if($self->output_file)
			{
				print OUTFILE join("\t", $variant_key, $dbsnp_status, $var_allele_freq) . "\n";				
			}

		}

	}


	print $stats{'num_variants'} . " unique variants\n";
	print $stats{'num_INS'} . " insertions\n";
	print $stats{'num_DEL'} . " deletions\n";
	print $stats{'num_SNP'} . " SNPs\n";
	print $stats{'num_SNP_known'} . " positions (" . $stats{'num_SNP_known_genotypes'} . ") were known to dbSNP\n";
	print $stats{'num_SNP_novel'} . " positions (" . $stats{'num_SNP_novel_genotypes'} . ") were novel\n";
#	foreach my $key (sort keys %stats)
#	{
#		print $stats{$key} . " $key\n";
#	}

	close(OUTFILE) if($self->output_file);


}




################################################################################################
# byGeneTranscript - sort by gene name ASC, tx name DESC, aa_position ASC
################################################################################################

sub getMutationSamples
{
	my $maf_string = shift(@_);

	my @maf_lines = split(/\n/, $maf_string);

	my $samples = "";
	my %included = ();

	foreach my $maf_line (@maf_lines)
	{
		my @lineContents = split(/\t/, $maf_line);
		my $tumor_sample = $lineContents[15];
		if(!$included{$tumor_sample})
		{
			$samples .= "\n" if($samples);
			$samples .= $tumor_sample;
			$included{$tumor_sample} = 1;
		}
	}

	return($samples);
}



1;

