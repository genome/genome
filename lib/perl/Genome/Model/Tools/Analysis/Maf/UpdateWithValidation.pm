
package Genome::Model::Tools::Analysis::Maf::UpdateWithValidation;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# UpdateWithValidation - Perform a proximity analysis on mutations in the MAF file.
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	08/24/2010 by D.K.
#	MODIFIED:	08/24/2010 by D.K.
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

class Genome::Model::Tools::Analysis::Maf::UpdateWithValidation {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file	=> { is => 'Text', doc => "Original MAF file", is_input => 1 },
		validation_file	=> { is => 'Text', doc => "Validation data in chrom,start,stop,ref,var,tumor_sample,review_code format", is_input => 1 },		
		output_file	=> { is => 'Text', doc => "Output name for updated MAF file", is_output => 1 },
		verbose		=> { is => 'Text', doc => "Print verbose output", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Updates a MAF with validation status"                 
}

sub help_synopsis {
    return <<EOS
This command updates a MAF file with validation information
EXAMPLE:	gmt analysis maf update-with-validation --maf-file myMaf.tsv --validation-file myVal.tsv --output-file myMafUpdated.tsv
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

	my $output_file = $self->output_file;

	my %validation = load_validation($self->validation_file);

	## Column index for fields in MAF file ##
	
	my %column_index = ();
	my @columns = ();

	my %snp_is_dnp = ();

	## open outfile ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";

	my $input = new FileHandle ($maf_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;
		
		if($lineCounter <= 2 && $line =~ "Chrom")
		{
			
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
			
			print OUTFILE "$line\n";
		}
		elsif($lineCounter < 2)
		{
			print OUTFILE "$line\n";
		}
		elsif($lineCounter > 2 && !@columns)
		{
			die "No Header in MAF file!\n";
		}
		elsif($lineCounter > 2 && @columns)
		{
			$stats{'num_mutations'}++;
			
			## Build a record for this line, assigning all values to respective fields ##
			
			my %record = ();

			foreach my $column_name (keys %column_index)
			{
				my $column_number = $column_index{$column_name};
				$record{$column_name} = $lineContents[$column_number];
			}			


			## Here's how to parse out information for this record ##
			
			my $chrom = $record{'Chromosome'};
			my $chr_start = $record{'Start_position'};
			my $chr_stop = $record{'End_position'};			
			my $tumor_sample = $record{'Tumor_Sample_Barcode'};
			my $variant_type = $record{'Variant_Type'};
			my $ref = $record{'Reference_Allele'};			
			my $var = $record{'Tumor_Seq_Allele2'};
			$var = $record{'Tumor_Seq_Allele1'} if($var eq $ref);

			my @temp = split(/\-/, $tumor_sample);
			my $patient = join("-", "TCGA", $temp[1], $temp[2]);

			my $mutation_key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var, $patient);

			$stats{'num_' . $variant_type}++;

			if($validation{$mutation_key})
			{
				my $val_code = $validation{$mutation_key};
				$stats{'num_' . $variant_type . '_val_' . $val_code}++;
				
				my $somatic_status = my $val_status = my $tumor_val_allele1 = my $tumor_val_allele2 = my $normal_val_allele1 = my $normal_val_allele2 = "";
				
				if($val_code eq "S")
				{
					$somatic_status = "Somatic";
					$val_status = "Valid";
					$normal_val_allele1 = $ref;
					$normal_val_allele2 = $ref;
					$tumor_val_allele1 = $ref;
					$tumor_val_allele2 = $var;					
				}
				elsif($val_code eq "LOH")
				{
					$somatic_status = "LOH";
					$val_status = "Valid";
					$normal_val_allele1 = $ref;
					$normal_val_allele2 = $var;
					$tumor_val_allele1 = $var;
					$tumor_val_allele2 = $var;										
				}
				elsif($val_code eq "G")
				{
					$somatic_status = "Germline";
					$val_status = "Valid";
					$normal_val_allele1 = $ref;
					$normal_val_allele2 = $var;
					$tumor_val_allele1 = $ref;
					$tumor_val_allele2 = $var;										
				}
				elsif($val_code eq "O" || $val_code eq "WT")
				{
					$somatic_status = "None";
					$val_status = "Wildtype";
					$normal_val_allele1 = $ref;
					$normal_val_allele2 = $ref;
					$tumor_val_allele1 = $ref;
					$tumor_val_allele2 = $ref;										
				}
				
				my $new_line = "";
				
				for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
				{
					$new_line .= "\t" if($colCounter > 0);
					
					if($columns[$colCounter] eq 'Tumor_Validation_Allele1')
					{
						$new_line .= $tumor_val_allele1;
					}
					elsif($columns[$colCounter] eq 'Tumor_Validation_Allele2')
					{
						$new_line .= $tumor_val_allele2;
					}					
					elsif($columns[$colCounter] eq 'Match_Norm_Validation_Allele1')
					{
						$new_line .= $normal_val_allele1;
					}					
					elsif($columns[$colCounter] eq 'Match_Norm_Validation_Allele2')
					{
						$new_line .= $normal_val_allele2;
					}
					elsif($columns[$colCounter] eq 'Mutation_Status')
					{
						$new_line .= $somatic_status;
					}
					elsif($columns[$colCounter] eq 'Validation_Status')
					{
						$new_line .= $val_status;
					}					
					elsif($columns[$colCounter] eq 'Validation_Method')
					{
						$new_line .= "Illumina_Capture_gDNA";
					}
					else
					{
						$new_line .= $lineContents[$colCounter];
					}
					
				}
				
				print OUTFILE "$new_line\n";
				$line = $new_line;
			}
			else
			{
				print OUTFILE "$line\n";

			}


			
		}

	}

	close($input);
	
	close(OUTFILE);

	foreach my $key (sort keys %stats)
	{
		print $stats{$key} . "\t" . $key . "\n";
	}
	
	return(1);
}




################################################################################################
# Load all SNV positions
#
################################################################################################

sub load_validation
{                               # replace with real execution logic.
	my $FileName = shift(@_);

	my %review = ();
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my ($chrom, $chr_start, $chr_stop, $ref, $var, $tumor_sample, $review_code) = split(/\t/, $line);

		my @temp = split(/\-/, $tumor_sample);
		my $patient = join("-", "TCGA", $temp[1], $temp[2]);
		
		my $variant_key = my $mutation_key = "";

		if($ref eq '*' || $ref eq '0' || $ref eq '-' || length($var) > 1)
		{
			## INSERTION ##
			$ref = "-";
#			$chr_start++;
#			$chr_stop++;
			$chr_start-- if($chr_stop eq $chr_start);
			$variant_key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
			$mutation_key = join("\t", $variant_key, $patient);
		}
		elsif($var eq '*' || $var eq '0' || $var eq '-' || length($ref) > 1)
		{
			## DELETION ##
			$var = "-";
			$chr_start++;
			$variant_key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
			$mutation_key = join("\t", $variant_key, $patient);                                

		}
		else
		{
			## SNV ##
			$chr_start = $chr_stop if($chr_stop > $chr_start);
			$variant_key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
			$mutation_key = join("\t", $variant_key, $patient);

		}
		
		if($review_code ne "NA" && $review_code ne "A")
		{
			$review{$variant_key} = $review_code;
			$review{$mutation_key} = $review_code;
		}		
		
	}
	
	close($input);
	return(%review);

}




1;

