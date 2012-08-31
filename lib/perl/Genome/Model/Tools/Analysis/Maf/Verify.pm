
package Genome::Model::Tools::Analysis::Maf::Verify;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# VerifyMaf - Align reads with SSAHA2 or other aligner
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	02/25/2009 by D.K.
#	MODIFIED:	02/25/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

my %file1_results = my %file2_results = ();

class Genome::Model::Tools::Analysis::Maf::Verify {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file	=> { is => 'Text', doc => "Original MAF file" },
		output_file	=> { is => 'Text', doc => "Original MAF file", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Finds inconsistent and/or missing information in maf file "                 
}

sub help_synopsis {
    return <<EOS
This command finds inconsistent or missing information in a maf file
EXAMPLE:	gt analysis maf verify --maf-file original.maf --output-file fixed.maf
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
		die "Error: MAF file 1 not found!\n";
	}

	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
	}

	my %stats = ();

	## Column index for fields in MAF file ##
	
	my %column_index = ();
	my @columns = ();

	## Parse the MAF file ##
	
	my $input = new FileHandle ($maf_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my @lineContents = split(/\t/, $line);
	
		if($lineCounter == 1 && !($line =~ "Chrom"))
		{
			## Print the version number ##
			if($self->output_file)
			{
				print OUTFILE "$line\n";
			}
		}
		elsif($lineCounter <= 2 && $line =~ "Chrom")
		{
			if($self->output_file)
			{
				print OUTFILE "$line\n";
			}
			
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
				$columns[$column_index{$column}] = $column;	## Save the column order ##
			}
		}
		elsif($lineCounter == 1)
		{
			die "No Header in MAF file!\n";
		}
		else
		{
			$stats{'num_mutations'}++;
			## Build a record for this line, assigning all values to respective fields ##
			
			my %record = ();

			foreach my $column_name (keys %column_index)
			{
				my $column_number = $column_index{$column_name};
				$record{$column_name} = $lineContents[$column_number];
			}			


			## Get the ref and variant alleles ##
			
			my $ref_allele = $record{'Reference_Allele'};
			my $var_allele = "";
			
			$var_allele = $record{'Tumor_Seq_Allele1'} if($record{'Tumor_Seq_Allele1'} ne $ref_allele);
			$var_allele = $record{'Tumor_Seq_Allele2'} if($record{'Tumor_Seq_Allele2'} ne $ref_allele);


			if($ref_allele && $var_allele && $ref_allele ne $var_allele)
			{
				## Check matched normal alleles ##
				
				if(!$record{'Match_Norm_Seq_Allele1'} || ($record{'Match_Norm_Seq_Allele1'} ne $ref_allele && $record{'Match_Norm_Seq_Allele1'} ne $var_allele))
				{
#					print "WARNING: No normal allele1: ";
#					print $record{'Chromosome'} . "\t" . $record{'Start_position'} . "\t" . $record{'End_position'} . "\t";
#					print $record{'Variant_Type'} . "\t" . $record{'Validation_Status'} . "\n";					
					$record{'Match_Norm_Seq_Allele1'} = $ref_allele;
				}
				
				if(!$record{'Match_Norm_Seq_Allele2'} || ($record{'Match_Norm_Seq_Allele2'} ne $ref_allele && $record{'Match_Norm_Seq_Allele2'} ne $var_allele))
				{
#					print "WARNING: No normal allele2: ";
#					print $record{'Chromosome'} . "\t" . $record{'Start_position'} . "\t" . $record{'End_position'} . "\t";
#					print $record{'Variant_Type'} . "\t" . $record{'Validation_Status'} . "\n";
					$record{'Match_Norm_Seq_Allele2'} = $var_allele;
				}
				
				## Check Mutation type ##
				
				if($ref_allele eq "*" || $var_allele eq "*")
				{
					print "WARNING: invalid ref/var allele: $ref_allele/$var_allele\n";
				}
				
				## Check validation alleles ##
				if(lc($record{'Validation_Status'}) ne "unknown")
				{
#					if(!$record{'Match_Norm_Validation_Allele1'} || !$record{'Match_Norm_Validation_Allele2'} || !$record{'Tumor_Validation_Allele1'} || !$record{'Tumor_Validation_Allele2'})
#					{
#						print "WARNING: Missing Validation alleles: ($record{'Match_Norm_Validation_Allele1'} || $record{'Match_Norm_Validation_Allele2'} || $record{'Tumor_Validation_Allele1'} || $record{'Tumor_Validation_Allele2'})";
#						print $record{'Chromosome'} . "\t" . $record{'Start_position'} . "\t" . $record{'End_position'} . "\t";
#						print $record{'Variant_Type'} . "\t" . $record{'Validation_Status'} . "\t";
#						print "\n";
#					}
#					elsif($record{'Validation_Status'} eq "Somatic")
#					{
#						## Normal should be homozygous-reference; Tumor should be het or hom variant
#						if($record{'Match_Norm_Validation_Allele1'} eq $ref_allele && $record{'Match_Norm_Validation_Allele2'} eq $ref_allele)
#						{
#							if($record{'Tumor_Validation_Allele1'} eq $var_allele || $record{'Tumor_Validation_Allele2'} eq $var_allele)
#							{
#								
#							}
#						}
#					}
					if($record{'Validation_Status'} eq "Valid" && $record{'Mutation_Status'} eq "Somatic")
					{
						$record{'Match_Norm_Validation_Allele1'} = $ref_allele;
						$record{'Match_Norm_Validation_Allele2'} = $ref_allele;
						$stats{'validation_alleles_fixed'}++ if(!$record{'Tumor_Validation_Allele1'} || !$record{'Tumor_Validation_Allele2'});
						$record{'Tumor_Validation_Allele1'} = $record{'Tumor_Seq_Allele1'};
						$record{'Tumor_Validation_Allele2'} = $record{'Tumor_Seq_Allele2'};
					}
					elsif($record{'Validation_Status'} eq "Valid" && $record{'Mutation_Status'} eq "LOH")
					{
						$record{'Match_Norm_Validation_Allele1'} = $record{'Match_Norm_Seq_Allele1'} if(!$record{'Match_Norm_Validation_Allele1'} || ($record{'Match_Norm_Validation_Allele1'} ne $ref_allele && $record{'Match_Norm_Validation_Allele1'} ne $var_allele));
						$record{'Match_Norm_Validation_Allele2'} = $record{'Match_Norm_Seq_Allele2'} if(!$record{'Match_Norm_Validation_Allele2'} || ($record{'Match_Norm_Validation_Allele2'} ne $ref_allele && $record{'Match_Norm_Validation_Allele2'} ne $var_allele));
						$record{'Tumor_Validation_Allele1'} = $ref_allele if(!$record{'Tumor_Validation_Allele1'} || ($record{'Tumor_Validation_Allele1'} ne $ref_allele && $record{'Tumor_Validation_Allele1'} ne $var_allele));
						$record{'Tumor_Validation_Allele2'} = $var_allele if(!$record{'Tumor_Validation_Allele2'} || ($record{'Tumor_Validation_Allele2'} ne $ref_allele && $record{'Tumor_Validation_Allele2'} ne $var_allele));												
					}
					elsif($record{'Validation_Status'} eq "Valid" && $record{'Mutation_Status'} eq "Germline")
					{
						$record{'Match_Norm_Validation_Allele1'} = $record{'Tumor_Seq_Allele1'};
						$record{'Match_Norm_Validation_Allele2'} = $record{'Tumor_Seq_Allele2'};
						$record{'Tumor_Validation_Allele1'} = $record{'Tumor_Seq_Allele1'};
						$record{'Tumor_Validation_Allele2'} = $record{'Tumor_Seq_Allele2'};
					}
					elsif(lc($record{'Validation_Status'}) eq "wildtype")
					{
						$record{'Match_Norm_Validation_Allele1'} = $ref_allele;
						$record{'Match_Norm_Validation_Allele2'} = $ref_allele; 
						$record{'Tumor_Validation_Allele1'} = $ref_allele;
						$record{'Tumor_Validation_Allele2'} = $ref_allele;
					}
					else
					{
						print "WARNING: Unknown val status: ";
						print $record{'Chromosome'} . "\t" . $record{'Start_position'} . "\t" . $record{'End_position'} . "\t";
						print $record{'Variant_Type'} . "\t" . $record{'Validation_Status'} . "\t" . $record{'Mutation_Status'} . "\n";
					}

				}
			}
			else
			{
				print "WARNING: No allele2: ";
				print $record{'Chromosome'} . "\t" . $record{'Start_position'} . "\t" . $record{'End_position'} . "\t";
				print $record{'Variant_Type'} . "\t" . $record{'Validation_Status'} . "\t";
				print "Ref=$ref_allele\tVar=$var_allele\n";				
			}


			## Build newline ##
			
			my $newline = "";
			
			my $numColumns = @columns;
			
			for(my $colCounter = 0; $colCounter < $numColumns; $colCounter++)
			{
				my $column_name = $columns[$colCounter];
				my $column_value = $record{$column_name};
				$newline .= "$column_value\t";					
			}
			
			if($self->output_file)
			{
				print OUTFILE "$newline\n";
			}

		}

	}

	close($input);	
		

	print $stats{'num_mutations'} . " mutations in the MAF file\n";
	print $stats{'validation_alleles_fixed'} . " were missing validation alleles and were fixed\n";


	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




#############################################################
# parse_file - takes input file and parses it
#
#############################################################

sub check_allele
{
	(my $ref_base, my $var_base, my $this_base, my $other_base) = @_;

	if($this_base && ($this_base eq "A" || $this_base eq "C" || $this_base eq "G" || $this_base eq "T"))
	{
		return($this_base);
	}
	elsif($this_base)
	{
		if($other_base eq $ref_base)
		{
			return($var_base);
		}
		else
		{
			return($ref_base);
		}
	}
	else
	{
		return($var_base);
	}
		
}





#############################################################
# parse_file - takes input file and parses it
#
#############################################################

sub check_platform
{
	(my $platform) = @_;


	if($platform =~ "454" && $platform =~ "3730")
	{
		$platform = "Roche 454 / ABI 3730xl";
	}
	elsif($platform =~ "454")
	{
		$platform = "Roche 454";
	}
	elsif($platform =~ "Illumina")
	{
		$platform = "Illumina GAIIx";
	}
	elsif($platform =~ "3730")
	{
		$platform = "ABI 3730xl";
	}
	else
	{
		$platform = "Roche 454";
	}
	
	return($platform);
}







1;

