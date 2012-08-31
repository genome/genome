
package Genome::Model::Tools::Analysis::Maf::Compare;     # rename this when you give the module file a different name <--

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

my $maf_header = "";
my %file1_results = my %file2_results = ();

class Genome::Model::Tools::Analysis::Maf::Compare {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file1	=> { is => 'Text', doc => "Original MAF file" },
		maf_file2	=> { is => 'Text', doc => "Original MAF file" },
		shared_file	=> { is => 'Text', doc => "Output maf1 entry for shared entries", is_optional => 1 },
		merged_file	=> { is => 'Text', doc => "Output non-redundant merged file", is_optional => 1 },
		unique1_file	=> { is => 'Text', doc => "Output entries unique to file 1", is_optional => 1 },
		unique2_file	=> { is => 'Text', doc => "Output entries unique to file 2", is_optional => 1 },
		verbose	=> { is => 'Text', doc => "If set to 1, give verbose output", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Compares the contents and validation statuses of two MAF files"                 
}

sub help_synopsis {
    return <<EOS
This command compares the contents and validation statuses of two MAF files
EXAMPLE:	gt analysis maf compare --maf-file1 myFile.maf --maf-file2 --theirFile.maf --output-file fixed.maf
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
	my $maf_file1 = $self->maf_file1;
	my $maf_file2 = $self->maf_file2;

	if(!(-e $maf_file1))
	{
		die "Error: MAF file 1 not found!\n";
	}
	
	if(!(-e $maf_file2))
	{
		die "Error: MAF file 2 not found!\n";
	}	

	my %stats = ();

	my %maf1 = my %maf2 = ();
	print "Loading MAF file 1...\n";
	%maf1 = load_maf_file($maf_file1);
	print "Loading MAF file 2...\n";
	%maf2 = load_maf_file($maf_file2);
	
	


	if($self->merged_file)
	{
		open(MERGED, ">" . $self->merged_file) or die "Can't open outfile: $!\n";
		print MERGED "$maf_header\n";
	}

	if($self->shared_file)
	{
		open(SHARED, ">" . $self->shared_file) or die "Can't open outfile: $!\n";
		print SHARED "$maf_header\n";
	}

	if($self->unique1_file)
	{
		open(UNIQUE1, ">" . $self->unique1_file) or die "Can't open outfile: $!\n";
		print UNIQUE1 "$maf_header\n";
	}	

	if($self->unique2_file)
	{
		open(UNIQUE2, ">" . $self->unique2_file) or die "Can't open outfile: $!\n";
		print UNIQUE2 "$maf_header\n";
	}	


	## Find samples in common ##
	
	print "SampleID\tMAF1\tMAF2\tShared\tUnique1\tUnique2\n";

	
	foreach my $sample (keys %maf1)
	{
		if($maf2{$sample})
		{
			$stats{'SamplesInCommon'}++;

			$stats{'MutationsInCommon'} = 0;
			$stats{'Mutations1'} = 0;
			$stats{'Mutations2'} = 0;
			$stats{'MutationsUnique1'} = 0;
			$stats{'MutationsUnique2'} = 0;
			
			## Get mutations from each ##

			my %mutations1 = load_mutations($maf1{$sample});
			my %mutations2 = load_mutations($maf2{$sample});
			
			## Count shared and unique-to-maf1 mutations ##
			
			foreach my $mutation_key (keys %mutations1)
			{
				$stats{'Mutations1'}++;
				
				my ($gene, $gene_id, $center, $build, $chrom, $chr_start, $chr_stop) = split(/\t/, $mutations1{$mutation_key});
				my @mutationContents1 = split(/\t/, $mutations1{$mutation_key});
				
				my $val_status1 = join("\t", $mutationContents1[24], $mutationContents1[25]);
				
				if($mutations2{$mutation_key})
				{
					$stats{'MutationsInCommon'}++;
					if($self->shared_file)
					{
						print SHARED "$mutations1{$mutation_key}\n";
					}
					
					if($self->merged_file)
					{
						print MERGED "$mutations1{$mutation_key}\n";
					}
				}
				else
				{
					print "$mutation_key\tUniqueToFile1\n" if($self->verbose);

					$stats{'MutationsUnique1'}++;
					if($self->unique1_file)
					{
						print UNIQUE1 "$mutations1{$mutation_key}\n";
					}

					if($self->merged_file)
					{
						print MERGED "$mutations1{$mutation_key}\n";
					}
				}
			}
			
			## Count unique-to-maf2 mutations ##
			
			foreach my $mutation_key (keys %mutations2)
			{
				$stats{'Mutations2'}++;

				## Check for unique-to-file2 mutations ##
				
				if(!$mutations1{$mutation_key})
				{
					$stats{'MutationsUnique2'}++;

					if($self->unique2_file)
					{
						print UNIQUE2 "$mutations2{$mutation_key}\n";
					}

					print "$mutation_key\tUniqueToFile2\n" if($self->verbose);

					if($self->merged_file)
					{
						print MERGED "$mutations2{$mutation_key}\n";
					}
				}
				
			}

			## Print output ##
			
			print "$sample\t" . $stats{'Mutations1'} . "\t" . $stats{'Mutations2'} . "\t" . $stats{'MutationsInCommon'} . "\t" . $stats{'MutationsUnique1'} . "\t" . $stats{'MutationsUnique2'} . "\n";

			## Count mutations in common, etc. ##

			$stats{'TotalMutationsInCommon'} += $stats{'MutationsInCommon'};
			$stats{'TotalMutationsUnique1'} += $stats{'MutationsUnique1'};
			$stats{'TotalMutationsUnique2'} += $stats{'MutationsUnique2'};
		}
		else
		{
			my %mutations1 = load_mutations($maf1{$sample});

			foreach my $mutation_key (keys %mutations1)
			{
				if($self->unique1_file)
				{
					print UNIQUE1 "$mutations1{$mutation_key}\n";
				}				
			}

		}
		
		exit(0);
	}


	## Get sample entries unique to file 2 ##
	
	foreach my $sample (keys %maf2)
	{

		if(!$maf1{$sample})
		{
			my %mutations2 = load_mutations($maf2{$sample});
			foreach my $mutation_key (keys %mutations2)
			{
				if($self->unique2_file)
				{
					print UNIQUE2 "$mutations2{$mutation_key}\n";
				}				
			}
		}
	}

	
	$stats{'SamplesInCommon'} = 0 if(!$stats{'SamplesInCommon'});
	$stats{'TotalMutationsInCommon'} = 0 if(!$stats{'TotalMutationsInCommon'});
	$stats{'TotalMutationsUnique1'} = 0 if(!$stats{'TotalMutationsUnique1'});
	$stats{'TotalMutationsUnique2'} = 0 if(!$stats{'TotalMutationsUnique2'});
	
	print "$stats{'SamplesInCommon'} samples in common\n";
	print "$stats{'TotalMutationsInCommon'} mutations in common\n";
	print "$stats{'TotalMutationsUnique1'} mutations unique to MAF 1\n";
	print "$stats{'TotalMutationsUnique2'} mutations unique to MAF 2\n";
	

	if($self->merged_file)
	{
		close(MERGED);
	}

	if($self->shared_file)
	{
		close(SHARED);
	}

	if($self->unique1_file)
	{
		close(UNIQUE1);
	}	

	if($self->unique2_file)
	{
		close(UNIQUE2);
	}	

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}





################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_maf_file
{
	my $maf_file = shift(@_);

	my %maf_by_sample = ();
	
	## Column index for fields in MAF file ##
	
	my %column_index = ();
	my @columns = ();
	my $header_line = ();
	
	## Parse the MAF file ##
	
	my $input = new FileHandle ($maf_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my @lineContents = split(/\t/, $line);
	
		if($lineCounter <= 2 && $line =~ "Hugo_Symbol")
		{
			$maf_header = $line;
			
			my $numContents = @lineContents;
			
			for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
			{
				if($lineContents[$colCounter])
				{
					$column_index{$lineContents[$colCounter]} = $colCounter;
				}
			}
			
			$header_line = $line;
		}
		elsif($lineCounter > 1 && !%column_index)
		{
			die "No Header in MAF file!\n";
		}
		elsif(%column_index)
		{
			## Build a record for this line, assigning all values to respective fields ##
			
			my %record = ();

			foreach my $column_name (keys %column_index)
			{
				my $column_number = $column_index{$column_name};
				$record{$column_name} = $lineContents[$column_number];
			}			

			my @tempArray = split(/\-/, $record{'Tumor_Sample_Barcode'});
			my $sample = $tempArray[0] . "-" . $tempArray[1] . "-" . $tempArray[2];# . "-" . $tempArray[3];
			
			$maf_by_sample{$sample} = "$header_line\n" if(!$maf_by_sample{$sample});
			$maf_by_sample{$sample} .= $line . "\n";
#			print "$sample\n";

		}

	}

	close($input);
	
	return(%maf_by_sample);

}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_mutations
{
	my $mutation_list = shift(@_);
	my %mutations = ();

	my %column_index = ();
	my @columns = ();
	
	my @mutation_lines = split(/\n/, $mutation_list);
	my $lineCounter = 0;
	
	foreach my $line (@mutation_lines)
	{
		$lineCounter++;
		my @lineContents = split(/\t/, $line);
	
		if($lineCounter == 1 && $line =~ "Chrom")
		{
			my $numContents = @lineContents;
			
			for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
			{
				if($lineContents[$colCounter])
				{
					$column_index{$lineContents[$colCounter]} = $colCounter;
				}
			}
		}
		elsif($lineCounter == 1)
		{
			die "No Header in MAF file!\n";
		}
		else
		{
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

			my $mutation_key = $record{'Chromosome'} . "\t" . $record{'Start_position'} . "\t" . $record{'End_position'};
			$mutations{$mutation_key} = $line;
			
			
#			$mutation .= 	$record{'Variant_Type'} . "\t" . $ref_allele . "\t" . $var_allele . "\t" . $record{'Validation_Status'} . "\t";
#			print "$sample\n";

		}

	}
	
	return(%mutations);

}




1;

