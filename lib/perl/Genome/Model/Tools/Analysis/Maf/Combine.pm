
package Genome::Model::Tools::Analysis::Maf::Combine;     # rename this when you give the module file a different name <--

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
my %gene_ids = ();

class Genome::Model::Tools::Analysis::Maf::Combine {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file1	=> { is => 'Text', doc => "Original MAF file" },
		maf_file2	=> { is => 'Text', doc => "Original MAF file" },
		output_file	=> { is => 'Text', doc => "Original MAF file", is_optional => 1 },
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

	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
	}

	my %stats = ();

	my %maf1 = my %maf2 = ();
	print "Loading MAF file 1...\n";
	%maf1 = load_maf_file($maf_file1);
	print "Loading MAF file 2...\n";
	%maf2 = load_maf_file($maf_file2);
	
	my $maf_header = "";
	$maf_header = `head -1 $maf_file1`;
	chomp($maf_header);
	
	my $combined_maf = $maf_header . "\n";
	
	my @mafLines = ();
	my $mafLineCounter = 0;
	
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
				
				if($mutations2{$mutation_key})
				{
					$stats{'MutationsInCommon'}++;
				}
				else
				{
					$stats{'MutationsUnique1'}++;
				}

				## Print to combined MAF ##

				if($self->output_file)
				{
					$mafLines[$mafLineCounter] = $mutations1{$mutation_key};
					$mafLineCounter++;
#					print OUTFILE "$mutations1{$mutation_key}\n";
				}

			}
			
			## Count unique-to-maf2 mutations ##
			
			foreach my $mutation_key (keys %mutations2)
			{
				$stats{'Mutations2'}++;
				
				if(!$mutations1{$mutation_key})
				{
					$stats{'MutationsUnique2'}++;

					## Print new entries to combined MAF ##
	
					if($self->output_file)
					{
						$mafLines[$mafLineCounter] = $mutations2{$mutation_key};
						$mafLineCounter++;
#						print OUTFILE "$mutations2{$mutation_key}\n";
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
	}
	
	$stats{'SamplesInCommon'} = 0 if(!$stats{'SamplesInCommon'});
	$stats{'TotalMutationsInCommon'} = 0 if(!$stats{'TotalMutationsInCommon'});
	$stats{'TotalMutationsUnique1'} = 0 if(!$stats{'TotalMutationsUnique1'});
	$stats{'TotalMutationsUnique2'} = 0 if(!$stats{'TotalMutationsUnique2'});
	
	print "$stats{'SamplesInCommon'} samples in common\n";
	print "$stats{'TotalMutationsInCommon'} mutations in common\n";
	print "$stats{'TotalMutationsUnique1'} mutations unique to MAF 1\n";
	print "$stats{'TotalMutationsUnique2'} mutations unique to MAF 2\n";
	

	

	## Sort the output ##
	
	if($self->output_file)
	{
		print OUTFILE "$maf_header\n";

		@mafLines = sort byChrPos @mafLines;
		
		foreach my $mafLine (@mafLines)
		{
			print OUTFILE "$mafLine\n";
		}
	}

	sub byChrPos
	{
		my @temp = split(/\t/, $a);
		my $chrom_a = $temp[4];
		my $start_a = $temp[5];
		my $stop_a = $temp[6];

		@temp = split(/\t/, $b);
		my $chrom_b = $temp[4];
		my $start_b = $temp[5];
		my $stop_b = $temp[6];
	
		$chrom_a =~ s/X/23/;
		$chrom_a =~ s/Y/24/;
		$chrom_a =~ s/MT/24/;

		$chrom_b =~ s/X/23/;
		$chrom_b =~ s/Y/24/;
		$chrom_b =~ s/MT/24/;

		$chrom_a <=> $chrom_b
		or
		$start_a <=> $start_b
		or
		$stop_a <=> $stop_b;		
	}


	if($self->output_file)
	{
		close(OUTFILE);
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
			
			$header_line = $line;
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

			## Get gene ID if possible ##
			my $gene_name = $record{'Hugo_Symbol'};

			if($record{'Entrez_Gene_Id'} && $record{'Entrez_Gene_Id'} > 0)
			{
				my $gene_id = $record{'Entrez_Gene_Id'};
				$gene_ids{$gene_name} = $gene_id;
			}
			elsif($gene_ids{$gene_name})
			{
				## Grab the gene name ##
				my $gene_id = $gene_ids{$gene_name};
				
				my $search_string = "$gene_name\t0\t";
				my $replace_string = "$gene_name\t$gene_id\t";

				$line =~ s/$search_string/$replace_string/;
			}


			## Get the ref and variant alleles ##
			
			my $ref_allele = $record{'Reference_Allele'};
			my $var_allele = "";
			
			$var_allele = $record{'Tumor_Seq_Allele1'} if($record{'Tumor_Seq_Allele1'} ne $ref_allele);
			$var_allele = $record{'Tumor_Seq_Allele2'} if($record{'Tumor_Seq_Allele2'} ne $ref_allele);

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

