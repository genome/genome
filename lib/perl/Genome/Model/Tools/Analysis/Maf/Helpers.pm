package Genome::Model::Tools::Analysis::Maf::Helpers;

use strict;
use warnings;

use above 'Genome';

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    load_annotation
    load_aa_changes
    byGeneTranscript
    getMutationSamples
);

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

sub byGeneTranscript
{
	my ($gene1, $tx1, $pos1) = split(/\t/, $a);
	my ($gene2, $tx2, $pos2) = split(/\t/, $b);

	$gene1 cmp $gene2
	or
	$tx2 cmp $tx1
	or
	$pos1 <=> $pos2;
}

sub load_aa_changes
{
	my $maf_file = shift(@_);
	my $annotation_file = shift(@_);
    my $stats = shift(@_);
	## Parse the MAF file ##

	my %mutated_positions = ();
	
	## Column index for fields in MAF file ##
	
	my %column_index = ();
	my @columns = ();


	warn "Loading variant annotations...\n";
	my %annotation = load_annotation($annotation_file);


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
			$stats->{'mutations_in_file'}++;
			## Build a record for this line, assigning all values to respective fields ##
			
			my %record = ();

			foreach my $column_name (keys %column_index)
			{
				my $column_number = $column_index{$column_name};
				$record{$column_name} = $lineContents[$column_number];
			}			


			## Here's how to parse out information for this record ##
			
			my $hugo_name = $record{'Hugo_Symbol'};
			my $tumor_sample = $record{'Tumor_Sample_Barcode'};

			## Parse the gene name ##
			my $chrom = $record{'Chromosome'};
			my $chr_start = $record{'Start_position'};
			my $chr_stop = $record{'End_position'};
			my $ref_allele = $record{'Reference_Allele'};
			my $var_allele = $record{'Tumor_Seq_Allele2'};
			$var_allele = $record{'Tumor_Seq_Allele1'} if($var_allele eq $ref_allele);
			my $var_type = $record{'Variant_Type'};
			my $variant_key = join("\t", $chrom, $chr_start, $chr_stop, $ref_allele, $var_allele);

			## Parse the gene name ##
			
#			my $gene = $record{'gene_name'};
#			my $trv_type = $record{'trv_type'};
#			my $c_position = $record{'c_position'};
#			my $aa_change = $record{'amino_acid_change'};
#			$c_position =~ s/c\.// if($c_position);
#			$aa_change =~ s/p\.// if($aa_change);
			my $gene = my $transcript_name = my $trv_type = my $c_position = my $aa_change = "";
			my $aa_position = 0;
			my $tx_start = my $tx_stop = 0;
			my $aa_position_start = my $aa_position_stop = 0;
			my $inferred_aa_start = my $inferred_aa_stop = 0;
			my $aa_pos = my $inferred_aa_pos = 0;


			## If we have no external annotation, but annotation within the maf, use it ##

			if(!$annotation{$variant_key} && $record{'gene_name'} && $record{'trv_type'} && $record{'c_position'})
			{
				$stats->{'num_annotation_in_maf'}++;
				$gene = $record{'gene_name'};
				$transcript_name = $record{'transcript_name'};
				$trv_type = $record{'trv_type'};
				$c_position = $record{'c_position'};
				$aa_change = $record{'amino_acid_change'};			

				$c_position = "" if(!$c_position);
				$aa_change = "" if(!$aa_change);
				$annotation{$variant_key} = join("\t", $gene, $trv_type, $c_position, $aa_change, $transcript_name);
			}
			else
			{
				## Otherwise this annotation came from the external file ##
				$stats->{'num_annotation_in_external_file'}++;
			}


			## If we have annotation for this variant, use it ##
			
			if($annotation{$variant_key})
			{
				$stats->{'num_with_annotation'}++;
				
				## Flags for counting if we had a nonsilent, aa-change, valid-aa-pos annotation for this one ##
				my $flag_non_silent = my $flag_aa_change = my $flag_aa_position = 0;
				
				my @annotation_lines = split(/\n/, $annotation{$variant_key});

				foreach my $annotation_line (@annotation_lines)
				{
					($gene, $trv_type, $c_position, $aa_change, $transcript_name) = split(/\t/, $annotation_line);
					
					## Proceed if we have a non-silent variant ##
		
					if($trv_type ne "silent" && $trv_type ne "rna" && $trv_type ne "intronic" && $trv_type ne "5_prime_flanking_region" && $trv_type ne "3_prime_flanking_region")
					{
						$flag_non_silent = 1;
						## Parse out aa_change if applicable and not a splice site ##
						if($aa_change && $aa_change ne "NULL" && substr($aa_change, 0, 1) ne "e")
						{
							$aa_pos = $aa_change;
							$aa_pos =~ s/[^0-9]//g;
						}
						
						## Parse out c_position if applicable ##
						
						if($c_position && $c_position ne "NULL")
						{
							## If multiple results, parse both ##
							
							if($c_position =~ '_' && !($trv_type =~ 'splice'))
							{
								($tx_start, $tx_stop) = split(/\_/, $c_position);
								$tx_start =~ s/[^0-9]//g;
								$tx_stop =~ s/[^0-9]//g;
								
								if($tx_stop < $tx_start)
								{
									$inferred_aa_start = $tx_stop / 3;
									$inferred_aa_start = sprintf("%d", $inferred_aa_start) + 1 if($tx_stop % 3);
									$inferred_aa_stop = $tx_start / 3;
									$inferred_aa_stop = sprintf("%d", $inferred_aa_stop) + 1 if($tx_start % 3);							
								}
								else
								{
									$inferred_aa_start = $tx_start / 3;
									$inferred_aa_start = sprintf("%d", $inferred_aa_start) + 1 if($tx_start % 3);
									$inferred_aa_stop = $tx_stop / 3;							
									$inferred_aa_stop = sprintf("%d", $inferred_aa_stop) + 1 if($tx_stop % 3);
								}
		
							}
							else
							{
								(my $tx_pos) = split(/[\+\-\_]/, $c_position);
								$tx_pos =~ s/[^0-9]//g;
		
								$tx_start = $tx_stop = $tx_pos;
		
								if($tx_pos)
								{
									$inferred_aa_pos = $tx_pos / 3;
									$inferred_aa_pos = sprintf("%d", $inferred_aa_pos) + 1 if($tx_pos % 3);
									$inferred_aa_start = $inferred_aa_stop = $inferred_aa_pos;
								}
								else
								{
									# IGNORE NEGATIVE TX POSITIONS (5' flanking mutations) warn "Unable to parse tx pos from $c_position\n";
								}						
							}
		
						}
			
			
						## If we inferred aa start stop, proceed with it ##
						
						if($inferred_aa_start && $inferred_aa_stop)
						{
							$aa_position_start = $inferred_aa_start;
							$aa_position_stop = $inferred_aa_stop;
							$stats->{'aa_position_inferred'}++;
						}
						## Otherwise if we inferred aa position ##
						elsif($aa_pos)
						{
							$aa_position_start = $aa_pos;
							$aa_position_stop = $aa_pos;
							$stats->{'aa_position_only'}++;
						}
						## Otherwise we were unable to infer the info ##
						else
						{
							$stats->{'aa_position_not_found'}++;
							#warn "Miss on $chrom\t$chr_start\t$chr_stop\t$ref_allele\t$var_allele\t$gene\t$c_position\t$aa_change\t$aa_pos\t$inferred_aa_pos\n";					
						}
						
						
						## Proceed if we have aa_position_start and stop ##
						
						if($aa_position_start && $aa_position_stop)
						{
							$flag_aa_position = 1;

							for(my $this_aa_pos = $aa_position_start; $this_aa_pos <= $aa_position_stop; $this_aa_pos++)
							{
								my $key = "$gene\t$transcript_name\t$this_aa_pos";

								## If we have a mutation recorded here, make sure it's not for this MAF entry ##
								if($mutated_positions{$key})
								{
									if($mutated_positions{$key} =~ $line)
									{
										## Skip because we already have it ##
									}
									## If not, append to the list of mutations ##
									else
									{
										$mutated_positions{$key} .= "\n" . $line;
									}
								}
								## Otherwise save this mutation ##
								else
								{
									$mutated_positions{$key} = $line;
								}
							}
		
						}
		
		
					} ## else it's a silent or noncoding mutation

				} ## Goto next annotation line
			
				$stats->{'num_with_aa_pos'}++ if($flag_aa_position);
			} ## Else we have no annotation for this variant in both MAF and external file

		} ## Otherwise this line is non-recognizable in MAF.

	}

	close($input);	
	

	print $stats->{'mutations_in_file'} . " mutations in file\n";
	print $stats->{'num_with_annotation'} . " with annotation available\n";
	print $stats->{'num_with_aa_pos'} . " with valid amino acid position(s)\n";
#	print $stats->{'aa_position_inferred'} . " positions were inferred\n";
#	print $stats->{'aa_position_only'} . " were notated but couldn't be inferred\n" if($stats->{'aa_position_only'});
#	print $stats->{'aa_position_not_found'} . " positions couldn't be parsed due to missing info\n";

	return(%mutated_positions);

}

sub load_annotation
{
	my $FileName = shift(@_);
	my %annotation = ();
	
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $start, $stop, $ref, $var) = split(/\t/, $line);
		my @lineContents= split(/\t/, $line);
		my $gene_name = $lineContents[6];
		my $transcript_name = $lineContents[7];
		my $trv_type = $lineContents[13];
		my $c_position = $lineContents[14];
		my $aa_change = $lineContents[15];
		
		my $variant_key = join("\t", $chrom, $start, $stop, $ref, $var);
		
		## Save or append annotation information ##
		
		if($annotation{$variant_key})
		{
			$annotation{$variant_key} .= "\n" . join("\t", $gene_name, $trv_type, $c_position, $aa_change, $transcript_name);	
		}
		else
		{
			$annotation{$variant_key} = join("\t", $gene_name, $trv_type, $c_position, $aa_change, $transcript_name);
		}

	}
	
	close($input);
	
	return(%annotation);
}

1;
