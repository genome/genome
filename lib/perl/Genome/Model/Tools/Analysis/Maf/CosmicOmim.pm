
package Genome::Model::Tools::Analysis::Maf::CosmicOmim;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# CosmicOmim - Compare mutations in a MAF file to COSMIC and OMIM databases
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	08/24/2010 by D.K.
#	MODIFIED:	08/27/2010 by D.K.
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

class Genome::Model::Tools::Analysis::Maf::CosmicOmim {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file	=> { is => 'Text', doc => "Original MAF file" },
		annotation_file	=> { is => 'Text', doc => "Full annotation for variants in MAF file", is_optional => 1 },
		aa_distance	=> { is => 'Text', doc => "Maximum aa distance between mutations [1]", is_optional => 1, default => 1},
		output_cosmic	=> { is => 'Text', doc => "Output file for COSMIC hits", is_optional => 1},
		output_omim	=> { is => 'Text', doc => "Output file for OMIM hits", is_optional => 1},
		path_to_cosmic	=> { is => 'Text', doc => "Path to compiled COSMIC mutations", is_optional => 1, default => "/gscmnt/sata180/info/medseq/biodb/shared/cosmic/cosmic_will/Cosmic_Database.tsv"},
		path_to_omim	=> { is => 'Text', doc => "Path to compiled OMIM mutations", is_optional => 1, default => "/gscmnt/200/medseq/analysis/software/resources/OMIM/OMIM_Will/OMIM_aa_will.csv"},
		verbose		=> { is => 'Text', doc => "Print verbose output", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Compares mutations in a MAF file to COSMIC/OMIM databases"                 
}

sub help_synopsis {
    return <<EOS
This command compares mutations in a MAF file to COSMIC/OMIM databases
EXAMPLE:	gmt analysis maf cosmic-omim --maf-file original.maf --annotation-file original.maf.all.annotation --output-file cosmic-omim.maf
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
	my $aa_distance = $self->aa_distance;
	my $maf_file = $self->maf_file;
	my $annotation_file = $self->annotation_file;
	my $path_to_omim = $self->path_to_omim;
	my $path_to_cosmic = $self->path_to_cosmic;
	
	## Check file existence ##

	die "Error: MAF file not found at $maf_file\n" if(!(-e $maf_file));
	die "Error: COSMIC database not found at $path_to_cosmic\n" if(!(-e $path_to_cosmic));
	die "Error: OMIM database not found at $path_to_omim\n" if(!(-e $path_to_omim));


	## Load the databases ##
	
	warn "Loading the COSMIC database...\n";
	my %cosmic = load_cosmic($path_to_cosmic);

	warn "Loading the OMIM database...\n";
	my %omim = load_omim($path_to_omim);

	warn "Loading variant annotations...\n";
	my %annotation = load_annotation($annotation_file) if($annotation_file);

	## save some stats ##
	
	my %stats = ();
	$stats{'num_mutations'} = $stats{'num_with_annotation'} = $stats{'num_annotation_in_maf'} = 0;
	my %omim_gene_counts = my %cosmic_gene_counts = my %both_gene_counts = ();


	## Open output files ##
	
	if($self->output_cosmic)
	{
		open(OUTCOSMIC, ">" . $self->output_cosmic) or die "Can't open output file " . $self->output_cosmic . ": $!\n";
	}
	
	if($self->output_omim)
	{
		open(OUTOMIM, ">" . $self->output_omim) or die "Can't open output file " . $self->output_omim . ": $!\n";
	}

	## Column index for fields in MAF file ##
	
	my %column_index = ();
	my @columns = ();

	## Set up string to contain invalid annotation info ##
	
	my $invalid_annotation = "";


	## Parse the MAF file ##

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
			## Print header to output files if needed ##
			
			if($self->output_cosmic)
			{
				print OUTCOSMIC "$line\tCOSMICgene\tCOSMICchrom\tCOSMICstart\tCOSMICstop\tCOSMICamino\tCOSMICnucleotide\tCOSMICstatus\t\n";
			}
			
			if($self->output_omim)
			{
				print OUTOMIM "$line\tOMIMgene\tOMIMentry\tOMIMresidue\tOMIMaa1\tOMIMaa2\tOMIMdisorders\t\n";
			}
			
			## Parse the MAF header line to determine field locations ##	
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
			## Build a record for this line, assigning all values to respective fields ##
			
			my %record = ();

			foreach my $column_name (keys %column_index)
			{
				my $column_number = $column_index{$column_name};
				$record{$column_name} = $lineContents[$column_number];
			}			


			$stats{'num_mutations'}++;

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

			## Reset variables and counting flags ##

			my $gene = my $trv_type = my $c_position = my $aa_change = "";
			my $flag_non_silent = my $flag_aa_change = my $flag_aa_position = my $flag_matched_cosmic = my $flag_matched_omim = 0;		
			my $cosmic_match = my $omim_match = "";
			my $annotation_matching_cosmic = my $annotation_matching_omim = "";

			## If we have no external annotation, but annotation within the maf, use it ##

			if(!$annotation{$variant_key} && $record{'gene_name'} && $record{'trv_type'} && $record{'c_position'})
			{
				$stats{'num_annotation_in_maf'}++;
				$gene = $record{'gene_name'};
				$trv_type = $record{'trv_type'};
				$c_position = $record{'c_position'};
				$aa_change = $record{'amino_acid_change'};			

				$c_position = "" if(!$c_position);
				$aa_change = "" if(!$aa_change);
				$annotation{$variant_key} = join("\t", $gene, $trv_type, $c_position, $aa_change);
			}
			else
			{
				## Otherwise this annotation came from the external file ##
				$stats{'num_annotation_in_external_file'}++;
			}


			## If we have annotation for this variant, use it ##
			
			if($annotation{$variant_key})
			{
				$stats{'num_with_annotation'}++;
				
				## Attempt to match to COSMIC by chromosome and position alone ##
				my $cosmic_key = "";				
					
				if($var_type eq "SNP" || $var_type eq "DNP")
				{
					$cosmic_key = "$chrom\t$chr_start\t$chr_start\tSNP";
				}
				elsif($var_type eq "INS" || $var_type eq "DEL")
				{
					$cosmic_key = "$chrom\t$chr_start\t$chr_stop\t$var_type";				
				}
				
				## MATCH COSMIC BY CHROMOSOME AND POSITION ##
				
				if($cosmic{$cosmic_key})
				{
					$flag_matched_cosmic = 1;
					$cosmic_match = $cosmic{$cosmic_key};
				}
				else
				{
					## Loop through the positions trying to match to cosmic ##
					
					for(my $this_position = $chr_start; $this_position <= $chr_stop; $this_position++)
					{
						my $this_key = "$chrom\t$this_position";
						if($cosmic{$this_key})
						{
							$cosmic_match = $cosmic{$this_key} if(!$flag_matched_cosmic);
							$flag_matched_cosmic = 2 if(!$flag_matched_cosmic);
						}
					}
				}
				
				my @annotation_lines = split(/\n/, $annotation{$variant_key});
				
				foreach my $annotation_line (@annotation_lines)
				{
					($gene, $trv_type, $c_position, $aa_change) = split(/\t/, $annotation_line);
					
					## Proceed if we have a non-silent variant ##
		
					if($trv_type ne "silent" && $trv_type ne "rna" && $trv_type ne "intronic" && $trv_type ne "5_prime_flanking_region" && $trv_type ne "3_prime_flanking_region")
					{
						$flag_non_silent = 1;
						
						## Proceed if we got codon or AA change info ##
						if($c_position || $aa_change)
						{
							$flag_aa_change = 1;
							my ($aa_start, $aa_stop) = get_aa_positions($c_position, $aa_change);
							
							if($aa_start && $aa_stop)
							{
								$flag_aa_position = 1;
								
								## Go through each possible AA position ##
								
								for(my $this_aa_position = $aa_start - $aa_distance; $this_aa_position <= $aa_stop + $aa_distance; $this_aa_position++)
								{
									my $aa_key = join("\t", $gene, $this_aa_position);
									
									if($omim{$aa_key})
									{
										$flag_matched_omim = 1;
										if(!$omim_match)
										{
											$omim_match = $omim{$aa_key};
										}
										elsif($omim_match && $omim_match ne $omim{$aa_key})
										{
											$omim_match .= "\n" . $omim{$aa_key};	
										}
										
									}
									
									if($cosmic{$aa_key})
									{
										$cosmic_match = $cosmic{$aa_key} if(!$flag_matched_cosmic);
										$flag_matched_cosmic = 3 if(!$flag_matched_cosmic);
									}

								}
								
							}
							else
							{
								$invalid_annotation .= "$variant_key\t$gene\t$trv_type\t$c_position\t$aa_change\n";
							}
						}
					}
					
					
					
				}
			}
		

			## Count the gene matched to OMIM ##
			
			if($flag_matched_omim)
			{
				print OUTOMIM "$line\t$omim_match\n";
				my @result = split(/\t/, $omim_match);
				my $omim_gene = $result[0];
				$omim_gene_counts{$omim_gene}++;
			}
			
			
			## Count the gene matched to COSMIC ##
			
			if($flag_matched_cosmic)
			{
				print OUTCOSMIC "$line\t$cosmic_match\n";
				my @result = split(/\t/, $cosmic_match);
				my $cosmic_gene = $result[0];
				$cosmic_gene_counts{$cosmic_gene}++;
			}

			if($flag_matched_cosmic && $flag_matched_omim)
			{
				(my $gene) = split(/\t/, $cosmic_match);
				$both_gene_counts{$gene}++;
			}

			$stats{'num_non_silent'}++ if($flag_non_silent);
			$stats{'num_with_aa_change'}++ if($flag_aa_change);
			$stats{'num_with_aa_position'}++ if($flag_aa_position);
			$stats{'num_matched_omim'}++ if($flag_matched_omim);
			$stats{'num_matched_cosmic'}++ if($flag_matched_cosmic);
			$stats{'num_matched_both'}++ if($flag_matched_omim && $flag_matched_cosmic);
			$stats{'num_matched_cosmic_by_variant'}++ if($flag_matched_cosmic && $flag_matched_cosmic == 1);
			$stats{'num_matched_cosmic_by_coordinate'}++ if($flag_matched_cosmic && $flag_matched_cosmic == 2);
			$stats{'num_matched_cosmic_by_position'}++ if($flag_matched_cosmic && $flag_matched_cosmic == 3);
		}

	}

	close($input);	


	## Close output tfiles ##

	if($self->output_cosmic)
	{
		close(OUTCOSMIC);
	}
	
	if($self->output_omim)
	{
		close(OUTOMIM);
	}

	
	if($self->verbose)
	{
		warn "Invalid amino acid positions were the following:\n$invalid_annotation\n";		
	}


	print $stats{'num_mutations'} . " mutations in the MAF file\n";
	print $stats{'num_with_annotation'} . " had annotation information\n";
	print "\t" . $stats{'num_annotation_in_external_file'} . " had annotation in the external file\n";
	print "\t" . $stats{'num_annotation_in_maf'} . " instead had top annotation in the MAF file\n";

	print $stats{'num_non_silent'} . " non-silent mutations\n";
	print $stats{'num_with_aa_change'} . " with amino acid change\n";
	print $stats{'num_with_aa_position'} . " had valid amino acid position (turn on verbose to see failures)\n";

	print $stats{'num_matched_cosmic'} . " matched mutations in COSMIC\n";
	if($self->verbose)
	{
		print "\t" . $stats{'num_matched_cosmic_by_variant'} . " matched coordinates and variant type\n";
		print "\t" . $stats{'num_matched_cosmic_by_coordinate'} . " matched coordinates\n";
		print "\t" . $stats{'num_matched_cosmic_by_position'} . " matched amino acid position\n";
	}

	print $stats{'num_matched_omim'} . " matched mutations in OMIM\n";

	print $stats{'num_matched_both'} . " matched both databases\n";

	if($self->verbose)
	{
		## Print matched COSMIC genes ##
		
		my $cosmic_counts = "";
		foreach my $gene (sort keys %cosmic_gene_counts)
		{
			$cosmic_counts .= ", " if($cosmic_counts);
			$cosmic_counts .= "$cosmic_gene_counts{$gene} $gene";
		}
		
		print "COSMIC Gene Counts: ";	
		print "$cosmic_counts\n\n";


		## Print matched OMIM genes ##
		
		my $omim_counts = "";
		foreach my $gene (sort keys %omim_gene_counts)
		{
			$omim_counts .= ", " if($omim_counts);
			$omim_counts .= "$omim_gene_counts{$gene} $gene";
		}

		print "OMIM Gene Counts: ";		
		print "$omim_counts\n";



		## Print matched OMIM+COSMIC genes ##
		
		my $both_counts = "";
		foreach my $gene (sort keys %both_gene_counts)
		{
			$both_counts .= ", " if($both_counts);
			$both_counts .= "$both_gene_counts{$gene} $gene";
		}

		print "COSMIC+OMIM Gene Counts: ";		
		print "$both_counts\n";	

	}



}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub get_aa_positions
{
	my ($c_position, $aa_change) = @_;
	
	## Remove unnecessary characters from start of these guys ##
	$c_position =~ s/c\.// if($c_position);
	$aa_change =~ s/p\.// if($aa_change);
	
	my $inferred_aa_start = my $inferred_aa_stop = 0;
	
	## Proceed using codon positino ##
	
	if($c_position && $c_position ne "NULL")
	{
		if(substr($c_position, 0, 1) eq '-')
		{
			return(0, 0);	## Ignore upstream stuff ##
		}
		
		## If multiple results, parse both ##
		
		if($c_position =~ '_')
		{
			my ($tx_start, $tx_stop) = split(/\_/, $c_position);
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

			if($tx_pos)
			{
				my $inferred_aa_pos = $tx_pos / 3;
				$inferred_aa_pos = sprintf("%d", $inferred_aa_pos) + 1 if($tx_pos % 3);
				$inferred_aa_start = $inferred_aa_stop = $inferred_aa_pos;
			}
			else
			{
				warn "Unable to parse tx pos from $c_position\n";
			}						
		}

	}
	
	## If we don't have codon but have AA change information, make use of it ##
	
	elsif($aa_change && $aa_change ne "NULL" && substr($aa_change, 0, 1) ne "e")
	{
		(my $aa_pos) = split(/[\_\-]/, $aa_change);
		$aa_pos =~ s/[^0-9]//g;
		
		if($aa_pos)
		{
			$inferred_aa_start = $inferred_aa_stop = $aa_pos;
		}
	}
	
	return($inferred_aa_start, $inferred_aa_stop);
}





#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

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
		my $trv_type = $lineContents[13];
		my $c_position = $lineContents[14];
		my $aa_change = $lineContents[15];
		
		my $variant_key = join("\t", $chrom, $start, $stop, $ref, $var);
		
		## Save or append annotation information ##
		
		if($annotation{$variant_key})
		{
			$annotation{$variant_key} .= "\n" . join("\t", $gene_name, $trv_type, $c_position, $aa_change);	
		}
		else
		{
			$annotation{$variant_key} = join("\t", $gene_name, $trv_type, $c_position, $aa_change);
		}

	}
	
	close($input);
	
	return(%annotation);
}




#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub load_cosmic
{
	my $FileName = shift(@_);
	my %cosmic = ();
	
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		if($lineCounter == 1 && $ARGV[1])
		{
		}
		elsif($lineCounter > 1)
		{
			my ($gene, $chrom, $chr_start, $chr_stop, $aa_change, $codon_change, $somatic_status) = split(/\t/, $line);
			my $variant_type = "";
			my $allele1 = my $allele2 = "";

			$stats{'cosmic_entries'}++;

			## Make a key for the aa_position ##
			$aa_change =~ s/p\.//;
			(my $aa_position) = split(/[\-\_\.]/, $aa_change);
			$aa_position =~ s/[^0-9]//g;

			my $aa_key = "$gene\t$aa_position";
			$cosmic{$aa_key} = $line;

			my $variant_key = "";
			
			if($codon_change && $codon_change =~ '>')
			{
				$codon_change =~ s/c\.//;
				my $alleles = $codon_change;
				$alleles =~ s/[0-9]//g;
				($allele1, $allele2) = split(/\>/, $alleles);
				$variant_type = "SNP";				
			}
			elsif($codon_change =~ 'ins')
			{
				$variant_type = "INS";
			}
			elsif($codon_change =~ 'del')
			{
				$variant_type = "DEL";
			}

			if($variant_type)
			{
				$variant_key = join("\t", $chrom, $chr_start, $chr_stop, $variant_type);
				$cosmic{$variant_key} = $line;
			}
			
			if($chrom && $chr_start)
			{
				$cosmic{"$chrom\t$chr_start"} = $line;
			}

#			print "$chrom\t$chr_start\t$chr_stop\t$allele1\t$allele2\n";
		}
	}
	
	close($input);
	
	return(%cosmic);
}






#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub load_omim
{
	my $FileName = shift(@_);
	my %omim = ();
	
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		if($lineCounter == 1 && $ARGV[1])
		{
		}
		if($lineCounter > 1)
		{
			my ($gene, $entry, $position, $aa1, $aa2) = split(/\t/, $line);
			my @lineContents = split(/\t/, $line);
			my $omim_result = "";
			$omim_result = $lineContents[5] if($lineContents[5]);
			$omim_result = $lineContents[6] if($lineContents[6]);

			## Filter the omim result to remove junk ##
			
			my @result = split(/\s+/, $omim_result);
			my $brief_result = "";
			
			foreach my $word (@result)
			{
				my $test = $word;
				$test =~ s/[0-9]//g;
				if(!$test)
				{
					## Skip this junk ##
				}
				else
				{
					$test =~ s/\(//;
					$test =~ s/\)//;
	    
					if($test && $test ne ";")
					{
						$brief_result .= " " if($brief_result);
						$brief_result .= $word;											
					}
					elsif($test eq ";")
					{
						$brief_result .= ";";
					}
				}
			}
			
			## Fix comma-semicolon intersections ##
			my $str = ",;";
			$brief_result =~ s/$str/\;/g;

			## Remove square brackets ##
			
			$brief_result =~ s/\[//g;
			$brief_result =~ s/\]//g;

			## Remove braces ##
			
			$brief_result =~ s/\{//g;
			$brief_result =~ s/\}//g;

			
			## Fix the problem ##
			$line = join("\t", $gene, $entry, $position, $aa1, $aa2, $brief_result);


			## Save AA position, gene, AA1, AA2 ##
			
			my $key = join("\t", $gene, $position, $aa1, $aa2);
			$omim{$key} = $line;

			## Save just AA position and gene too ##

			$key = join("\t", $gene, $position);
			$omim{$key} = $line;


		}
	}
	
	close($input);
	
	return(%omim);
}




1;

