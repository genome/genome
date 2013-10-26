package Genome::Model::Tools::Vcf::ExtractSampleVariants; 

use strict;
use warnings;

use FileHandle;

use Genome;

## Declare global statistics hash ##

my %stats = ();


## Pre-define a ranking system for VEP annotation, where higher = more severe ##
my %vep_class_rank = ();
$vep_class_rank{'-'} = 				0;
$vep_class_rank{'NMD_TRANSCRIPT'} = 		0;
$vep_class_rank{'INTERGENIC'} = 		0;
$vep_class_rank{'UPSTREAM'} = 			1;
$vep_class_rank{'DOWNSTREAM'} = 		2;
$vep_class_rank{'INTRONIC'} = 			3;
$vep_class_rank{'COMPLEX_INDEL'} = 			3;
$vep_class_rank{'5PRIME_UTR'} = 		4;
$vep_class_rank{'3PRIME_UTR'} = 		5;
$vep_class_rank{'WITHIN_NON_CODING_GENE'} = 	6;
$vep_class_rank{'WITHIN_MATURE_miRNA'} = 	7;
$vep_class_rank{'PARTIAL_CODON'} = 		7;
$vep_class_rank{'CODING_UNKNOWN'} = 	        7;
$vep_class_rank{'SYNONYMOUS_CODING'} = 		8;
$vep_class_rank{'STOP_LOST'} = 			9;
$vep_class_rank{'SPLICE_SITE'} = 		10;
$vep_class_rank{'ESSENTIAL_SPLICE_SITE'} = 	11;
$vep_class_rank{'NON_SYNONYMOUS_CODING'} = 	12;
$vep_class_rank{'STOP_GAINED'} = 		13;
$vep_class_rank{'FRAMESHIFT_CODING'} = 		13;

class Genome::Model::Tools::Vcf::ExtractSampleVariants {
	is => 'Command',                       
	
	has => [
		vcf_file	=> { is => 'Text', doc => "Input VCF File" , is_optional => 0},
		vep_file	=> { is => 'Text', doc => "Input VCF File" , is_optional => 0},		
		output_file	=> { is => 'Text', doc => "Output file with sample variants" , is_optional => 0},
		genes		=> { is => 'Text', doc => "Desired genes in comma-separated list" , is_optional => 1},		
		chromosome	=> { is => 'Text', doc => "Specify a single chromosome if desired" , is_optional => 1},		
		sample_name	=> { is => 'Text', doc => "Only output variants for a specific sample" , is_optional => 1},		
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Uses VEP and VCF information to extract samples with certain variants"                 
}

sub help_synopsis {
    return <<EOS
This command uses VEP and VCF information to extract samples with certain variants.
It was developed for the BRC susceptibility project for the purpose of identifying
samples which had deleterious variants in known cancer susceptibility genes.
EXAMPLE:	gmt vcf extract-sample-variants --vcf-file my.vcf --vep-file my.vep.output --output-file my.sample.variants.TP53 --genes TP53
EOS
}

sub help_detail {
    return <<EOS 

EOS
}

sub execute {
	my $self = shift;

        my $vcf_file = $self->vcf_file;
        my $vep_file = $self->vep_file;	
	my $output_file = $self->output_file;

	my %vep = load_vep($self, $vep_file);
	
	## Open output file ##
	
	open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";
        print OUTFILE "chrom\tchr_start\tchr_stop\tref\tvar\tdbsnp_id\tdbsnp_info\tgene\ttranscript\tclass\tcdna_pos\tprotein_pos\tamino_acids\tpolyphen\tsift\tcondel\tsample_name\tsample_genotype\n";
	
	my @sample_names = ();
	my $numSamples = 0;

	$stats{'num_variants'} = $stats{'num_variants_passed'} = 0;

	## Parse the file ##

	my $input = new FileHandle ($vcf_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		if(substr($line, 0, 1) eq '#')
		{
			## HEader lines. Ignore unless samples ##
			
			if(substr($line, 0, 6) eq '#CHROM')
			{
				## Print Beagle header line ##
							
				## Get the sample names ##
				my @lineContents = split(/\t/, $line);
				my $numContents = @lineContents;
				
				for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
				{
					my $sample = $lineContents[$colCounter];
					$sample_names[$numSamples] = $sample;
					$numSamples++;

				}
				
			
				print "$numSamples samples in VCF file\n";
			}
		}
		else
		{
	                my ($chrom, $position, $id, $ref, $alt, $qual, $filter, $info, $format) = split(/\t/, $line);
			my @lineContents = split(/\t/, $line);
			my $numContents = @lineContents;
			
			if(!$self->chromosome || $chrom eq $self->chromosome)
			{
				$stats{'num_variants'}++;
				
				if($filter eq '.' || $filter eq 'PASS')
				{
					$stats{'num_variants_passed'}++;
	
					## Get the anticipated format for genotypes ##
					
					my @formatContents = split(/\:/, $format);
					my %genotype_column = ();
					my $numFormatContents = @formatContents;
					for(my $colCounter = 0; $colCounter < $numFormatContents; $colCounter++)
					{
						my $column_name = $formatContents[$colCounter];
						$genotype_column{$column_name} = $colCounter;
					}

					my @alts = split(/\,/, $alt);
					my $varAlleleCounter = 0;

					foreach my $var (@alts)
					{
						$varAlleleCounter++;
						my $variant_type = "snp";
						my $variant_key = "";
						my $chr_start = my $chr_stop = 0;
						my $allele1 = my $allele2 = "";
						
						## Check and adjust for indel ##
						if(length($ref) > 1 || length($var) > 1)
						{
							$variant_type = "indel";

							if(length($var) > length($ref))
							{
								$variant_type = "ins";
								## Insertion ##
								$allele2 = $var;
								$allele2 =~ s/$ref//;
								$allele1 = "-";
								$chr_start = $position;
								$chr_stop = $position + 1;
							}
							else
							{
								$variant_type = "del";
								## Deletion ##
								$allele1 = $ref;
								$allele1 =~ s/$var//;
								$allele2 = "-";
								my $indel_size = length($allele1);
								$chr_start = $position + 1;
								$chr_stop = $chr_start + $indel_size - 1;
							}
							
							$variant_key = join("\t", $chrom, $chr_start, $chr_stop, $allele1, $allele2);
						}
						else
						{
							## We have a SNP ##
							$chr_start = $chr_stop = $position;
							$allele1 = $ref;
							$allele2 = $var;
							$variant_key = join("\t", $chrom, $position, $position, $ref, $var);
						}
							

						if($vep{$variant_key})
						{
							$stats{'num_variants_passed_annotated'}++;
							my %allele_counts = ();
							my $numSamples = 0;

							for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
							{
								my $sample_name = $sample_names[$numSamples];
								$numSamples++;
								
								if($sample_name && (!$self->sample_name || $self->sample_name eq $sample_name))
								{
									my @genotypeContents = split(/\:/, $lineContents[$colCounter]);
									my $genotype = $genotypeContents[$genotype_column{'GT'}];
									my $coverage = $genotypeContents[$genotype_column{'DP'}];
									my $filter = $genotypeContents[$genotype_column{'FT'}];
									my $var_depth = $genotypeContents[$genotype_column{'AD'}];
									my $freq_alt = $genotypeContents[$genotype_column{'FA'}];
			
									## Get the highest var freq ##
									my $var_freq = 0;
									
									if($var_depth)
									{
										my @var_depths = split(/\,/, $var_depth);
										@var_depths = sort numericallyDesc @var_depths;
										$var_depth = $var_depths[0];
										$var_depth = 0 if($var_depth && $var_depth eq ".");
										$var_freq = $var_depth / $coverage if($coverage);				
									}
			
									
									## Only process the genotype if it has a value and is either unfiltered or for the control sample ##
									if(length($genotype) > 2 && $genotype ne '.' && ($filter eq 'PASS' || $filter eq '.'))
									{
										$stats{'sample_variants'}++;
										$genotype = convert_genotype($ref, $alt, $genotype) if($genotype ne '.');
										
										my $gt = "Missing";
										
										if(is_reference($allele1, $allele2, $genotype))
										{
											if($var_freq >= 0.05)
											{
												## Leave as missing, because we don't know ##
												$stats{'sample_variants_undercalled'}++;
											}
											else
											{
												$gt = "Ref";
											}
										}
										elsif(is_heterozygous($allele1, $allele2, $genotype))
										{
											$gt = "Het";
											$stats{'sample_variants_het'}++;												
										}
										elsif(is_homozygous($allele1, $allele2, $genotype))
										{
											$gt = "Hom";
											$stats{'sample_variants_hom'}++;												
										}
										else
										{
											warn "Unable to determine genotype ref=$ref var=$var gt=$genotype\n";
										}
										
										if($gt eq "Het" || $gt eq "Hom")
										{
											print OUTFILE join("\t", $variant_key, $id, $info, $vep{$variant_key}, $sample_name, $genotype) . "\n";											
										}

									}										
								}
								

							}

							
						}
						
						
						## Determine if we have annotation for this variant ##
						
						
						

					}
	
				}				
			}
			

			

		}

                
		

	}

	print "$stats{'target_genes_specified'} target genes specified\n" if($stats{'target_genes_specified'});	
	print "$stats{'num_variants'} variants in VCF file\n";
	print "$stats{'num_variants_passed'} passed filters\n";
	print "$stats{'num_variants_passed_annotated'} had qualifying VEP annotation\n";

	print "$stats{'sample_variants'} sample variants\n";
	print "$stats{'sample_variants_het'} were heterozygous\n";
	print "$stats{'sample_variants_hom'} were homozygous\n";
	print "$stats{'sample_variants_undercalled'} were reference but might be heterozygous\n" if($stats{'sample_variants_undercalled'});

	close($input);

	close(OUTFILE);

	return 1;
}

sub is_reference
{
	my ($ref, $var, $genotype) = @_;
	
	if(length($genotype) == 2)
	{
		if($genotype eq $ref . $ref)
		{
			return(1);
		}
	}
	else
	{
		if($genotype eq $ref . '/' . $ref)
		{
			return(1);
		}
	}
	
	return(0);
}

sub is_heterozygous
{
	my ($ref, $var, $genotype) = @_;
	
	if(length($genotype) == 2)
	{
		if($genotype eq $ref . $var || $genotype eq $var . $ref)
		{
			return(1);
		}
	}
	else
	{
		if($genotype eq $ref . '/' . $var || $genotype eq $var . '/' . $ref)
		{
			return(1);
		}
	}
	
	return(0);
}

sub is_homozygous
{
	my ($ref, $var, $genotype) = @_;
	
	if(length($genotype) == 2)
	{
		if($genotype eq $var . $var)
		{
			return(1);
		}
	}
	else
	{
		if($genotype eq $var . '/' . $var)
		{
			return(1);
		}
	}
	
	return(0);
}

sub load_vep
{
        my $self= shift(@_);
        my $annotation_file = shift(@_);
        my %annotation = ();

	my %desired_genes = ();        
	if($self->genes)
	{
		my @genes = split(/\,/, $self->genes);
		foreach my $gene (@genes)
		{
			$desired_genes{$gene} = 1;
			$stats{'target_genes_specified'}++;
		}
	}
	
        my $input = new FileHandle ($annotation_file);
        my $lineCounter = 0;
        
        while (<$input>)
        {
            chomp;
            my $line = $_;
            
            if(substr($line, 0, 1) eq '#')
            {
                ## Ignore header line ##
            }
            else
            {
                my @lineContents = split(/\t/, $line);
                my ($string) = split(/\t/, $line);
                my ($chrom, $position, $ref, $var) = split(/[\_\/]/, $string);
                my $key = join("\t", $chrom, $position, $position, $ref, $var);
        
                if($ref eq '-' || $var eq '-' || length($ref) > 1 || length($var) > 1)
                {
                    ## Load indels correctly ##
                    if($ref eq '-' || length($var) > 1)
                    {
                            ## INSERTION ##
                            my $chr_start = $position - 2;
                            my $chr_stop = $position - 1;
                            $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
                    }
                    else
                    {
                            ## DELETION ##
                            my $indel_size = length($ref);
                            my $chr_start = $position;
                            my $chr_stop = $position + $indel_size - 1;
                            $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);                       
                    }
                }
        
        
        
                my $ens_gene = $lineContents[3];
                my $class = $lineContents[6];
                my $cdna_pos = $lineContents[7];
                my $cds_pos = $lineContents[8];
                my $protein_pos = $lineContents[9];
                my $amino_acids = $lineContents[10];
                my $extra = $lineContents[13];
        
                ## Reset extra variables
                my $gene = my $polyphen = my $sift = my $condel = "";
        
                my @extraContents = split(/\;/, $extra);
                foreach my $entry (@extraContents)
                {
                        my ($key, $value) = split(/\=/, $entry);
        
                        $gene = $value if($key eq 'HGNC');
                        $polyphen = $value if($key eq 'PolyPhen');
                        $sift = $value if($key eq 'SIFT');
                        $condel = $value if($key eq 'Condel');  
                }
        
		## Determine if this gene meets the required criteria ##
		if($gene && (!$self->genes || $desired_genes{$gene}))
		{
			my @classes = split(/\,/, $class);
			foreach my $class (@classes)
			{
				if($polyphen || $sift || $condel)
				{
					$annotation{$key} .= "\n" if($annotation{$key});
					$annotation{$key} .= join("\t", $gene, $ens_gene, $class, $cdna_pos, $protein_pos, $amino_acids, $polyphen, $sift, $condel)
		#                               print join("\t", $chrom, $position, $alleles, $gene, $ens_gene, $class, $protein_pos, $amino_acids, $polyphen, $sift, $condel) . "\n";
				}
				else
				{
					$protein_pos = "-" if(!$protein_pos);
					$amino_acids = "-" if(!$amino_acids);
					$polyphen = "-" if(!$polyphen);
					$sift = "-" if(!$sift);
					$condel = "-" if(!$condel);
					$annotation{$key} .= "\n" if($annotation{$key});
					$annotation{$key} .= join("\t", $gene, $ens_gene, $class, $cdna_pos, $protein_pos, $amino_acids, $polyphen, $sift, $condel);
				}                               
			}			
		}

        
        
            }
        
        }
        
        close($input);
        
        
        ## Go through each key that has annotation and choose the top result ##
        
        foreach my $key (keys %annotation)
        {
            my @vepResults = split(/\n/, $annotation{$key});
            @vepResults = sort bySeverity @vepResults;
            my $top_result = $vepResults[0];
            $annotation{$key} = $top_result;
        }
        
        
        return(%annotation);
}

sub bySeverity
{
	my ($gene_a, $ens_gene_a, $class_a, $cdna_pos_a, $protein_pos_a, $amino_acids_a, $polyphen_a, $sift_a, $condel_a) = split(/\t/, $a);

	$polyphen_a = 0 if(!$polyphen_a || $polyphen_a eq '-');
	$sift_a = 0 if(!$sift_a || $sift_a eq '-');
	$condel_a = 0 if(!$condel_a || $condel_a eq '-');
	
	if($polyphen_a)
	{
		my @temp = split(/[\(\)]/, $polyphen_a);
		$polyphen_a = $temp[1];
	}
	if($sift_a)
	{
		my @temp = split(/[\(\)]/, $sift_a);
		$sift_a = $temp[1];
	}
	if($condel_a)
	{
		my @temp = split(/[\(\)]/, $condel_a);
		$condel_a = $temp[1];
	}
	
	my ($gene_b, $ens_gene_b, $class_b, $cdna_pos_b, $protein_pos_b, $amino_acids_b, $polyphen_b, $sift_b, $condel_b) = split(/\t/, $b);

	$polyphen_b = 0 if(!$polyphen_b || $polyphen_b eq '-');
	$sift_b = 0 if(!$sift_b || $sift_b eq '-');
	$condel_b = 0 if(!$condel_b || $condel_b eq '-');

	if($polyphen_b)
	{
		my @temp = split(/[\(\)]/, $polyphen_b);
		$polyphen_b = $temp[1];
	}
	if($sift_b)
	{
		my @temp = split(/[\(\)]/, $sift_b);
		$sift_b = $temp[1];
	}
	if($condel_b)
	{
		my @temp = split(/[\(\)]/, $condel_b);
		$condel_b = $temp[1];
	}
	
	my $fxn_code_a = fxn_class_code($class_a);
	my $fxn_code_b = fxn_class_code($class_b);
        
        die "Got no code for $class_a\n" if(!$class_a);
        die "Got no code for $class_b\n" if(!$class_b);

	## Sort by function code severity first ##
	$fxn_code_b <=> $fxn_code_a
	or
	$polyphen_b <=> $polyphen_a
	or
	$sift_b <=> $sift_a
	or
	$condel_b <=> $condel_a
}

sub fxn_class_code
{
	my $class = shift(@_);
	
	my @classes = split(/\,/, $class);
	my $num_classes = @classes;
	
	if($num_classes > 1)
	{
		@classes = sort byCode (@classes);
		$class = $classes[0];
	}
	
	foreach my $test_class (keys %vep_class_rank)
	{
		return($vep_class_rank{$test_class}) if($class eq $test_class);
	}
	
        die "No Rank provided for $class\n";
	return(0);
}

sub is_damaging
{
	my ($polyphen, $sift, $condel) = @_;
	if($polyphen && $polyphen =~ 'damaging')
	{
		return(1);
	}
	if($sift && $sift =~ 'deleterious')
	{
		return(1);
	}
	if($condel && $condel =~ 'deleterious')
	{
		return(1);
	}
	return(0);
}

sub code_to_allele
{
	my ($ref, $alt, $code) = @_;
	
	my @alt = split(/\,/, $alt);
	
	## Empty ##
	return("?") if($code eq '.');
	
	## Reference ##
	return($ref) if($code eq "0");

	## Variant ##
	return($alt[$code - 1]);	
}

sub convert_genotype
{
	my ($ref, $var, $genotype) = @_;

	return("NN") if($genotype eq '.');
	
	## Determine type of variant we're dealing with ##
	
	my $variant_type = "snp";
	(my $test_var) = split(/\,/, $var);
	
	$variant_type = "del" if(length($ref) > 1);
	$variant_type = "ins" if(length($test_var) > 1);

	## Proceed based on type of variant ##
	
	if($var =~ '\,')
	{
		my @vars = split(/\,/, $var);
		
		my ($gt1, $gt2) = split(/\//, $genotype);
		
		if($gt1 == 0)
		{
			$genotype = $ref;
		}
		else
		{
			$genotype = $vars[$gt1 - 1];
		}

		if($gt2 == 0)
		{
			$genotype .= $ref;
		}
		else
		{
			$genotype .= $vars[$gt2 - 1];
		}
	
		return($genotype);
	}
	else
	{
		if($variant_type eq 'ins')
		{
			$var = substr($var, 1, 99);
			$ref = "-";
			return($ref . '/' . $ref) if($genotype eq '0/0');
			return($ref . '/' . $var) if($genotype eq '0/1');
			return($var . '/' . $var) if($genotype eq '1/1');					
		}
		elsif($variant_type eq 'del')
		{
			$ref = substr($ref, 1, 99);
			$var = "-";
			return($ref . '/' . $ref) if($genotype eq '0/0');
			return($ref . '/' . $var) if($genotype eq '0/1');
			return($var . '/' . $var) if($genotype eq '1/1');								
		}
		else
		{
			return($ref . $ref) if($genotype eq '0/0');
			return($ref . $var) if($genotype eq '0/1');
			return($var . $var) if($genotype eq '1/1');					
		}

	}

	warn "Unable to convert $ref $var $genotype so setting to ??\n";		
	return("??");
}

sub numericallyDesc
{
	$b <=> $a;
}

1;

