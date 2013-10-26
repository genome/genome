
package Genome::Model::Tools::Analysis::Mendelian::RecessiveAlleles;

use strict;
use warnings;

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

my %damaging_classes = ();
$damaging_classes{'missense'} = 1;
$damaging_classes{'nonsense'} = 1;
$damaging_classes{'splice_site'} = 1;
$damaging_classes{'splice_site_ins'} = 1;
$damaging_classes{'splice_site_del'} = 1;
$damaging_classes{'nonstop'} = 1;
$damaging_classes{'frame_shift_ins'} = 1;
$damaging_classes{'frame_shift_del'} = 1;

use FileHandle;

use Genome;
use Genome::Model::Tools::Capture::Helpers qw(
    byChrPos
);

my $num_affected_called = my $affecteds_missing = my $unaffecteds_variant = my $affecteds_variant = my $affecteds_ambiguous = 0;
my %genotypes = ();

class Genome::Model::Tools::Analysis::Mendelian::RecessiveAlleles {
	is => 'Command',                       
	
	has => [		
        vcf_file	=> { is => 'Text', doc => "Input in VCF format", is_optional => 0, is_input => 1},
		output_file	=> { is => 'Text', doc => "Output file for VCF of recessive candidate variants", is_optional => 1, is_input => 1},
		output_variants	=> { is => 'Text', doc => "Output file for recessive candidate variants with annotation", is_optional => 1, is_input => 1},
		control_samples	=> { is => 'Text', doc => "Comma-separated list of control sample names", is_optional => 1, is_input => 1},
		male_samples	=> { is => 'Text', doc => "Comma-separated list of male sample names", is_optional => 1, is_input => 1},		
		inheritance_model	=> { is => 'Text', doc => "Presumed model of mendelian inheritance [autosomal-recessive]", is_optional => 1, is_input => 1, default => 'autosomal-recessive'},
		transcript_annotation_file   => { is => 'Text', doc => "WU Transcript annotation in its native format", is_input => 1},
		vep_annotation_file   => { is => 'Text', doc => "VEP annotation in its native format", is_input => 1},
		min_coverage_to_refute	=> { is => 'Text', doc => "Minimum coverage to refute a possible variant in an affected", is_optional => 1, is_input => 1, default => 10},
		max_frequency_to_refute	=> { is => 'Text', doc => "Maximum observed variant allele frequency to refute a possible variant in an affected [5]", is_optional => 1, is_input => 1, default => 10},
		min_affected_variant	=> { is => 'Text', doc => "Minimum number of affecteds with variant to include", is_optional => 1, is_input => 1, default => 1},
		max_unaffected_variant	=> { is => 'Text', doc => "Maximum number of unaffecteds with variant to include", is_optional => 1, is_input => 1, default => 0},
		remove_mendel_failures	=> { is => 'Text', doc => "If set to 1, will put Mendelian inconsistencies into the filtered file", is_optional => 1, is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {    "Filters a VCF for variants that pass Mendelian rules of inheritance"                 
}

sub help_synopsis {
    return <<EOS
This command filters a VCF for variants that pass Mendelian rules of inheritance
EXAMPLE:	gmt analysis mendelian recessive-alleles --vcf-file myVCF.vcf --output-file myVCF.pass.vcf
EOS
}

sub help_detail {
    return <<EOS 

EOS
}

sub execute {
	my $self = shift;
	my $vcf_file = $self->vcf_file;

	my %stats = ();

	my %control_sample = ();
	if($self->control_samples)
	{
		my @samples = split(/\,/, $self->control_samples);
		foreach my $sample (@samples)
		{
			$control_sample{$sample} = 1;
			$stats{'num_samples_control'}++;
		}
	}

	my %male_sample = ();
	if($self->male_samples)
	{
		my @samples = split(/\,/, $self->male_samples);
		foreach my $sample (@samples)
		{
			$male_sample{$sample} = 1;
			$stats{'num_samples_male'}++;
		}
	}

	## Load the annotation ##
	warn "Loading Transcript annotation...\n";
	my %transcript_annotation = load_annotation($self->transcript_annotation_file);
	warn "Loading VEP annotation...\n";
	my %vep_annotation = load_vep($self->vep_annotation_file);
	warn "Loading Gene expression...\n";


	## Open outfile to be able to print the header ##

	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";		
	}
	
	if($self->output_variants)
	{
		open(OUTVARS, ">" . $self->output_variants) or die "Can't open outfile: $!\n";		
	}
		
	

	warn "Parsing VCF...\n";
	

	my %affected_genes = ();
	my %denovo_genes = ();
	my %denovo_variants = ();
	my %gene_variants = ();
	my %control_genes = ();
	my %vcf_lines_by_gene = ();
	my %final_annotations = ();
	
	## Parse the VCF file ##
	
	my $input = new FileHandle ($vcf_file);
	my $lineCounter = 0;
	
	my @sample_names = ();
	my $num_samples = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		warn "$lineCounter lines parsed...\n" if(!($lineCounter % 100000));
		
		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;
		
		if(substr($line, 0, 1) eq '#')
		{
			print OUTFILE "$line\n";
			## Parse out samples ##
			if(substr($line, 0, 6) eq '#CHROM')
			{
				for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
				{
					my $sample = $lineContents[$colCounter];
					if($sample =~ 'Samtools' || $sample =~ 'Varscan')
					{
						## Do nothing ##
					}
					else
					{
						$num_samples++;						
					}

					$sample_names[$colCounter] = $lineContents[$colCounter];
				}
				$stats{'num_samples'} = $num_samples;
			}

		}
		else
		{
			my ($chrom, $position, $id, $ref, $var, $score, $filter, $info, $format) = split(/\t/, $line);
			$stats{'vcf_lines'}++;
			
			## Only take sites passing filters ##
			if($filter eq '.' || uc($filter) eq 'PASS')
			{
				## Only take sites with a variant allele #
				if($ref && $var && $var ne '.')
				{
					$stats{'variants'}++;

					## PART 2: GET THE DBSNP STATUS ##
					
					my $dbsnp_status = "novel";
					my $rs_number = "";
					my $pop_probability = 1.0;
					
					if($id && $id ne ".")
					{
						$stats{'variants_dbsnp_known'}++;
						## We have a dbSNP ##
						$rs_number = $id;
						$dbsnp_status = "known";
						my @infoContents = split(/\;/, $info);
						my %info_values = ();
						foreach my $info_field (@infoContents)
						{
							if($info_field =~ '=')
							{
								my ($name, $value) = split(/\=/, $info_field);
								$info_values{$name} = $value;
							}
							else
							{
								$info_values{$info_field} = 1;
							}
						}
						
						## Common variant is marked as "G5" ##
						
						if($info_values{'G5'} || $info_values{'G5A'})
						{
							## Global MAF of >5% in one or all populations ##
							$dbsnp_status = "common";
						}
						else
						{
							if($info_values{'GMAF'})
							{
								if($info_values{'GMAF'} >= 0.05)
								{
									$dbsnp_status = "common";
								}
								elsif($info_values{'GMAF'} >= 0.01)
								{
									$dbsnp_status = "uncommon";
								}
								else
								{
									$dbsnp_status = "rare";
								}
							}
							else
							{
								## If a very rare mutation ##
								if($info_values{'MUT'} || $info_values{'CLN'})
								{
									$dbsnp_status = "mutation";
								}
							}
						}

					}
					else
					{
						## NO dbSNP ##
						$stats{'variants_dbsnp_novel'}++;
					}

					## PART 2: GET THE MENDEL STATUS ##
					
					my $mendel_status = "FAIL";
					my %variant_case_samples = ();
					my %variant_control_samples = ();
					
					## Get the anticipated format for genotypes ##
					
					my @formatContents = split(/\:/, $format);
					my %genotype_column = ();
					my $numFormatContents = @formatContents;
					for(my $colCounter = 0; $colCounter < $numFormatContents; $colCounter++)
					{
						my $column_name = $formatContents[$colCounter];
						$genotype_column{$column_name} = $colCounter;
					}
					
					my %allele_counts = ();
					my $readable_genotypes = "";
					## Go through and get each sample genotype ##
					
					for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
					{
						my $sample_name = $sample_names[$colCounter];
						
						## Adjust for separate VarScan/SAMtools columns ##
						if($sample_name =~ 'Samtools' || $sample_name =~ 'Varscan')
						{
							## Skip these ##
						}
						else
						{
							my @genotypeContents = split(/\:/, $lineContents[$colCounter]);
							my $genotype = $genotypeContents[$genotype_column{'GT'}];
							my $coverage = $genotypeContents[$genotype_column{'DP'}];
							my $filter = $genotypeContents[$genotype_column{'FT'}];
							my $var_depth = $genotypeContents[$genotype_column{'AD'}];
							my $freq_alt = $genotypeContents[$genotype_column{'FA'}];
	
							## Get the highest var freq ##
							my $var_freq = 0;
							
							if($var_depth && $var_depth ne ".")
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
	#							warn "Trying to convert $genotype in column $colCounter at line $lineCounter\n";
								$genotype = convert_genotype($ref, $var, $genotype) if($genotype ne '.');
	
								my $allele1 = substr($genotype, 0, 1);
								my $allele2 = substr($genotype, 1, 1);
								$allele_counts{$allele1}++;
								$allele_counts{$allele2}++;
		
								my $gt = "Missing";
								
								if($genotype eq $ref . $ref)
								{
									if($var_freq < 0.05)
									{
										$gt = "Ref";
										$stats{'wildtype_without_allele_depth'}++;
									}
									else
									{
										$stats{'wildtype_with_allele_depth'}++;
									}
								}
								elsif($genotype eq $ref . $var)
								{
									$gt = "Het";
								}
								elsif($genotype eq $var . $var)
								{
									$gt = "Hom";
								}
							
								## Only proceed if we have a case sample and they're het/hom ##
								
								if(!$control_sample{$sample_name} && ($gt eq "Het" || $gt eq "Hom"))
								{
									$variant_case_samples{$sample_name} = $gt;
								}
								elsif($control_sample{$sample_name} && ($gt eq "Het" || $gt eq "Hom"))
								{
									$variant_control_samples{$sample_name} = $gt;
								}
													
							}
							else
							{
							}							
						}


						
					}

					
					## PART TWO: GET VARIANT ANNOTATION ##
					my $variant_key = "";
					my $variant_type = "snp";
					if($var =~ ',')
					{
						## If multiple variant alleles, choose most prevalent one ##
						my @vars = split(/\,/, $var);
						my $chosen_var = my $chosen_var_count = 0;
						foreach my $this_var (@vars)
						{
							if(!$chosen_var || $allele_counts{$this_var} && $allele_counts{$this_var} > $chosen_var_count)
							{
								$chosen_var = $this_var;
								$chosen_var_count = $allele_counts{$this_var};
								$chosen_var_count = 0 if(!$chosen_var_count); 
							}
						}
						$variant_key = join("\t", $chrom, $position, $position, $ref, $chosen_var);
						
						## Check and adjust for indel ##
						if(length($ref) > 1 || length($chosen_var) > 1)
						{
							$variant_type = "indel";
							my $chr_start = my $chr_stop = 0;
							my $allele1 = my $allele2 = "";
							if(length($chosen_var) > length($ref))
							{
								$variant_type = "ins";
								## Insertion ##
								$allele2 = $chosen_var;
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
								$allele1 =~ s/$chosen_var//;
								$allele2 = "-";
								my $indel_size = length($allele1);
								$chr_start = $position + 1;
								$chr_stop = $chr_start + $indel_size - 1;
							}
							
							$variant_key = join("\t", $chrom, $chr_start, $chr_stop, $allele1, $allele2);							
						}

					}
					elsif(length($ref) > 1 || length($var) > 1)
					{
						$variant_type = "indel";
						my $chr_start = my $chr_stop = 0;
						my $allele1 = my $allele2 = "";
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
						$variant_key = join("\t", $chrom, $position, $position, $ref, $var);
					}

					$stats{'variants_type_' . $variant_type}++;
					
					if($transcript_annotation{$variant_key})
					{
						$stats{'variants_type_' . $variant_type . '_had_annot_tx'}++;

						my ($tx_gene, $tx_transcript, $tx_trv_type, $tx_c_position, $tx_aa_change, $tx_ucsc_cons, $tx_domain, $tx_errors) = split(/\t/, $transcript_annotation{$variant_key});

						## Get VEP annotation if available ##
						my ($vep_ens_gene, $vep_gene, $vep_class, $vep_cdna_pos, $vep_protein_pos, $vep_amino_acids, $vep_polyphen, $vep_sift, $vep_condel) = split(/\t/, $vep_annotation{$variant_key}) if($vep_annotation{$variant_key});
						$stats{'variants_type_' . $variant_type . '_had_annot_vep'}++ if($vep_annotation{$variant_key});						

						## Determine if we should fail due to an annotation error that VEP does not correct ##
						
						if($tx_errors ne "no_errors" && $tx_errors ne "-" && (!$vep_annotation{$variant_key} || $vep_class eq "INTERGENIC" || convert_vep_class($vep_class, $variant_type) eq $tx_trv_type))
						{
							$stats{'variants_type_' . $variant_type . '_removed_annot_error'}++;
						}
						else
						{
							## If no errors in annotation or VEP rescue, proceed ##
							my $variant_class = my $variant_gene = my $variant_change = my $our_annot = my $vep_annot = "";
							## Get our Annotation ##
							$our_annot = join(",", $tx_gene, $tx_trv_type, $tx_c_position, $tx_aa_change, $tx_ucsc_cons, $tx_errors);

							## Get VEP Annotation ##
							if($vep_class)
							{
								if($vep_gene)
								{
									$vep_annot = join(",", $vep_gene, $vep_class, $vep_cdna_pos, $vep_protein_pos . $vep_amino_acids);									
									$vep_annot .= ",Polyphen:" . $vep_polyphen if($vep_polyphen && length($vep_polyphen) > 1);
									$vep_annot .= ",SIFT:" . $vep_sift if($vep_sift && length($vep_sift) > 1);
									$vep_annot .= ",Condel:" . $vep_condel if($vep_condel && length($vep_condel) > 1);									
								}
								else
								{
									$vep_annot = $vep_class;
								}

							}
							
							## If our annotation has no errors, use it ##
							if($tx_errors eq "no_errors" || $tx_errors eq "-")#convert_vep_class($vep_class, $variant_type) ne $tx_trv_type
							{
								## Use Transcript annotation ##
								$variant_class = $tx_trv_type;
								$variant_gene = $tx_gene;
								$variant_change = join(",", $tx_c_position, $tx_aa_change, $tx_ucsc_cons);
							}
							else
							{
								## Use VEP annotation as tx annotation is bad ##
								$variant_class = convert_vep_class($vep_class, $variant_type);
								$variant_gene = $vep_gene;
							}

							## Proceed with analysis ##

							if($variant_class && $variant_gene && $damaging_classes{$variant_class})
							{
								$stats{'variants_type_' . $variant_type . '_with_annot'}++;
								
								## Require rareness ##
								
								if($dbsnp_status ne "common" && $dbsnp_status ne "uncommon")
								{
									$vcf_lines_by_gene{$variant_gene} .= "\n" if($vcf_lines_by_gene{$variant_gene});
									$vcf_lines_by_gene{$variant_gene} .= $line;
	
									## First go through any controls with the variant. If one is homozygous, rule it out ##
									my $is_recessive = 1;
									my $is_denovo = 1;
	
									foreach my $sample_name (sort keys %variant_control_samples)
									{
										my $gt = $variant_control_samples{$sample_name};
										my $key = join("\t", $variant_gene, $sample_name);
										$is_denovo = 0 if($gt eq "Het" || $gt eq "Hom");
	
										if($gt eq "Hom")
										{
											$is_recessive = 0;
											$gene_variants{$variant_gene} .= "\n" if($gene_variants{$variant_gene});
											$gene_variants{$variant_gene} .= join("\t", $sample_name, $gt, $chrom, $position);										
										}
										elsif($is_recessive)
										{
											$gene_variants{$variant_gene} .= "\n" if($gene_variants{$variant_gene});
											$gene_variants{$variant_gene} .= join("\t", $sample_name, $gt, $chrom, $position);										
										}
	
									}
	
									## If still possibly recessive and found in an affected, save this gene and affected genotypes ##
									if($is_recessive)
									{
										## Save the result for this one ##
										foreach my $sample_name (sort keys %variant_case_samples)
										{
											my $gt = $variant_case_samples{$sample_name};
											my $key = join("\t", $variant_gene, $sample_name);
											$affected_genes{$variant_gene} = 1;
											$gene_variants{$variant_gene} .= "\n" if($gene_variants{$variant_gene});
											$gene_variants{$variant_gene} .= join("\t", $sample_name, $gt, $chrom, $position);
										}
										
										$final_annotations{$chrom . "\t" . $position} = join("\t", $variant_class, $variant_change, $vep_annot, $dbsnp_status);
									}
									
									if($is_denovo)
									{
										my $num_cases_het = 0;
										foreach my $sample_name (sort keys %variant_case_samples)
										{
											my $gt = $variant_case_samples{$sample_name};
											if($gt eq "Het")
											{
												$num_cases_het++;
											}
										}
										
										if($num_cases_het > 0)
										{
											$affected_genes{$variant_gene} = 1;
											$denovo_genes{$variant_gene} = 1;
											my $variant_key = join("\t", $chrom, $position);
											$denovo_variants{$variant_key} = 1;
											$final_annotations{$chrom . "\t" . $position} = join("\t", $variant_class, $variant_change, $vep_annot, $dbsnp_status);
											
										}
										else
										{
											$is_denovo = 0;
										}
									}
									
								}


								
								
								
							}
							else
							{
								$stats{'variants_type_' . $variant_type . '_removed_annot_unknown'}++;
							}

						}

					}
					else
					{
						$stats{'variants_type_' . $variant_type . '_removed_annot_missing'}++;
					}
				}
				else
				{
					$stats{'non_variant_lines'}++;
				}
			}
			else
			{
				$stats{'lines_filter_' . $filter}++;							
			}
			
			
		} ## End If-else for header ##
		

	} ## Finish parsing file ##
	
	close($input);


	## Determine recessive genes ##



	my %recessive_genes = ();
	my @vcfLines = ();
	my $vcfLineCounter = 0;

	foreach my $gene (sort keys %affected_genes)
	{
		$stats{'num_genes_variant'}++;
		
		my %het_controls = my %hom_controls = ();
		my %het_cases = my %hom_cases = ();
		my %num_case_hets = ();
		
		## Go through all variants saved for this gene ##
		
		my $gene_is_recessive = 0;
		my $control_was_homozygous = 0;
		my %recessive_variants = ();
		
		my $gene_is_denovo = 0;
		$gene_is_denovo = 1 if($denovo_genes{$gene});
		
		my $controls_het = my $controls_hom = my $cases_het = my $cases_hom = 0;
		
		my @results = split(/\n/, $gene_variants{$gene});
		foreach my $result (@results)
		{
			my ($sample_name, $gt, $chrom, $position) = split(/\t/, $result);
			my $variant_key = join("\t", $chrom, $position);

			## Sample was a control ##
			if($control_sample{$sample_name})
			{
				if($gt eq "Het")
				{
					$controls_het++;
					$het_controls{$variant_key} .= "\n" if($het_controls{$variant_key});
					$het_controls{$variant_key} .= $sample_name;
				}
				elsif($gt eq "Hom")
				{
					$controls_hom++;
					$hom_controls{$variant_key} .= "\n" if($hom_controls{$variant_key});
					$hom_controls{$variant_key} .= $sample_name;
					$control_was_homozygous = 1;
				}
			}
			else
			{
				## An affected case, so save the status of this variant ##
				if($gt eq "Het")
				{
					$cases_het++;
					$het_cases{$variant_key} .= "\n" if($het_cases{$variant_key});
					$het_cases{$variant_key} .= $sample_name;
				}
				elsif($gt eq "Hom")
				{
					$cases_hom++;
					$hom_cases{$variant_key} .= "\n" if($hom_cases{$variant_key});
					$hom_cases{$variant_key} .= $sample_name;
				}
			}


		}
		
#		print join("\t", $controls_het, $controls_hom, $cases_het, $cases_hom) . "\n";
		
		## Go through any variants homozygous in an affected ##
		foreach my $variant_key (keys %hom_cases)
		{
			if(!$hom_controls{$variant_key})
			{
				$gene_is_recessive = 1;				
				$recessive_variants{$variant_key} = 1;
			}
		}
		
		my %shared_hets = ();
		
		## Go through any variants heterozygous in an affected ##
		foreach my $variant_key (keys %het_cases)
		{
			my @hetCases = split(/\n/, $het_cases{$variant_key});
			if($het_controls{$variant_key})
			{
				my @hetControls = split(/\n/, $het_controls{$variant_key});
	
				## Save all instances of hets shared by an affected and an unaffected parent ##
	
				foreach my $case_sample (@hetCases)
				{
					foreach my $control_sample(@hetControls)
					{
						$recessive_variants{$variant_key} = 1;
						my $pair_key = join("\t", $case_sample, $control_sample);
						$shared_hets{$pair_key}++;
					}
				}				
			}

		}

		## Count the number of shared hets inherited from unique parents ##

		my $num_shared_hets = 0;
		foreach my $pair_key (keys %shared_hets)
		{
			$num_shared_hets++;
		}
		
#		print join("\t", $gene, $num_shared_hets) . "\n";
		
		## If we got at least one het from each parent, this fits compound het. Note this may need adjustment if multiple affecteds ##
		
		my $recessive_type = "unknown";
		if($num_shared_hets >= 2 && !$control_was_homozygous)
		{
			$gene_is_recessive = 1;
			$stats{'num_genes_variant_recessive_compound_het'}++;
			$recessive_type = "compound-het";
#			print join("\t", $gene, "CompoundHet") . "\n";
		}
		elsif($gene_is_recessive && !$control_was_homozygous)
		{
			$stats{'num_genes_variant_recessive_rare_homozygote'}++;
#			print join("\t", $gene, "RareHom") . "\n";
			$recessive_type = "homozygous";
		}
		
		if($control_was_homozygous)
		{
			$stats{'num_genes_variant_control_homozygote'}++;
			$gene_is_recessive = 0;
		}

		
		
		if($gene_is_recessive)
		{
			$stats{'num_genes_variant_recessive'}++;
			$recessive_genes{$gene} = 1;

			my @geneLines = split(/\n/, $vcf_lines_by_gene{$gene});
			foreach my $gene_line (@geneLines)
			{
				my ($chrom, $position, $id, $ref, $var, $score, $filter, $info, $format) = split(/\t/, $gene_line);
				my @lineContents = split(/\t/, $gene_line);
				my $numContents = @lineContents;
				
				my $genotype_string = "";
				for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
				{
					$genotype_string .= "\t" if($genotype_string);
					$genotype_string .= $lineContents[$colCounter];
				}
				
				my $variant_key = join("\t", $chrom, $position);
				if($recessive_variants{$variant_key})
				{
					
					print OUTVARS join("\t", $chrom, $position, $ref, $var, $gene, $recessive_type, $final_annotations{$chrom . "\t" . $position}, $genotype_string) . "\n";
									
					
					$stats{'num_recessive_variants'}++;
					$vcfLines[$vcfLineCounter] = $gene_line;
					$vcfLineCounter++;					
				}

			}

		}
		
		if($gene_is_denovo)
		{
			$stats{'num_genes_variant_denovo'}++;

			my @geneLines = split(/\n/, $vcf_lines_by_gene{$gene});
			foreach my $gene_line (@geneLines)
			{
				my ($chrom, $position, $id, $ref, $var, $score, $filter, $info, $format) = split(/\t/, $gene_line);
				my @lineContents = split(/\t/, $gene_line);
				my $numContents = @lineContents;
				
				my $genotype_string = "";
				for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
				{
					$genotype_string .= "\t" if($genotype_string);
					$genotype_string .= $lineContents[$colCounter];
				}
				
				my $variant_key = join("\t", $chrom, $position);
				if($denovo_variants{$variant_key})
				{
					print OUTVARS join("\t", $chrom, $position, $ref, $var, $gene, "denovo", $final_annotations{$chrom . "\t" . $position}, $genotype_string) . "\n";				
					$stats{'num_denovo_variants'}++;
					$vcfLines[$vcfLineCounter] = $gene_line;
					$vcfLineCounter++;					
				}

			}			
		}
	}

	## Sort all VCF lines by chrom/position ##
	@vcfLines = sort byChrPos (@vcfLines);
	foreach my $line (@vcfLines)
	{
		print OUTFILE "$line\n";
	}

	if($self->output_file)
	{
		close(OUTFILE);
	}

	foreach my $key (sort keys %stats)
	{
		print "$stats{$key}\t$key\n";
	}
	
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

sub load_annotation
{
    my $annotation_file = shift(@_);
    my %annotation = ();
    
    my $input = new FileHandle ($annotation_file);
    my $lineCounter = 0;

    while (<$input>)
    {
        chomp;
        my $line = $_;
        
        my ($chrom, $chr_start, $chr_stop, $ref, $var, $vartype, $gene, $transcript, $species, $source, $version, $strand, $status, $trv_type, $c_position, $aa_change, $ucsc_cons, $domain, $all_domains, $del_subs, $errors) = split(/\t/, $line);

	my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
	$gene = "-" if(!$gene);
	$transcript = "-" if(!$transcript);
	$c_position = "-" if(!$c_position);
	$aa_change = "-" if(!$aa_change);
	$domain = "-" if(!$domain);
	$ucsc_cons = "-" if(!$ucsc_cons);
	$annotation{$key} = join("\t", $gene, $transcript, $trv_type, $c_position, $aa_change, $ucsc_cons, $domain, $errors);

    }
    
    close($input);

    return(%annotation);
}

sub numericallyDesc
{
	$a =~ s/\./0/;
	$b =~ s/\./0/;
	$b <=> $a;
}

sub load_vep
{
    my $annotation_file = shift(@_);
    my %annotation = ();
    
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

            my @classes = split(/\,/, $class);
            foreach my $class (@classes)
            {
                    if($polyphen || $sift || $condel)
                    {
                            $annotation{$key} .= "\n" if($annotation{$key});
                            $annotation{$key} .= join("\t", $ens_gene, $gene, $class, $cdna_pos, $protein_pos, $amino_acids, $polyphen, $sift, $condel)
#				print join("\t", $chrom, $position, $alleles, $ens_gene, $gene, $class, $protein_pos, $amino_acids, $polyphen, $sift, $condel) . "\n";
                    }
                    else
                    {
                            $annotation{$key} .= "\n" if($annotation{$key});
                            $annotation{$key} .= join("\t", $ens_gene, $gene, $class, $cdna_pos, "", "", "", "", "")				;
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
	my ($ens_gene_a, $gene_a, $class_a, $cdna_pos_a, $protein_pos_a, $amino_acids_a, $polyphen_a, $sift_a, $condel_a) = split(/\t/, $a);

	$polyphen_a = 0 if(!$polyphen_a);
	$sift_a = 0 if(!$sift_a);
	$condel_a = 0 if(!$condel_a);
	
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
	
	my ($ens_gene_b, $gene_b, $class_b, $cdna_pos_b, $protein_pos_b, $amino_acids_b, $polyphen_b, $sift_b, $condel_b) = split(/\t/, $b);

	$polyphen_b = 0 if(!$polyphen_b);
	$sift_b = 0 if(!$sift_b);
	$condel_b = 0 if(!$condel_b);

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

sub convert_vep_class
{
	my $vep_class = shift(@_);
	my $variant_type = shift(@_);
	
	if($variant_type eq "ins" || $variant_type eq "del")
	{
		return("in_frame_" . $variant_type) if($vep_class eq "NON_SYNONYMOUS_CODING");
		return("frame_shift_" . $variant_type) if($vep_class eq "FRAMESHIFT_CODING");
	}
	
	
	return("rna") if($vep_class eq 'NMD_TRANSCRIPT');
	return("nonstop") if($vep_class eq 'PARTIAL_CODON');
	return("-") if($vep_class eq 'INTERGENIC');
	return("5_prime_flanking_region") if($vep_class eq 'UPSTREAM');
	return("3_prime_flanking_region") if($vep_class eq 'DOWNSTREAM');
	return("intronic") if($vep_class eq 'INTRONIC');
	return("5_prime_untranslated_region") if($vep_class eq '5PRIME_UTR');
	return("3_prime_untranslated_region") if($vep_class eq '3PRIME_UTR');
	return("rna") if($vep_class eq 'WITHIN_NON_CODING_GENE');
	return("rna") if($vep_class eq 'WITHIN_MATURE_miRNA');
	return("silent") if($vep_class eq 'CODING_UNKNOWN');
	return("silent") if($vep_class eq 'SYNONYMOUS_CODING');
	return("nonstop") if($vep_class eq 'STOP_LOST');
	return("splice_region") if($vep_class eq 'SPLICE_SITE');
	return("splice_site") if($vep_class eq 'ESSENTIAL_SPLICE_SITE');
	return("missense") if($vep_class eq 'NON_SYNONYMOUS_CODING');
	return("nonsense") if($vep_class eq 'STOP_GAINED');
	return("-") if($vep_class eq '-');
	warn "Unable to convert vep class $vep_class\n";
	return("-");
}


sub byCode
{
	my $code_a = my $code_b = 0;

	foreach my $test_class (keys %vep_class_rank)
	{
		$code_a = $vep_class_rank{$test_class} if($a eq $test_class);
		$code_b = $vep_class_rank{$test_class} if($b eq $test_class);
	}
	
	$code_b <=> $code_a;
}

sub get_code
{
	my $class = shift(@_);
	
	my @classes = split(/\,/, $class);
	my $num_classes = @classes;
	
	if($num_classes > 1)
	{
		warn "Warning: multiple VEP classes sent for VEP code: $class\n";
	}
	
	foreach my $test_class (keys %vep_class_rank)
	{
		return($vep_class_rank{$test_class}) if($class eq $test_class);
	}
	
	return(0);
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

1;
