
package Genome::Model::Tools::Analysis::Mendelian::PrioritizeVcf;     # rename this when you give the module file a different name <--

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

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

my $num_affected_called = my $affecteds_missing = my $unaffecteds_variant = my $affecteds_variant = my $affecteds_ambiguous = 0;
my %genotypes = ();

class Genome::Model::Tools::Analysis::Mendelian::PrioritizeVcf {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		vcf_file	=> { is => 'Text', doc => "Input in VCF format", is_optional => 0, is_input => 1},
		output_file	=> { is => 'Text', doc => "Output file for pass-filter sites", is_optional => 1, is_input => 1},
		filtered_file	=> { is => 'Text', doc => "Output file for fail-filter sites", is_optional => 1, is_input => 1},
		readable_file	=> { is => 'Text', doc => "A simplified tab-delimited readable file for collaborators", is_optional => 1, is_input => 1},
		control_samples	=> { is => 'Text', doc => "Comma-separated list of control sample names", is_optional => 1, is_input => 1},
		carrier_samples	=> { is => 'Text', doc => "Comma-separated list of carrier sample names (e.g. parents in auto recessive)", is_optional => 1, is_input => 1},
		male_samples	=> { is => 'Text', doc => "Comma-separated list of male sample names", is_optional => 1, is_input => 1},		
		inheritance_model	=> { is => 'Text', doc => "Presumed model of mendelian inheritance [autosomal-dominant]", is_optional => 1, is_input => 1, default => 'autosomal-dominant'},
		transcript_annotation_file   => { is => 'Text', doc => "WU Transcript annotation in its native format", is_input => 1},
		vep_annotation_file   => { is => 'Text', doc => "VEP annotation in its native format", is_input => 1},
		gene_expression_file   => { is => 'Text', doc => "Two-column file of gene and expression FPKM value", is_input => 1, is_optional => 1},
		shared_ibd_file   => { is => 'Text', doc => "ibdregions file from BEAGLE, chrom start stop markers samples", is_input => 1, is_optional => 1},
		priority_gene_list   => { is => 'Text', doc => "List of gene symbols to flag in scored output", is_input => 1, is_optional => 1},
		min_coverage_to_refute	=> { is => 'Text', doc => "Minimum coverage to refute a possible variant in an affected", is_optional => 1, is_input => 1, default => 10},
		max_frequency_to_refute	=> { is => 'Text', doc => "Maximum observed variant allele frequency to refute a possible variant in an affected [5]", is_optional => 1, is_input => 1, default => 10},
		min_affected_variant	=> { is => 'Text', doc => "Minimum number of affecteds with variant to include", is_optional => 1, is_input => 1, default => 1},
		max_unaffected_variant	=> { is => 'Text', doc => "Maximum number of unaffecteds with variant to include", is_optional => 1, is_input => 1, default => 0},
		remove_mendel_failures	=> { is => 'Text', doc => "If set to 1, will put Mendelian inconsistencies into the filtered file", is_optional => 1, is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Filters a VCF for variants that pass Mendelian rules of inheritance"                 
}

sub help_synopsis {
    return <<EOS
This command filters a VCF for variants that pass Mendelian rules of inheritance
EXAMPLE:	gmt analysis mendelian filter-vcf --vcf-file myVCF.vcf --output-file myVCF.pass.vcf --filtered-file myVCF.fail.vcf
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

	my %carrier_sample = ();
	if($self->carrier_samples)
	{
		my @samples = split(/\,/, $self->carrier_samples);
		foreach my $sample (@samples)
		{
			$carrier_sample{$sample} = 1;
			$stats{'num_samples_carrier'}++;
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


	my %gene_expression = load_expression($self->gene_expression_file) if($self->gene_expression_file);

	my %priority_genes = load_priority_genes($self->priority_gene_list) if($self->priority_gene_list);	
	
	my %shared_ibd = load_shared_ibd($self->shared_ibd_file) if($self->shared_ibd_file);

	warn "Parsing VCF...\n";
	
	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";		
	}
	
	if($self->filtered_file)
	{
		open(FILTERFILE, ">" . $self->filtered_file) or die "Can't open filtered outfile: $!\n";				
	}
	
	if($self->readable_file)
	{
		open(READABLEFILE, ">" . $self->readable_file) or die "Can't open readable outfile: $!\n";				
	}
	

	my @candidate_variants = ();
	my $num_candidates = 0;
	my %candidate_variants_readable = ();
	
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
		
		warn "$lineCounter lines parsed...\n" if(!($lineCounter % 20000));
		
		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;
		
		if(substr($line, 0, 1) eq '#')
		{
			if(substr($line, 0, 6) eq '#CHROM')
			{
				## Print our additional filter/info fields ##
				print OUTFILE qq{##INFO=<ID=RPscore,Number=A,Type=Float,Description="Combined RP probability score (higher=better): MendelScore * PopScore * EffectScore * ExprScore">\n};
				print OUTFILE qq{##INFO=<ID=Annot,Number=1,Type=String,Description="Our internal annotation with gene, trv_type, codon, aa_change, conservation, errors">\n};
				print OUTFILE qq{##INFO=<ID=VEP,Number=1,Type=String,Description="VEP annotation with gene, trv_type, codon, aa_change, Polyphen/SIFT/Condel">\n};

				print OUTFILE qq{##INFO=<ID=MendelScore,Number=A,Type=Float,Description="Inheritance score. 0.20 if affected(s) wildtype, 0.05 if control(s) had variant. Otherwise 1 - (affecteds homozygous * 0.20) or 1 - (affecteds not called * 0.10)">\n};
				print OUTFILE qq{##INFO=<ID=PopScore,Number=A,Type=Float,Description="Rareness score. 1 if novel to dbSNP. 0.95 if dbSNP mutation. 0.60 if dbSNP but no MAFs. Otherwise, 0.20 if MAF<1%, 0.02 if 1%<MAF<5%, or 0.001 if MAF>5%.">\n};
				print OUTFILE qq{##INFO=<ID=EffectScore,Number=A,Type=Float,Description="Predicted deleteriousness score. 1 if truncating, 0.95 if damaging missense, 0.80 if missense, 0.20 if splice_region, 0.05 if synonymous or inframe">\n};
				print OUTFILE qq{##INFO=<ID=ExprScore,Number=A,Type=Float,Description="Gene expression rank score. 1 if in top 25%, 0.95 if top 50%, 0.80 if bottom 25%. 0.90 if rank 50-75% or not available">\n};			

				print OUTFILE qq{##FILTER=<ID=AnnotationError,Description="Annotation could be erroneous (pseudogene, no stop codon, etc.)">\n};
				print OUTFILE qq{##FILTER=<ID=NoAnnotation,Description="No transcript/VEP annotation found for variant. This shouldn't happen.">\n};
				print OUTFILE qq{##FILTER=<ID=NotVariant,Description="No variant alleles observed. This shouldn't happen.">\n};
				print OUTFILE qq{##FILTER=<ID=NoVariantClass,Description="Unable to determine variant class from annotation. This shouldn't happen.">\n};


				## Print our additional filter/info fields ##
				print FILTERFILE qq{##INFO=<ID=RPscore,Number=A,Type=Float,Description="Combined RP probability score (higher=better): MendelScore * PopScore * EffectScore * ExprScore">\n};
				print FILTERFILE qq{##INFO=<ID=Annot,Number=1,Type=String,Description="Our internal annotation with gene, trv_type, codon, aa_change, conservation, errors">\n};
				print FILTERFILE qq{##INFO=<ID=VEP,Number=1,Type=String,Description="VEP annotation with gene, trv_type, codon, aa_change, Polyphen/SIFT/Condel">\n};

				print FILTERFILE qq{##INFO=<ID=MendelScore,Number=A,Type=Float,Description="Inheritance score. 0.20 if affected(s) wildtype, 0.05 if control(s) had variant. Otherwise 1 - (affecteds homozygous * 0.20) or 1 - (affecteds not called * 0.10)">\n};
				print FILTERFILE qq{##INFO=<ID=PopScore,Number=A,Type=Float,Description="Rareness score. 1 if novel to dbSNP. 0.95 if dbSNP mutation. 0.60 if dbSNP but no MAFs. Otherwise, 0.20 if MAF<1%, 0.02 if 1%<MAF<5%, or 0.001 if MAF>5%.">\n};
				print FILTERFILE qq{##INFO=<ID=EffectScore,Number=A,Type=Float,Description="Predicted deleteriousness score. 1 if truncating, 0.95 if damaging missense, 0.80 if missense, 0.20 if splice_region, 0.05 if synonymous or inframe">\n};
				print FILTERFILE qq{##INFO=<ID=ExprScore,Number=A,Type=Float,Description="Gene expression rank score. 1 if in top 25%, 0.95 if top 50%, 0.80 if bottom 25%. 0.90 if rank 50-75% or not available">\n};			

				print FILTERFILE qq{##FILTER=<ID=AnnotationError,Description="Annotation could be erroneous (pseudogene, no stop codon, etc.)">\n};
				print FILTERFILE qq{##FILTER=<ID=NoAnnotation,Description="No transcript/VEP annotation found for variant. This shouldn't happen.">\n};
				print FILTERFILE qq{##FILTER=<ID=NotVariant,Description="No variant alleles observed. This shouldn't happen.">\n};
				print FILTERFILE qq{##FILTER=<ID=NoVariantClass,Description="Unable to determine variant class from annotation. This shouldn't happen.">\n};

				if($self->readable_file)
				{
					print READABLEFILE "#CHROM\tPOS\tREF\tALT\tPRIORITY_FLAG\tRP_SCORE\tMENDEL_SCORE\tPOP_SCORE\tEFFECT_SCORE\tEXPR_SCORE\tIBD_SCORE\t";
					print READABLEFILE "MENDEL_STATUS\tDBSNP_STATUS\tDBSNP_ID\tDBSNP_INFO\t";
					print READABLEFILE "VARIANT_CLASS\tVARIANT_GENE\tWU_ANNOTATION\tVEP_ANNOTATION";
				}

				for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
				{
					my $sample = $lineContents[$colCounter];
					if($sample =~ 'Samtools' || $sample =~ 'Varscan')
					{
						## Do nothing ##
					}
					else
					{
						print READABLEFILE "\t$sample-genotype\t$sample-depth\t$sample-varfreq"  if($self->readable_file);
						$num_samples++;						
					}

					$sample_names[$colCounter] = $lineContents[$colCounter];
				}
				print READABLEFILE "\n" if($self->readable_file);

				$stats{'num_samples'} = $num_samples;
			}

			print OUTFILE "$line\n" if($self->output_file);
			print FILTERFILE "$line\n" if($self->filtered_file);
			

		}
		else
		{
			my ($chrom, $position, $id, $ref, $var, $score, $filter, $info, $format) = split(/\t/, $line);
			$stats{'vcf_lines'}++;
			
			## Only take sites passing filters ##
			if($filter eq '.' || uc($filter) eq 'PASS')
			{
#				$stats{'lines_in_vcf'}++;
				
				## Only take sites with a variant allele #
				if($ref && $var && $var ne '.')
				{
					$stats{'variants'}++;
					
					## Get the dbSNP Status ##
					
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
								if($info_values{'MUT'} || $info_values{'CLN'} || $info_values{'PM'})
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


					## Determine population probability ##

					if($dbsnp_status eq "common")
					{
						## common means it is almost definitely not causal ##
						$pop_probability = 0.001;
					}
					elsif($dbsnp_status eq "uncommon")
					{
						## Freq of 1-5% means virtually no chance of causality ##
						$pop_probability = 0.02;
					}
					elsif($dbsnp_status eq "rare")
					{
						## Verifiably rare, but seen in populations, so about 3.23% chance of being causale ##
						$pop_probability = 0.20;							
					}
					elsif($dbsnp_status eq "known")
					{
						## About 1/3 chance that it was pulled in from OMIM ##
						$pop_probability = 0.60;							
					}
					elsif($dbsnp_status eq "mutation")
					{
						## Slight penalty since we'd likely not be screening sample ##
						$pop_probability = 0.95;	
					}
					elsif($dbsnp_status eq "novel")
					{
						$pop_probability = 1.0;
					}
					else
					{
						warn "Warning: dbsnp status $dbsnp_status not scored!\n";
					}

					
					
					## PART 1: GET THE MENDEL STATUS ##
					
					my $mendel_status = "PASS";
					my $total_affected = my $total_control = 0;
					my $num_affected_called = my $num_control_called = my $num_control_called_variant = my $num_affected_called_variant = 0;
					my $num_affected_called_hom = my $num_affected_called_wildtype = my $num_male_affected_X_het = 0;
					my $num_carrier_called = my $num_carrier_called_het = my $num_carrier_called_hom = my $num_carrier_called_wildtype = 0;
					my $num_carrier_called_lowfreq = my $num_control_called_lowfreq = my $num_affected_called_lowfreq = 0;
					
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
							my $filter = ".";
							my $freq_alt = ".";
							$filter = $genotypeContents[$genotype_column{'FT'}] if($genotype_column{'FT'});
							my $var_depth = $genotypeContents[$genotype_column{'AD'}];
							
							if($genotype_column{'FA'})
							{
								$freq_alt = $genotypeContents[$genotype_column{'FA'}];								
							}
							elsif($genotype_column{'FREQ'})
							{
								$freq_alt = $genotypeContents[$genotype_column{'FREQ'}];																
							}

	
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
	
							if($control_sample{$sample_name})
							{
								$total_control++;
							}
							else
							{
								$total_affected++;	
							}
							
							## Only process the genotype if it has a value and is either unfiltered or for the control sample ##
							if(length($genotype) > 2 && $genotype ne '.' && $genotype ne "./." && ($filter eq 'PASS' || $filter eq '.' || $control_sample{$sample_name}))
							{
							
	#							warn "Trying to convert $genotype in column $colCounter at line $lineCounter\n";
								$genotype = convert_genotype($ref, $var, $genotype) if($genotype ne '.');
	
								$readable_genotypes .= "\t" if($readable_genotypes);
								$readable_genotypes .= join("\t", $genotype, $coverage, $freq_alt);
	
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
										$gt = "LowFreq";
										$stats{'wildtype_with_allele_depth'}++;
									}
								}
								elsif($var =~ ",")
								{
									my @vars = split(/\,/, $var);
									foreach my $this_var (@vars)
									{
										if($genotype eq $ref . $ref)
										{
											if($var_freq < 0.05)
											{
												$gt = "Ref";
												$stats{'wildtype_without_allele_depth'}++;
											}
											else
											{
												$gt = "LowFreq";
												$stats{'wildtype_with_allele_depth'}++;
											}
										}
										elsif($genotype eq $ref . $this_var)
										{
											$gt = "Het";
										}
										elsif($genotype eq $this_var . $this_var)
										{
											$gt = "Hom";
										}
										elsif($genotype =~ '\/')
										{
											my ($a1, $a2) = split(/\//, $genotype);
											if($a1 eq $a2)
											{
												$gt = "Hom";
											}
											elsif($a1 eq "-" || $a2 eq "-")
											{
												$gt = "Het";
											}
										}

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
								elsif($genotype =~ '\/')
								{
									my ($a1, $a2) = split(/\//, $genotype);
									if($a1 eq $a2)
									{
										$gt = "Hom";
									}
									elsif($a1 eq "-" || $a2 eq "-")
									{
										$gt = "Het";
									}
									else
									{
										warn "Genotype $genotype could not be matched to alleles $ref/$var 3\n";
									}
								}
								else
								{						
									warn "Genotype $genotype could not be matched to alleles $ref/$var 4\n";
								}
	
							
								if($gt ne "Missing")
								{
									if($self->inheritance_model eq "autosomal-dominant")							
									{
										if($control_sample{$sample_name})
										{
											## Control sample ##
											$num_control_called++;
											
											if($gt eq 'Het' || $gt eq 'Hom')
											{
												$num_control_called_variant++;
												$mendel_status = "Control_Was_Variant";												
											}
											elsif($gt eq "Ref")
											{
#												$mendel_status = "Control_Was_Wildtype";	
											}
										}
										else
										{
											## Affected sample ##
											$num_affected_called++;
											
											if($gt eq 'Ref' && $coverage >= $self->min_coverage_to_refute)
											{
												$num_affected_called_wildtype++;
												$mendel_status = "Affected_Was_Wildtype";
											}
											elsif($gt eq 'Het' || $gt eq 'Hom')
											{
												## Pass ##
												$num_affected_called_variant++;
												if($gt eq 'Hom')
												{
													if($male_sample{$sample_name} && ($chrom eq "X" || $chrom eq "Y"))
													{
														## No penalty for males on sex chromosomes
													}
													else
													{
														$num_affected_called_hom++;												
													}
												}
												else
												{
													## Penalize male samples that are heterozygous ##
													if($male_sample{$sample_name} && $chrom eq "X")
													{
														$num_male_affected_X_het++;
													}
													## Het variant in auto-recessive ##
													if($self->inheritance_model eq "autosomal-recessive" && $gt eq "Het")
													{
														$mendel_status = "Affected_Was_Heterozygous";
													}
												}
		
											
											}
										}										
									}
									elsif($self->inheritance_model eq "autosomal-recessive")
									{
										if($carrier_sample{$sample_name})
										{
											## Control sample ##
											$num_carrier_called++;											
											if($gt eq 'Ref' && $coverage >= $self->min_coverage_to_refute)
											{
												$num_carrier_called_wildtype++;
												$mendel_status = "Carrier_Was_Wildtype";
											}
											elsif($gt eq "Het")
											{
												$num_carrier_called_het++;
#												$mendel_status = "PASS";
											}
											elsif($gt eq "LowFreq")
											{
												$num_carrier_called_lowfreq++;
											}
											elsif($gt eq "Hom")
											{
												$num_carrier_called_hom++;
												$mendel_status = "Carrier_Was_Homozygous";
											}											
										}										
										elsif($control_sample{$sample_name})
										{
											## Control sample ##
											$num_control_called++;											
											if($gt eq 'Ref' && $coverage >= $self->min_coverage_to_refute)
											{
#												$mendel_status = "PASS";
											}
											elsif($gt eq "Het")
											{
												#$mendel_status = "Control_Was_Heterozygous";
											}
											elsif($gt eq "Hom")
											{
												$mendel_status = "Control_Was_Homozygous";
											}
											elsif($gt eq "LowFreq")
											{
												$num_control_called_lowfreq++;
											}											
										}
										else
										{
											## Affected sample ##
											$num_affected_called++;
											
											if($gt eq 'Ref' && $coverage >= $self->min_coverage_to_refute)
											{
												$num_affected_called_wildtype++;
												$mendel_status = "Affected_Was_Wildtype";
											}
											elsif($gt eq 'Ref')
											{
												## Called WT but at lower coverage ##
												$mendel_status = "Affected_Not_Called";
											}
											elsif($gt eq "Het")
											{
												$mendel_status = "Affected_Was_Het";
											}
											elsif($gt eq "Hom")
											{
												$num_affected_called_hom++;
#												$mendel_status = "PASS";
											}
											elsif($gt eq "LowFreq")
											{
												$num_affected_called_lowfreq++;
											}
										}
										
										

								
										
									}


								}
								else
								{
									## Genotype is missing ##
								}
								

							}
							else
							{
#								print "Not converting because gt=$genotype and filter=$filter\n";
								## Genotype missing or was filtered out ####
								$genotype = "NN";
								$readable_genotypes .= "\t" if($readable_genotypes);
								$coverage = "-" if(!$coverage);
								$freq_alt = "-" if(!$freq_alt);
								$readable_genotypes .= join("\t", $genotype, $coverage, $freq_alt);

							}															



						}


						
					}


					my $mendel_probability = 0;
					
					## Determine probability based on mendel segregation status, where any error must be a wrong variant call ##
					if($self->inheritance_model eq "autosomal-dominant")
					{
						if($mendel_status eq "PASS" || $mendel_status eq "Control_Was_Wildtype")
						{
							if($num_affected_called_hom > 0 || $num_male_affected_X_het > 0)
							{
								## Homozygous in an affected. This is unexpected but not a deal breaker ##
								if($self->inheritance_model eq "autosomal-dominant")
								{
									$mendel_probability = 1 - ($num_affected_called_hom * 0.20);								
								}
	
								## Heterozygous in a male on chromosome X. This indicates an ugly variant call ##
								$mendel_probability -= ($num_male_affected_X_het * 0.20);
								$mendel_probability = 0.20 if($mendel_probability < 0.20);
							}
							elsif($num_affected_called_variant == $total_affected)
							{
								## Present in every affected - no penalty ##
								$mendel_probability = 1.0;
							}
							else #if($num_affected_called_variant >= $self->min_affected_variant)
							{
								my $samples_missing = $total_affected - $num_affected_called_variant;
								## Present in perhaps all but one affected. Still promising, but less so. For each sample where missed, reduce mendel score expecting 90% coverage ##
								$mendel_probability = 1 - ($samples_missing * 0.10);
								$mendel_probability = 0.20 if($mendel_probability < 0.20);
							}
						}
						elsif($mendel_status eq "Affected_Was_Wildtype")
						{
							## Here we'd have To have missed a variant. It does happen. Assume 80% sensitivity
							## Another possibility is that there are two causal variants in a family pedigree. 
							$mendel_probability = 0.20;
						}
						elsif($mendel_status eq "Control_Was_Variant")
						{
							## Here we'd have to have miscalled a variant in a control. Assume that happens rarely ##
							$mendel_probability = 0.05;
						}
						else
						{
							warn "Mendel probability unknown for $mendel_status\n";
						}						
					}
					elsif($self->inheritance_model eq "autosomal-recessive")
					{
						if($mendel_status eq "PASS")
						{
							$mendel_probability = 1.0;
							## Ideally, we'd like all carriers to be hets and all affecteds to be homozygous ##

							## Penalty if affecteds lowfreq ##
							
							if($num_affected_called_lowfreq)
							{
								$mendel_probability = $mendel_probability - ($num_affected_called_lowfreq * 0.40);
							}

							if($num_control_called_lowfreq)
							{
								$mendel_probability = $mendel_probability - ($num_control_called_lowfreq * 0.20);
							}
							
							if($num_carrier_called_lowfreq)
							{
								$mendel_probability = $mendel_probability - ($num_carrier_called_lowfreq * 0.10);
							}

							## Penalty if not all carriers were het ##
							if($num_carrier_called_het < $num_carrier_called)
							{
								## If this occurs, then we called the carrier(s) wildtype but they lacked enough coverage ##
								my $num_carriers_missing = $num_carrier_called - $num_carrier_called_het - $num_carrier_called_lowfreq;
								$mendel_probability = $mendel_probability - ($num_carriers_missing * 0.20);
							}

							$mendel_probability = 0.20 if($mendel_probability < 0.20);							

						}
						elsif($mendel_status eq "Affected_Was_Wildtype")
						{
							## Here we'd have To have missed a variant. It does happen but rarely for homozygous variants. Assume 95% sensitivity
							## Another possibility is that there are two causal variants in a family pedigree. 
							$mendel_probability = 0.05;
						}
						elsif($mendel_status eq "Affected_Was_Het")
						{
							## Here we'd have To have missed a variant. It does happen. Assume 80% sensitivity
							## Another possibility is that there are two causal variants in a family pedigree. 
							$mendel_probability = 0.40;
						}
						elsif($mendel_status eq "Affected_Not_Called")
						{
							## HEre the affected was not called variant but we have low coverage. Assume 80% sensitivity ##
							$mendel_probability = 0.20;
						}
						elsif($mendel_status eq "Carrier_Was_Wildtype")
						{
							## A should-be carrier was wildtype for this variant, which means it was either missed or not inherited from both ##
							## Assume 80% sensitivity ##
							$mendel_probability = 0.20;
						}
						elsif($mendel_status eq "Carrier_Was_Homozygous")
						{
							## Carrier was homozygous, which is not what we expect for recessive. It could be a mis-call, however unlikely ##
							$mendel_probability = 0.60;
						}						
						elsif($mendel_status eq "Control_Was_Heterozygous")
						{
							## It is possible that we have an unaffected who's a carrier, although it's unlikely ##
							$mendel_probability = 0.80;
						}
						elsif($mendel_status eq "Control_Was_Homozygous")
						{
							## It is 
							$mendel_probability = 0.20;
						}						
					}


					
					
					## Get Shared IBD Probability ##
					
					my $shared_ibd_probability = 1;
					
					if($self->shared_ibd_file && $shared_ibd{$chrom})
					{
						my $ibd_score = get_shared_ibd($position, $shared_ibd{$chrom});
						
						if($ibd_score)
						{
							$shared_ibd_probability = $ibd_score;
						}
						else
						{
							$shared_ibd_probability = 0.5;
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
							my $newline = join("\t", $chrom, $position, $id, $ref, $var, $score, "Annotation_Error", $info, $format);
							for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
							{
								$newline .= "\t" . $lineContents[$colCounter];
							}
							print FILTERFILE "$newline\n" if($self->filtered_file);
							## Exclude variant with transcript-annotation errors ##
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

							if($variant_class)
							{
								$stats{'variants_type_' . $variant_type . '_with_annot'}++;
								$stats{'variants_type_' . $variant_type . '_with_annot_mendel_' . $mendel_status}++;
								
								## VARIANT CLASS PROBABILITY ##

								my $class_probability = 0.01;
								## Get damaging ##
								
								if($variant_class eq "nonsense" || $variant_class eq "nonstop" || $variant_class eq "frame_shift_ins" || $variant_class eq "frame_shift_del")
								{
									## Could be causal, though these are rare. ##
									$class_probability = 1.00;
								}
								elsif($variant_class eq "missense")
								{
									if(is_damaging($vep_polyphen, $vep_sift, $vep_condel))
									{
										## Damaging missense - very likely. No penalty ##
										$class_probability = 0.95;
									}
									else
									{
										## Not damaging missense. In HGMD, only 7.38% of variants look neutral to all three - slight penality of 20% ##
										$class_probability = 0.80;
									}
								}
								elsif($variant_class eq "splice_site")
								{
									## Not damaging, but less likely to cause disease ##
									$class_probability = 1.00; 
								}
								elsif($variant_class eq "splice_region")
								{
									## Splice region is promising, but less likely ##
									$class_probability = 0.20;
								}
								elsif($variant_class eq "silent" || $variant_class eq "in_frame_ins" || $variant_class eq "in_frame_del")
								{
									$class_probability = 0.05;
								}


								## Get Gene expression probability ##
#								my $gene_probability = 0.90;
								my $gene_probability = 0.50; # Set a default value at 50% rank for unknown genes #
								my $gene_expr_rank = "-";
								if($gene_expression{$variant_gene})
								{
									$stats{'variants_type_' . $variant_type . '_with_annot_had_expr'}++;
									my ($fpkm, $rank) = split(/\t/, $gene_expression{$variant_gene});
									## ADjust rank to value between zero and one ##
									$rank = $rank / 100;								
									$gene_expr_rank= sprintf("%.2f", $rank);

									## Set gene probability to 1 - rank ##
									
									$gene_probability = 1 - $rank;
									
#									if($rank < 0.25)
#									{
#										## If in the top 25%, let's up-rank this one ##
#										$gene_probability = 1.0;
#									}
#									elsif($rank < 0.50)
#									{
#										## If in the top 50%, let's up-rank this one ##
#										$gene_probability = 0.95;
#									}
#									elsif($rank > 0.75)
#									{
#										# Relatively very low expression, so let's penalize ##
#										$gene_probability = 0.80;
#									}
								}

								## IF we have a novel Mendel-passed variant, count it up ##
								if($mendel_status eq "PASS" && $dbsnp_status eq "novel")
								{
									$stats{'variants_type_' . $variant_type . '_with_annot_mendel_' . $mendel_status . '_dbsnp_' . $dbsnp_status}++;
									$stats{'variants_type_' . $variant_type . '_with_annot_mendel_' . $mendel_status . '_dbsnp_' . $dbsnp_status . '_class_' . $variant_class}++;
								}


								
								## We no longer remove mendelian variants. ##
								
								if($self->remove_mendel_failures && $mendel_status ne "PASS")
								{
									my $newline = join("\t", $chrom, $position, $id, $ref, $var, $score, $mendel_status, $info, $format);
									for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
									{
										$newline .= "\t" . $lineContents[$colCounter];
									}
									print FILTERFILE "$newline\n" if($self->filtered_file);
									$stats{'variants_type_' . $variant_type . '_with_annot_removed_for_mendel_failure'}++;								
								}
								else
								{
									## Determine the probability score ##
									my $final_probability = $mendel_probability * $pop_probability * $class_probability * $gene_probability * $shared_ibd_probability;
									## Create a new info field ##
									my $new_info = join(";", "RPscore=$final_probability", "Gene=$variant_gene");
									my $score_info = join(";", "MendelScore=$mendel_probability", "PopScore=$pop_probability", "EffectScore=$class_probability", "ExprScore=$gene_probability", "IbdScore=$shared_ibd_probability"); 
									$new_info .= ";" . "Annot=" . $our_annot if($our_annot);
									$new_info .= ";" . "VEP=" . $vep_annot if($vep_annot);
									$new_info .= ";" . $score_info if($score_info);
									$new_info .= ";" . $info if($info && $info ne ".");

									$candidate_variants[$num_candidates] = join("\t", $final_probability, $chrom, $position, $id, $ref, $var, $score, $filter, $new_info, $format);
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
											$candidate_variants[$num_candidates] .= "\t" . $lineContents[$colCounter];											
										}
									}
									$num_candidates++;
									
									## SAVE Candidate variant in REadable form ##
									my $priority_flag = 0;
									$priority_flag = $priority_genes{$variant_gene} if($priority_genes{$variant_gene});
									
									my $candidate_key = join("\t", $chrom, $position, $ref, $var);
									$candidate_variants_readable{$candidate_key} = join("\t", $priority_flag, $final_probability, $mendel_probability, $pop_probability, $class_probability, sprintf("%.2f", $gene_probability), $shared_ibd_probability);
									$candidate_variants_readable{$candidate_key} .= "\t" . join("\t", $mendel_status, $dbsnp_status, $id, $info);
									$candidate_variants_readable{$candidate_key} .= "\t" . join("\t", $variant_class, $variant_gene, $our_annot, $vep_annot);
									$candidate_variants_readable{$candidate_key} .= "\t$readable_genotypes";

								}


							}
							else
							{
								my $newline = join("\t", $chrom, $position, $id, $ref, $var, $score, "NoVariantClass", $info, $format);
								for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
								{
									$newline .= "\t" . $lineContents[$colCounter];
								}
								print FILTERFILE "$newline\n" if($self->filtered_file);
								
								print $transcript_annotation{$variant_key} . "\n";
								print $vep_annotation{$variant_key} . "\n";
								$stats{'variants_type_' . $variant_type . '_removed_annot_unknown'}++;
							}

						}

					}
					else
					{
						my $newline = join("\t", $chrom, $position, $id, $ref, $var, $score, "NoAnnotation", $info, $format);
						for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
						{
							$newline .= "\t" . $lineContents[$colCounter];
						}
						print FILTERFILE "$newline\n" if($self->filtered_file);
						$stats{'variants_type_' . $variant_type . '_removed_annot_missing'}++;

					}



				}
				else
				{
					my $newline = join("\t", $chrom, $position, $id, $ref, $var, $score, "NotVariant", $info, $format);
					for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
					{
						$newline .= "\t" . $lineContents[$colCounter];
					}
					print FILTERFILE "$newline\n" if($self->filtered_file);
					$stats{'non_variant_lines'}++;
				}
			}
			else
			{
				my $newline = join("\t", $chrom, $position, $id, $ref, $var, $score, $filter, $info, $format);
				for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
				{
					$newline .= "\t" . $lineContents[$colCounter];
				}
				print FILTERFILE "$newline\n" if($self->filtered_file);
				$stats{'lines_filter_' . $filter}++;							
			}
			
			
		} ## End If-else for header ##
		

	}
	
	close($input);



	foreach my $key (sort keys %stats)
	{
		print "$stats{$key}\t$key\n";
	}

	print "$num_candidates candidate variants\n";
	
	@candidate_variants = sort byProbability @candidate_variants;
	my %printed = ();
	foreach my $candidate (@candidate_variants)
	{
		my @lineContents = split(/\t/, $candidate);
		my $numContents = @lineContents;
		my $new_line = "";
		my $chrom = $lineContents[1];
		my $pos = $lineContents[2];
		my $ref = $lineContents[4];
		my $var = $lineContents[5];
		my $candidate_key = join("\t", $chrom, $pos, $ref, $var);
		
		if(!$printed{$candidate_key})
		{
			## Print every column except the first which was prob score for sorting #
			for(my $colCounter = 1; $colCounter < $numContents; $colCounter++)
			{
				$new_line .= "\t" if($new_line);
				$new_line .= $lineContents[$colCounter];
			}
			print OUTFILE "$new_line\n";
			
			print READABLEFILE join("\t", $candidate_key, $candidate_variants_readable{$candidate_key}) . "\n" if($self->readable_file);
			
			$printed{$candidate_key} = 1;
		}

	}
	
	for(my $canCounter = 0; $canCounter < 20; $canCounter++)
	{
		my @lineContents = split(/\t/, $candidate_variants[$canCounter]);
		my $score = $lineContents[0];
		my $chrom = $lineContents[1];
		my $pos = $lineContents[2];
		my $ref = $lineContents[3];
		my $var = $lineContents[4];
		my $new_info = $lineContents[8];
#		print join("\t", $chrom, $pos, $ref, $var, $new_info) . "\n";
	}
	
	sub byProbability
	{
		my ($prob_a) = split(/\t/, $a);
		my ($prob_b) = split(/\t/, $b);
		$prob_b <=> $prob_a;
	}
	

	if($self->output_file)
	{
		close(OUTFILE);
	}
	
	if($self->filtered_file)
	{
		close(FILTERFILE);
	}
	
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub get_shared_ibd
{
	my ($position, $lines) = @_;
	
	my @lines = split(/\n/, $lines);
	foreach my $line (@lines)
	{
		my ($chr_start, $chr_stop, $num_mark, $score) = split(/\t/, $line);
		if($chr_start <= $position && $chr_stop >= $position)
		{
			return($score);
		}
	}
	
	return(0);
}

################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

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


################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub load_expression
{
	my $annotation_file = shift(@_);
	my %annotation = ();
	
	my $num_genes = `cat $annotation_file | wc -l`; chomp($num_genes);
	
	my $input = new FileHandle ($annotation_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($gene, $fpkm) = split(/\t/, $line);
		
		my $expression_rank = $lineCounter / $num_genes * 100;
		
		$annotation{$gene} = join("\t", $fpkm, $expression_rank);

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

################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub load_priority_genes
{
	my $annotation_file = shift(@_);
	my %annotation = ();
	
	my $input = new FileHandle ($annotation_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($gene) = split(/\t/, $line);
		
		$annotation{$gene} = 1;
	}
	
	close($input);
	
	return(%annotation);
}




################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub load_shared_ibd
{
	my $annotation_file = shift(@_);
	my %annotation = ();
	
	my $input = new FileHandle ($annotation_file);
	my $lineCounter = 0;

	my $max_num_samples = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $chr_start, $chr_stop, $num_markers, $num_samples) = split(/\t/, $line);

		if($chrom && $chrom ne "chrom")
		{
			$max_num_samples = $num_samples if($num_samples > $max_num_samples);
		}

	}
	
	close($input);
	
	
	## Parse file again to calculate fraction of samples with shared iBD ##
	
	$input = new FileHandle ($annotation_file);
	$lineCounter = 0;


	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $chr_start, $chr_stop, $num_markers, $num_samples) = split(/\t/, $line);

		if($chrom && $chrom ne "chrom")
		{
			my $pct_max_samples = $num_samples / $max_num_samples;
			$annotation{$chrom} .= "\n" if($annotation{$chrom});
			$annotation{$chrom} .= join("\t", $chr_start, $chr_stop, $num_markers, $pct_max_samples);
		}

	}
	
	close($input);	
	
	return(%annotation);
}



################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub load_vep
{
    my $annotation_file = shift(@_);
    my %annotation = ();
    
    my %has_canonical = ();
    
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
	    my $canonical = 0;

            my @extraContents = split(/\;/, $extra);
            foreach my $entry (@extraContents)
            {
                    my ($key, $value) = split(/\=/, $entry);

                    $gene = $value if($key eq 'HGNC');
                    $polyphen = $value if($key eq 'PolyPhen');
                    $sift = $value if($key eq 'SIFT');
                    $condel = $value if($key eq 'Condel');
		    $canonical = $value if($key eq "CANONICAL");
            }

            my @classes = split(/\,/, $class);


		if($canonical && $canonical eq "YES")
		{
			## Reset any existing annotation ##
			$annotation{$key} = "";
			
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
			$has_canonical{$key} = 1;	
		}
		elsif($has_canonical{$key})
		{
			## Ignore this annotation, because the variant already has a canonical one ##
		}
		else
		{
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



#############################################################
# load_vep_results - parses the file
#
#############################################################

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


#############################################################
# load_vep_results - parses the file
#
#############################################################

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



#############################################################
# load_vep_results - parses the file
#
#############################################################

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

#############################################################
# convert_vep_class - convert to our trv_type
#
#############################################################
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

#############################################################
# get_code - get a VEP class rank for this annotation 
#
#############################################################

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





################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub convert_genotype
{
	my ($ref, $var, $genotype) = @_;

	return("NN") if($genotype eq '.' || $genotype eq "./." || $genotype eq "NN");
	
	## Determine type of variant we're dealing with ##
	
	my $variant_type = "snp";
	(my $test_var) = split(/\,/, $var);
	
	$variant_type = "del" if(length($ref) > 1);
	$variant_type = "ins" if(length($test_var) > 1);

#	print "$ref\t$var\t$genotype\t$variant_type\n";

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


sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;


