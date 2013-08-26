
package Genome::Model::Tools::Analysis::Mendelian::RareHetRuleOutVcf;     # rename this when you give the module file a different name <--

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


use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

my $num_affected_called = my $affecteds_missing = my $unaffecteds_variant = my $affecteds_variant = my $affecteds_ambiguous = 0;
my %genotypes = ();

class Genome::Model::Tools::Analysis::Mendelian::RareHetRuleOutVcf {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		vcf_file	=> { is => 'Text', doc => "Input in VCF format", is_optional => 0, is_input => 1},
		output_basename	=> { is => 'Text', doc => "Output basename for files", is_optional => 0, is_input => 1},
		control_samples	=> { is => 'Text', doc => "Comma-separated list of control sample names", is_optional => 1, is_input => 1},
		male_samples	=> { is => 'Text', doc => "Comma-separated list of male sample names", is_optional => 1, is_input => 1},
		ignore_samples	=> { is => 'Text', doc => "Comma-separated list of sample names to ignore", is_optional => 1, is_input => 1},
		inheritance_model	=> { is => 'Text', doc => "Presumed model of mendelian inheritance [autosomal-dominant]", is_optional => 1, is_input => 1, default => 'autosomal-dominant'},
		min_coverage	=> { is => 'Text', doc => "Minimum coverage to refute a possible variant in an affected", is_optional => 1, is_input => 1, default => 20},
		min_call_rate	=> { is => 'Text', doc => "Minimum callrate for affecteds to include a variant", is_optional => 1, is_input => 1, default => 0.50},
		max_frequency_to_refute	=> { is => 'Text', doc => "Maximum observed variant allele frequency to refute a possible variant in an affected [5]", is_optional => 1, is_input => 1, default => 10},
		min_affected_variant	=> { is => 'Text', doc => "Minimum number of affecteds with variant to include", is_optional => 1, is_input => 1, default => 1},
		max_unaffected_variant	=> { is => 'Text', doc => "Maximum number of unaffecteds with variant to include", is_optional => 1, is_input => 1, default => 0},
		plot_results	=> { is => 'Text', doc => "If set to 1, generate per-chromosome plots of results", is_optional => 1, is_input => 1, default => 0},
		centromere_file	=> { is => 'Text', doc => "A UCSC 0-based BED file of centromere locations per chromosome", is_optional => 0, is_input => 1, default => '/gscmnt/sata809/info/medseq/dkoboldt/SNPseek2/ucsc/hg19/gapTable.centromere.nochr.txt'},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Attempts to use haplotype and inheritance rules to include/exclude regions of shared haplotypes among affecteds"                 
}

sub help_synopsis {
    return <<EOS
This tool attempts to use haplotype and inheritance rules to include/exclude regions of shared haplotypes among affecteds
EXAMPLE:	gmt analysis mendelian rare-het-rule-out-vcf --vcf-file myVCF.vcf --output-basename myVCF.rhro
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 
This tool attempts to use haplotype and inheritance rules to include/exclude regions of shared haplotypes among affecteds
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

	my %ignore_sample = ();
	if($self->ignore_samples)
	{
		my @samples = split(/\,/, $self->ignore_samples);
		foreach my $sample (@samples)
		{
			$ignore_sample{$sample} = 1;
			$stats{'num_samples_ignored'}++;
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



	warn "Parsing VCF...\n";
	
	
	

	my @candidate_variants = ();
	my $num_candidates = 0;
	my %candidate_variants_readable = ();
	
	my $output_file = $self->output_basename . ".counts.tsv";
	my $windows_file = $self->output_basename . ".windows.tsv";

	if($self->plot_results && $self->plot_results == 2)
	{
		do_plot_results($output_file, $windows_file, $self);
		do_single_plot($output_file, $windows_file, $self);
		exit(0);
	}

	my %centromeres = load_centromeres($self);

	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	print OUTFILE "chrom\tposition\tref\tvar\tdbsnp_status\tmendel_score\tis_hom_diff\tpct_affected_het\tcases\tcases_ref\tcases_het\tcases_hom\tcases_na\tcontrols\tcontrols_ref\tcontrols_het\tcontrols_hom\tcontrols_na\n";
	
	open(WINDOWFILE, ">$windows_file") or die "Can't open outfile: $!\n";
	print WINDOWFILE "chrom\tchr_start\tchr_stop\tnum_variants\n";
	
	my %window = ();
	

	
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
				for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
				{
					my $sample = $lineContents[$colCounter];
					$sample_names[$colCounter] = $lineContents[$colCounter];
					$num_samples++;
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


					
					## PART 1: GET THE MENDEL STATUS ##
					
					my $mendel_status = "PASS";
					my %marker_counts = ();
					
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
						my @genotypeContents = split(/\:/, $lineContents[$colCounter]);
						my $genotype = $genotypeContents[$genotype_column{'GT'}];
						my $coverage = 0;
						$coverage = $genotypeContents[$genotype_column{'DP'}];
						$coverage = 0 if($coverage && $coverage eq '.');
						my $filter = $genotypeContents[$genotype_column{'FT'}];

						## Get the highest var freq ##
						my $var_depth = $genotypeContents[$genotype_column{'AD'}];
						my $var_freq = 0;
						
						if($var_depth)
						{
							my @var_depths = split(/\,/, $var_depth);
							@var_depths = sort numericallyDesc @var_depths;
							$var_depth = $var_depths[0];
							$var_depth = 0 if($var_depth && $var_depth eq ".");
							$var_freq = $var_depth / $coverage if($coverage);				
						}



						my $freq_alt = $genotypeContents[$genotype_column{'FA'}];
						my $sample_status = "";

						if($control_sample{$sample_name})
						{
							$sample_status = "control";
						}
						else
						{
							$sample_status = "case";
						}
						
						$marker_counts{$sample_status}++;
						
						## Only process the genotype if it has a value and is either unfiltered or for the control sample ##
						if(length($genotype) > 2 && $genotype ne '.' && ($filter eq 'PASS' || $filter eq '.' || $control_sample{$sample_name}) && !$ignore_sample{$sample_name})
						{
#							warn "Trying to convert $genotype in column $colCounter at line $lineCounter\n";
							$genotype = convert_genotype($ref, $var, $genotype) if($genotype ne '.');


							my $allele1 = substr($genotype, 0, 1);
							my $allele2 = substr($genotype, 1, 1);
							$allele_counts{$allele1}++;
							$allele_counts{$allele2}++;
	
							my $gt = "Missing";
							
							if($coverage >= $self->min_coverage)
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
										$stats{'wildtype_with_allele_depth'}++;
									}

								}
								elsif(is_homozygous($genotype))
								{
									if($var_freq > 0.90)
									{
										$gt = "Hom";										
									}
									else
									{
										$gt = "Het";
									}

								}
								elsif(is_heterozygous($genotype))
								{
									$gt = "Het";
								}								
							}



							$marker_counts{$sample_status . "_" . $gt}++;
							
							if($gt ne "Missing")
							{
								$marker_counts{$sample_status . "_called"}++;
							}

						}
						else
						{
							## Genotype missing or was filtered out ####
							$genotype = "NN";
							$readable_genotypes .= "\t" if($readable_genotypes);
							$coverage = "-" if(!$coverage);
							$freq_alt = "-" if(!$freq_alt);
							$readable_genotypes .= join("\t", $genotype, $coverage, $freq_alt);
						}

						
					}


					$marker_counts{'case'} = 0 if(!$marker_counts{'case'});
					$marker_counts{'control'} = 0 if(!$marker_counts{'control'});
					$marker_counts{'case_called'} = 0 if(!$marker_counts{'case_called'});
					$marker_counts{'control_called'} = 0 if(!$marker_counts{'control_called'});					

					$marker_counts{'case_Ref'} = 0 if(!$marker_counts{'case_Ref'});
					$marker_counts{'control_Ref'} = 0 if(!$marker_counts{'control_Ref'});
					$marker_counts{'case_Het'} = 0 if(!$marker_counts{'case_Het'});
					$marker_counts{'control_Het'} = 0 if(!$marker_counts{'control_Het'});					
					$marker_counts{'case_Hom'} = 0 if(!$marker_counts{'case_Hom'});
					$marker_counts{'control_Hom'} = 0 if(!$marker_counts{'control_Hom'});					
					$marker_counts{'case_Missing'} = 0 if(!$marker_counts{'case_Missing'});
					$marker_counts{'control_Missing'} = 0 if(!$marker_counts{'control_Missing'});
					
					## Get centromere start/stop ##
					
					
					my $case_call_rate = 0;
					$case_call_rate = $marker_counts{'case_called'} / $marker_counts{'case'} if($marker_counts{'case'});
					
					my $pct_affected_het = "NA";
					$pct_affected_het = sprintf("%.3f", $marker_counts{'case_Het'} / $marker_counts{'case_called'}) if($marker_counts{'case_called'});
					
					## Determine basic mendel score ##
					my $basic_mendel_score = 0;
					if($pct_affected_het ne "NA" && $pct_affected_het == 1)
					{
						$basic_mendel_score = 0.5;
						
						if($marker_counts{'control_Het'} == 0)
						{
							$basic_mendel_score = 1;
						}
					}
					
					## Adjust for chromosome X ##
					if($chrom eq "X")
					{
						if(($marker_counts{'case_Het'} + $marker_counts{'case_Hom'}) == $marker_counts{'case_called'})
						{
							$basic_mendel_score = 0.5;
							if($marker_counts{'control_Het'} == 0 && $marker_counts{'control_Hom'} == 0)
							{
								$basic_mendel_score = 1;
							}
						}
					}
					
					## Determine Hom-diff ##
					
					my $is_hom_diff = 0;
					
					if($marker_counts{'case_Ref'} > 0 && $marker_counts{'case_Hom'} > 0)
					{
#						my $num_diff = $marker_counts{'case_Ref'};
#						$num_diff = $marker_counts{'case_Hom'} if($marker_counts{'case_Hom'} < $num_diff);					
#						$is_hom_diff = $num_diff / $marker_counts{'case_called'};



						$is_hom_diff = 1;
						## Try calculating this as the # of samples that are hom Diff ##
#						$is_hom_diff = $marker_counts{'case_Ref'};
#						$is_hom_diff = $marker_counts{'case_Hom'} if($marker_counts{'case_Hom'} < $is_hom_diff);

					}
					
					## Determine probability based on mendel segregation status, where any error must be a wrong variant call ##
			

					
					if($case_call_rate >= $self->min_call_rate && ($dbsnp_status eq "novel" || $is_hom_diff > 0)) #$dbsnp_status ne "common" && $dbsnp_status ne "uncommon")
					{
						$stats{'variants_included'}++;
						$stats{'variants_included_and_homdiff'}++ if($is_hom_diff);

						my $output_line = join("\t", $chrom, $position, $ref, $var, $dbsnp_status, $basic_mendel_score, $is_hom_diff, $pct_affected_het);
						$output_line .= "\t" . join("\t", $marker_counts{'case_called'}, $marker_counts{'case_Ref'}, $marker_counts{'case_Het'}, $marker_counts{'case_Hom'}, $marker_counts{'case_Missing'});
						$output_line .= "\t" . join("\t", $marker_counts{'control_called'}, $marker_counts{'control_Ref'}, $marker_counts{'control_Het'}, $marker_counts{'control_Hom'}, $marker_counts{'control_Missing'});
						print OUTFILE "$output_line\n";
						
						
						## Handle window stuff ##
						
						if($window{'chrom'} && $window{'chrom'} ne $chrom)
						{
							## Ended a chromosome, so report the window then reset ##
#							print WINDOWFILE join("\t", $window{'chrom'}, $window{'start'}, $window{'stop'}, $window{'variants'}) . "\n" if($window{'variants'} > 1 || ($window{'stop'} - $window{'start'} + 1) >= 1000000);
							print_window($self, $window{'chrom'}, $window{'start'}, $window{'stop'}, $window{'variants'}) if($window{'variants'} > 1 || ($window{'stop'} - $window{'start'} + 1) >= 1000000);
							%window = ();
						}
						

						
						if($is_hom_diff)
						{
							## Hit Homdiff, so report window and reset ##
							if($window{'chrom'})
							{
								$window{'stop'} = $position - 1;
#								print WINDOWFILE join("\t", $window{'chrom'}, $window{'start'}, $window{'stop'}, $window{'variants'}) . "\n" if($window{'variants'} > 1 || ($window{'stop'} - $window{'start'} + 1) >= 1000000);
								print_window($self, $window{'chrom'}, $window{'start'}, $window{'stop'}, $window{'variants'}) if($window{'variants'} > 1 || ($window{'stop'} - $window{'start'} + 1) >= 1000000);
								%window = ();
							}
							else
							{
								## Otherwise move start up to the last-homdiff position ##
								$window{'last_homdiff_chrom'} = $chrom;
								$window{'last_homdiff_pos'} = $position;
							}
						}
						elsif($basic_mendel_score == 1 || ($pct_affected_het && $pct_affected_het == 1))
						{
							## We have a Mendelian consistency or complete hetness among affecteds, so start or continue a window ##
							if($window{'chrom'})
							{
								$window{'stop'} = $position;
								$window{'variants'}++;
							}
							else
							{
								$window{'chrom'} = $chrom;
								## If we had a prior homdiff position on this chromosome, use that as the window start ##
								if($window{'last_homdiff_chrom'} && $window{'last_homdiff_chrom'} eq $chrom && $window{'last_homdiff_pos'} && $window{'last_homdiff_pos'} < $position)
								{
									$window{'start'} = $window{'last_homdiff_pos'} + 1;
								}
								## Otherwise start at this position ##
								else
								{
									$window{'start'} = $position;
								}

								$window{'stop'} = $position;
								$window{'variants'} = 1;
							}
						}
						
					}


				}
				## Otherwise not a variant

			}
			## Otherwise failed filter status
			
			
		} ## End If-else for header ##
		

	}
	
	close($input);


	if($window{'chrom'})
	{
		## Ended a chromosome, so report the window then reset ##
#		print WINDOWFILE join("\t", $window{'chrom'}, $window{'start'}, $window{'stop'}, $window{'variants'}) . "\n" if($window{'variants'} > 1 || ($window{'stop'} - $window{'start'} + 1) >= 1000000);
		print_window($self, $window{'chrom'}, $window{'start'}, $window{'stop'}, $window{'variants'}) if($window{'variants'} > 1 || ($window{'stop'} - $window{'start'} + 1) >= 1000000);
		%window = ();
	}


	close(OUTFILE);
	close(WINDOWFILE);

	if($self->plot_results)
	{
		do_plot_results($output_file, $windows_file, $self);
	}


	foreach my $key (sort keys %stats)
	{
		print "$stats{$key}\t$key\n";
	}
	
	
	sub print_window
	{
		my ($self, $chrom, $start, $stop, $variants) = @_;
		my %centromeres = load_centromeres($self);
		my $cen_start = my $cen_stop = 0;		
		($cen_start, $cen_stop) = split(/\t/, $centromeres{$chrom}) if($centromeres{$chrom});

		if($cen_start && $cen_stop && $start <= $cen_start && $stop >= $cen_stop)
		{
			## Window contains entire centromere, so adjust ##
			print WINDOWFILE join("\t", $chrom, $start, $cen_start, $variants) . "\n";
			print WINDOWFILE join("\t", $chrom, $cen_stop, $stop, $variants) . "\n";		
		}
		elsif($cen_start && $cen_stop && $start <= $cen_start && $stop >= $cen_start)
		{
			## WIndow runs into centromere ##
			print WINDOWFILE join("\t", $chrom, $start, $cen_start, $variants) . "\n";
		}
		elsif($cen_start && $cen_stop && $start <= $cen_stop && $stop >= $cen_stop)
		{
			## Window starts within centromere
			print WINDOWFILE join("\t", $chrom, $cen_stop, $stop, $variants) . "\n";
		}
		elsif($cen_start && $cen_stop && $start >= $cen_start && $stop <= $cen_stop)
		{
			## Window within centromere - don't print ##
		}
		else
		{

			print WINDOWFILE join("\t", $chrom, $start, $stop, $variants) . "\n"
		}

	}


	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub do_plot_results
{
	my $infile = shift(@_);
	my $windows_file = shift(@_);
	my $self = shift(@_);

	my %centromeres = load_centromeres($self);

	open(INDEX, ">" . $self->output_basename . ".plot.index.html") or die "Can't open outfile: $!\n";
	print INDEX "<HTML><BODY><TABLE CELLSPACING=0 CELLPADDING=5 BORDER=0 WIDTH=\"100%\">\n";
	print INDEX "<TR>\n";
	my $num_printed_in_column = 0;

	my $script_filename = $self->output_basename . ".R";
	open(SCRIPTFILE, ">" . $script_filename) or die "Can't open output file for R script: $!\n";

	my $plot_basename = $self->output_basename . ".plot";

	print SCRIPTFILE qq {count <- read.table("$infile", header=TRUE)\n};
	print SCRIPTFILE qq {windows <- read.table("$windows_file", header=TRUE)\n};	
	print SCRIPTFILE qq {library(ggplot2)\n};
	
	my %rhro_by_chrom = load_rhro_regions($windows_file);
	
	for(my $chrCounter = 1; $chrCounter <= 23; $chrCounter++)
	{
		my $chrom = $chrCounter;
		$chrom = "X" if($chrCounter == 23);
		
		my $cen_start = my $cen_stop = -1;
		($cen_start, $cen_stop) = split(/\t/, $centromeres{$chrom}) if($centromeres{$chrom});
		
#		print SCRIPTFILE qq{
#png("$plot_basename.$chrom.png", height=600, width=800)
#par(mar=c(4, 4, 1, 1) + 0.1)
#plot(count\$position[count\$chrom=="$chrom" & count\$is_hom_diff == 1], count\$is_hom_diff[count\$chrom=="$chrom" & count\$is_hom_diff == 1] - 0.01, pch=19, cex=0.5, type="h", col="red", xlim=c(0,max(count\$position[count\$chrom=="$chrom"]) + 1), xlab="Position on chromosome $chrom", ylab="Fraction of Affecteds", ylim=c(0,1), main="chr$chrom")
#points(count\$position[count\$chrom=="$chrom"], count\$pct_affected_het[count\$chrom=="$chrom"], pch=19, cex=0.5, col="lightskyblue")
#points(count\$position[count\$chrom=="$chrom" & count\$pct_affected_het == 1], count\$pct_affected_het[count\$chrom=="$chrom" & count\$pct_affected_het == 1], pch=19, cex=0.5, col="blue", type="h")
#points(count\$position[count\$chrom=="$chrom" & count\$mendel_score == 1], count\$mendel_score[count\$chrom=="$chrom" & count\$mendel_score == 1] + 0.01, pch=19, cex=0.5, col="green")
#};


#png("$plot_basename.$chrom.png", height=600, width=800)

################# THIS WAS THE ORIGINAL PLOTTING ######################
#		print SCRIPTFILE qq{
#png("$plot_basename.$chrom.png", height=200, width=800)
#par(mar=c(4, 4, 1, 1) + 0.1)
#plot(count\$position[count\$chrom=="$chrom"], count\$pct_affected_het[count\$chrom=="$chrom"], pch=19, cex=0.5, xlim=c(0,max(count\$position[count\$chrom=="$chrom"]) + 1), col="lightskyblue", xlab="Position on chromosome $chrom", ylab="Fraction of Affecteds", ylim=c(0,1), main="chr$chrom")
#points(count\$position[count\$chrom=="$chrom" & count\$pct_affected_het == 1], count\$pct_affected_het[count\$chrom=="$chrom" & count\$pct_affected_het == 1], pch=19, cex=0.5, col="blue")
#points(count\$position[count\$chrom=="$chrom" & count\$mendel_score == 1], count\$mendel_score[count\$chrom=="$chrom" & count\$mendel_score == 1] + 0.01, pch=19, cex=0.5, col="green")
#points(count\$position[count\$chrom=="$chrom" & count\$is_hom_diff == 1], count\$is_hom_diff[count\$chrom=="$chrom" & count\$is_hom_diff == 1] - 0.01, pch=19, cex=0.5, col="red")
#};
#		print SCRIPTFILE qq|\nif(length(windows\$chr_start[windows\$chrom=="$chrom"]) > 0) {\n
#segments(windows\$chr_start[windows\$chrom=="$chrom"], 1.02, windows\$chr_stop[windows\$chrom=="$chrom"], 1.02, col="black", lwd=2)\n
#}\n|;
#
#		print SCRIPTFILE qq|\nif(length(windows\$chr_start[windows\$chrom=="$chrom"]) > 0) {\n
#segments(windows\$chr_start[windows\$chrom=="$chrom"], 0.98, windows\$chr_stop[windows\$chrom=="$chrom"], 0.98, col="black", lwd=2)\n
#}\n|;
#
#		print SCRIPTFILE qq{abline(v=$cen_start, lty=2)
#abline(v=$cen_stop, lty=2);
#legend("bottomright", c("Some Het", "All Het", "Fits Mendel", "Rule Out"), pch=c(19, 19, 19), col=c("lightskyblue", "blue","green","red"))
#dev.off()
#};	

################### END ORIGINAL PLOTTING ##################################

################# GGPLOT INSTEAD ################
		print SCRIPTFILE qq{
png("$plot_basename.$chrom.png", height=200, width=800)
par(mar=c(4, 4, 1, 1) + 0.1)
p <- ggplot(count, aes(x = position[count\$chrom=="$chrom"], y = pct_affected_het[count\$chrom=="$chrom"])) + geom_point(color="blue") + xlab("Position on chr$chrom") + ylab("% Affecteds Het")
q <- geom_point(data=count,aes(position[count\$chrom=="$chrom" & count\$is_hom_diff >= 1],is_hom_diff[count\$chrom=="$chrom" & count\$is_hom_diff >= 1]),color="red")
m <- geom_point(data=count,aes(position[count\$chrom=="$chrom" & count\$mendel_score == 1],mendel_score[count\$chrom=="$chrom" & count\$mendel_score == 1]),color="green")
};

my $plot_line = "p + opts(title=\"chr$chrom\") + scale_x_continuous(limits = c(0,max(count\$position[count\$chrom==\"$chrom\"]))) + q + m ";
#my $plot_line = "p + opts(title=\"chr$chrom\") + scale_x_continuous(limits = c(30000000,130000000)) + q + m ";

if($rhro_by_chrom{$chrom})
{
	my $region_no = 0;
	my @rhro_lines = split(/\n/, $rhro_by_chrom{$chrom});
	
	
	foreach my $rhro_line (@rhro_lines)
	{
		$region_no++;
		my ($region_start, $region_stop) = split(/\t/, $rhro_line);

		print SCRIPTFILE qq{
rect$region_no <- data.frame (xmin=$region_start, xmax=$region_stop, ymin=-Inf, ymax=Inf)
r$region_no <- geom_rect(data=rect$region_no, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="lightskyblue", fill="lightskyblue", alpha=0.5, inherit.aes = FALSE)
};
		$plot_line .= "+ r$region_no ";
	}
}

if($cen_start)
{
	$plot_line .= "+ geom_vline(xintercept = $cen_start, linetype = \"longdash\")";
}

if($cen_stop)
{
	$plot_line .= "+ geom_vline(xintercept = $cen_stop, linetype = \"longdash\")";
}


print SCRIPTFILE "$plot_line\n";
print SCRIPTFILE "dev.off()\n";
#		print SCRIPTFILE qq{abline(v=$cen_start, lty=2)
#abline(v=$cen_stop, lty=2);
#legend("bottomright", c("Rare Het", "Rule Out"), pch=c(19, 19), col=c("lightskyblue", "red"))
#dev.off()
#};	


################## END GGPLOT ################

		my @temp = split(/\//, $self->output_basename);
		my $numContents = @temp;
		my $true_basename = $temp[$numContents - 1];
		my $image_filename = "$true_basename.plot.$chrom.png";

#		print INDEX "<TD><A HREF=\"$image_filename\"><IMG SRC=\"$image_filename\" HEIGHT=240 WIDTH=320 BORDER=0></A></TD>\n";
		print INDEX "<TD><A HREF=\"$image_filename\"><IMG SRC=\"$image_filename\" HEIGHT=100 WIDTH=400 BORDER=0></A></TD>\n";

		$num_printed_in_column++;

		if($num_printed_in_column >= 4)
		{
			print INDEX "</TR><TR>\n";
			$num_printed_in_column = 0;
		}

	}

	close(SCRIPTFILE);

	system("R --no-save < $script_filename");

	print INDEX "</TR></TABLE></BODY></HTML>\n";
	close(INDEX);
}



################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub do_single_plot
{
	my $infile = shift(@_);
	my $windows_file = shift(@_);
	my $self = shift(@_);

	my %centromeres = load_centromeres($self);

	open(INDEX, ">" . $self->output_basename . ".plot.single.html") or die "Can't open outfile: $!\n";
	print INDEX "<HTML><BODY><TABLE CELLSPACING=0 CELLPADDING=0 BORDER=0>\n";
	print INDEX "<TR>\n";
	my $num_printed_in_column = 0;

	my $script_filename = $self->output_basename . ".single.R";
	open(SCRIPTFILE, ">" . $script_filename) or die "Can't open output file for R script: $!\n";

	my $plot_basename = $self->output_basename . ".single.plot";

	print SCRIPTFILE qq {count <- read.table("$infile", header=TRUE)\n};
	print SCRIPTFILE qq {windows <- read.table("$windows_file", header=TRUE)\n};	
	print SCRIPTFILE qq {library(ggplot2)\n};
	
	my %rhro_by_chrom = load_rhro_regions($windows_file);
	
	for(my $chrCounter = 1; $chrCounter <= 23; $chrCounter++)
	{
		my $chrom = $chrCounter;
		$chrom = "X" if($chrCounter == 23);
		
		my $cen_start = my $cen_stop = -1;
		($cen_start, $cen_stop) = split(/\t/, $centromeres{$chrom}) if($centromeres{$chrom});
		

################# GGPLOT INSTEAD ################
		print SCRIPTFILE qq{
png("$plot_basename.$chrom.png", height=150, width=150)
par(mar=c(4, 4, 1, 1) + 0.1)
p <- ggplot(count, aes(x = position[count\$chrom=="$chrom"], y = pct_affected_het[count\$chrom=="$chrom"])) + geom_point(color="blue") + xlab("") + ylab("% Affecteds Het")
q <- geom_point(data=count,aes(position[count\$chrom=="$chrom" & count\$is_hom_diff >= 1],is_hom_diff[count\$chrom=="$chrom" & count\$is_hom_diff >= 1]),color="red")
m <- geom_point(data=count,aes(position[count\$chrom=="$chrom" & count\$mendel_score == 1],mendel_score[count\$chrom=="$chrom" & count\$mendel_score == 1]),color="green")
};

my $plot_line = "p ";
#if($chrom eq "1")
#{
#	$plot_line .= "+ opts(title=\"chr$chrom\", plot.margin = unit(c(0.01,0.01,0.01,0.01), \"cm\"), axis.ticks=theme_blank(), axis.text.x=theme_blank(),axis.title.x=theme_blank()) + q + m ";
#	
#}
#else
#{
	$plot_line .= "+ opts(title=\"chr$chrom\", plot.margin = unit(c(0.01,0.01,0.01,0.01), \"cm\"), axis.ticks=theme_blank(), axis.text.x=theme_blank(),axis.title.x=theme_blank(), axis.text.y=theme_blank(),axis.title.y=theme_blank()) + q + m ";
#}


if($rhro_by_chrom{$chrom})
{
	my $region_no = 0;
	my @rhro_lines = split(/\n/, $rhro_by_chrom{$chrom});
	foreach my $rhro_line (@rhro_lines)
	{
		$region_no++;
		my ($region_start, $region_stop) = split(/\t/, $rhro_line);

		print SCRIPTFILE qq{
rect$region_no <- data.frame (xmin=$region_start, xmax=$region_stop, ymin=-Inf, ymax=Inf)
r$region_no <- geom_rect(data=rect$region_no, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="lightskyblue", fill="lightskyblue", alpha=0.5, inherit.aes = FALSE)
};
		$plot_line .= "+ r$region_no ";
	}
}

if($cen_start)
{
	$plot_line .= "+ geom_vline(xintercept = $cen_start, linetype = \"longdash\")";
}

if($cen_stop)
{
	$plot_line .= "+ geom_vline(xintercept = $cen_stop, linetype = \"longdash\")";
}


print SCRIPTFILE "$plot_line\n";
print SCRIPTFILE "dev.off()\n";


		my @temp = split(/\//, $self->output_basename);
		my $numContents = @temp;
		my $true_basename = $temp[$numContents - 1];
		my $image_filename = "$true_basename.single.plot.$chrom.png";

#		print INDEX "<TD><A HREF=\"$image_filename\"><IMG SRC=\"$image_filename\" HEIGHT=240 WIDTH=320 BORDER=0></A></TD>\n";
		print INDEX "<TD><A HREF=\"$image_filename\"><IMG SRC=\"$image_filename\" HEIGHT=150 WIDTH=150 BORDER=0></A></TD>\n";

		$num_printed_in_column++;

		if($num_printed_in_column == 11)
		{
			print INDEX "<TD> &nbsp; </TD></TR><TR>\n";
#			$num_printed_in_column = 0;
		}

	}

	close(SCRIPTFILE);

	system("R --no-save < $script_filename");

	print INDEX "</TR></TABLE></BODY></HTML>\n";
	close(INDEX);
}





################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub load_rhro_regions
{
	my $FileName = shift;
	my %regions = ();
	## Parse the VCF file ##
	
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $chr_start, $chr_stop, $variants) = split(/\t/, $line);
#		if($variants >= 2)
#		{
			$regions{$chrom} .= "\n" if($regions{$chrom});
			$regions{$chrom} .= join("\t", $chr_start, $chr_stop, $variants);			
#		}

	}
	
	close($input);

	return(%regions);
}




################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub load_centromeres
{
	my $self = shift(@_);
	my %centromeres = ();
	## Parse the VCF file ##
	
	my $input = new FileHandle ($self->centromere_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $chr_start, $chr_stop) = split(/\t/, $line);
		$centromeres{$chrom} = join("\t", $chr_start, $chr_stop);
	}
	
	close($input);

	return(%centromeres);
}


sub numericallyDesc
{
	$b = 0 if($b eq '.');
	$a = 0 if($a eq '.');
	$b <=> $a;
}

################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub is_homozygous
{
	my $gt = shift(@_);
	my $a1 = substr($gt, 0, 1);
	my $a2 = substr($gt, 1, 1);
	return(1) if($a1 eq $a2);
	return(0);
}

################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub is_heterozygous
{
	my $gt = shift(@_);
	my $a1 = substr($gt, 0, 1);
	my $a2 = substr($gt, 1, 1);
	return(1) if($a1 ne $a2);
	return(0);
}

################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

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


sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;


