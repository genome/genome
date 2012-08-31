
package Genome::Model::Tools::Analysis::Mendelian::IdentityByDescent;     # rename this when you give the module file a different name <--

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

my %genotypes = ();

class Genome::Model::Tools::Analysis::Mendelian::IdentityByDescent {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variant_file	=> { is => 'Text', doc => "List of variants to consider (annotation format)", is_optional => 0, is_input => 1},
		affected_files	=> { is => 'Text', doc => "Consensus files for affected individuals", is_optional => 0, is_input => 1},
		unaffected_files	=> { is => 'Text', doc => "Consensus files for unaffected individuals", is_optional => 1, is_input => 1},
		inheritance_model	=> { is => 'Text', doc => "Presumed model of mendelian inheritance [autosomal-dominant]", is_optional => 1, is_input => 1},
		min_coverage_to_refute	=> { is => 'Text', doc => "Minimum coverage to refute a possible variant in an affected [10]", is_optional => 1, is_input => 1, default => 10},
		max_frequency_to_refute	=> { is => 'Text', doc => "Maximum observed variant allele frequency to refute a possible variant in an affected [5]", is_optional => 1, is_input => 1, default => 10},
		min_affected_variant	=> { is => 'Text', doc => "Minimum number of affecteds with variant to include", is_optional => 1, is_input => 1, default => 1},
		max_unaffected_variant	=> { is => 'Text', doc => "Maximum number of affecteds with variant to include", is_optional => 1, is_input => 1, default => 0},
		ucsc_cytoband	=> { is => 'Text', doc => "Path to the UCSC refGene.txt file", is_optional => 0, is_input => 1, default => '/gscuser/dkoboldt/SNPseek/SNPseek2/ucsc/cytoBand.txt'},
		output_file	=> { is => 'Text', doc => "Output file for QC result", is_optional => 1, is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Identifies regions that are identical by descent in affected samples"                 
}

sub help_synopsis {
    return <<EOS
This command identifies regions that are identical-by-descent in affected samples
EXAMPLE:	gmt analysis mendelian report-variants --variant-file mySNPs.tier1 --affected-files sample1.cns,sample2.cns --output-file variants.report.tsv
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

	my $variant_file = $self->variant_file;
	my $affected_files = $self->affected_files;
	my $unaffected_files = $self->unaffected_files if($self->unaffected_files);	

	## Get user-defined settings ##
	my $min_affecteds_variant = $self->min_affected_variant;
	my $max_unaffecteds_variant = $self->max_unaffected_variant;
	my $min_coverage_to_refute = $self->min_coverage_to_refute;
	my $max_frequency_to_refute = $self->max_frequency_to_refute;

	my $inheritance_model = $self->inheritance_model;

	## Load UCSC Cytogenetic band information ##

	my $ucsc_cytoband = $self->ucsc_cytoband;
	my %bands_by_chrom = load_bands($ucsc_cytoband);
	
	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";		
	}
	
	my %stats = ();
	$stats{'num_variants'} = $stats{'num_variants_complete'} = $stats{'num_variants_incomplete'} = 0;
	
	## Build an array of affected individuals' genotypes ##
	
	my @affected_array = ();
	my $num_affected = my $num_unaffected = 0;
	
	my @affected_files = split(/\,/, $affected_files);
	my @unaffected_files = split(/\,/, $unaffected_files) if($unaffected_files);
	
	print "Loading Affected samples...\n";
	
	## Count the files of each type and print the header ##
	my $header = join("\t", "chrom", "chr_start", "chr_stop", "ref", "var");
	foreach my $affected_file (@affected_files)
	{
		$num_affected++;
		$header .= "\t" if($header);
		$header .= $affected_file;
		load_consensus($affected_file);
	}

	print "$num_affected affected samples\n";

	my %band_score_sum = my %band_score_num = ();

	## Print the variants ##

	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $chrom, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);

		if(length($ref) > 1 || $var eq "-")
		{
			## Undo the adjustment made when formatting deletions for annotation.
			$chr_start--;
		}

		my $variant_type = "SNP";
		
		if(length($ref) > 1 || $var eq "-")
		{
			$variant_type = "DEL";
		}
		elsif(length($var) > 1 || $ref eq "-")
		{
			$variant_type = "INS";
		}

		if($lineCounter >= 0)
		{
			$stats{'num_variants'}++;

			## AFFECTED GENOTYPES ##
			
			my $affecteds_with_genotype = my $affecteds_variant = my $affecteds_wildtype = my $affecteds_missing = my $affecteds_ambiguous = my $unaffecteds_variant = 0;
			my $affecteds_homozygous = 0;
			
			## See how many affecteds carry it ##
			my %sample_gt_counts = ();
			my %genotypes_by_sample = ();
			
			
			foreach my $affected_file (@affected_files)
			{
#				my %genotypes = load_consensus($affected_file);
				my $sample_genotype = "-\t-\t-\t-";

				my $key = "$affected_file\t$chrom\t$chr_start";
				
				if($genotypes{$key})
				{
					$affecteds_with_genotype++;
					(my $sample_call, my $sample_reads1, my $sample_reads2, my $sample_freq) = split(/\t/, $genotypes{$key});
					$sample_call = code_to_genotype($sample_call);
					$sample_gt_counts{$sample_call}++;
					$genotypes_by_sample{$affected_file} = $sample_call;
				
					$sample_genotype = "$sample_call\t$sample_reads1\t$sample_reads2\t$sample_freq";
				}
				else
				{
					$affecteds_missing++;
				}

			}

			
			if($affecteds_with_genotype == $num_affected)
			{
				$stats{'num_variants_complete'}++;
			
				## Get cytogenetic band ##
				
				my $cyto_band = match_to_band($chrom, $chr_start, $bands_by_chrom{$chrom});

				my $identity_score = 0;
				
				## Calculate the score ##
				
				my %pair_counted = ();
				
				foreach my $sample1 (sort keys %genotypes_by_sample)
				{
					foreach my $sample2 (sort keys %genotypes_by_sample)
					{
						if($sample1 ne $sample2 && !$pair_counted{$sample1 . '_' . $sample2})
						{
							my $genotype1 = $genotypes_by_sample{$sample1};
							my $genotype2 = $genotypes_by_sample{$sample2};
							

							
							if($genotype1 eq $genotype2)
							{
								$identity_score++;
								if($genotype1 eq $ref . $ref)
								{
									## Get one point for both being wild-type ##
#									$identity_score += 1;
								}
								elsif(substr($genotype1, 0, 1) eq substr($genotype1, 1, 1))
								{
									## Get two points for both being homozygous ##
#									$identity_score += 2;
								}
								else
								{
									## Get three points for both being the same het ##
#									$identity_score += 3;
								}
							}
							else
							{
								$identity_score--;
								if($genotype1 ne $ref . $ref && $genotype2 ne $ref . $ref)
								{
#									$identity_score++;
									## Get one point for neither matching the reference, although they differ ##
									
#									$identity_score += 1;
								}
								else
								{
#									$identity_score--;
									## Subtract a point for being different ##
#									$identity_score -= 1;
								}
							}
							
							## Count all possible combinations ##
							$pair_counted{$sample1 . '_' . $sample2} = 1;
							$pair_counted{$sample2 . '_' . $sample1} = 1;
						}
					}
				}
				

				my $string = "";
				foreach my $genotype (keys %sample_gt_counts)
				{
					$string .= "$sample_gt_counts{$genotype} $genotype,";
				}
				
				$band_score_sum{$cyto_band} += $identity_score;
				$band_score_num{$cyto_band}++;
				
#				print join("\t", $chrom, $chr_start, $chr_stop, $ref, $var, $cyto_band, $string, $identity_score) . "\n";
			}
			else
			{
				$stats{'num_variants_incomplete'}++;
			}

		}		
		
	}
	
	close($input);
	
	if($self->output_file)
	{
		close(OUTFILE);
	}
	
	print $stats{'num_variants'} . " variants\n";
	print $stats{'num_variants_incomplete'} . " with incomplete genotype data\n";
	print $stats{'num_variants_complete'} . " with complete genotype data\n";

	my @band_scores = ();
	my $num_bands = 0;

	foreach my $band (sort keys %band_score_sum)
	{
		if($band_score_num{$band})
		{
			if($band_score_num{$band} > 1 && !($band =~ 'NT' || $band =~ 'MT'))
			{
				my $avg_score = $band_score_sum{$band} / $band_score_num{$band};
				$avg_score = sprintf("%.3f", $avg_score);
				$band_scores[$num_bands] = "$avg_score\t$band_score_num{$band}\t$band";
				$num_bands++;				
			}
		}
	}

	@band_scores = sort byScore @band_scores;
	
	sub byScore
	{
		my ($score_a) = split(/\t/, $a);
		my ($score_b) = split(/\t/, $b);
		$score_b <=> $score_a;
	}


	## Print the top 10 bands ##
	
	for (my $counter = 0; $counter < 10; $counter++)
	{
		print "$band_scores[$counter]\n";
	}

#	print "And $num_bands bands\n";
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




################################################################################################
# Load Genotypes
#
################################################################################################

sub load_consensus
{                               # replace with real execution logic.
	my $genotype_file = shift(@_);
#	my %genotypes = ();
	
	my $input = new FileHandle ($genotype_file);
	my $lineCounter = 0;
	my $gtCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $chrom, my $position, my $ref, my $cns, my $reads1, my $reads2, my $var_freq) = split(/\t/, $line);

		if($ref =~ /[0-9]/)
		{
			($chrom, $position, my $stop, $ref, $cns, $reads1, $reads2, $var_freq) = split(/\t/, $line);			
		}

		if(length($ref) > 1 || length($cns) > 1 || $ref eq "-" || $cns eq "-")
		{
			## If CNS is not formatted, do so ##

			if(!($cns =~ '/'))
			{
				$cns = "$ref/$cns";				
			}

		}
#		my $key = "$chrom\t$position";
		my $key = "$genotype_file\t$chrom\t$position";
		if($cns ne "N")
		{
			$genotypes{$key} = "$cns\t$reads1\t$reads2\t$var_freq";
		}

	}
	close($input);
                            # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)

#	return(%genotypes);
}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_bands
{
	my $FileName = shift(@_);
	my %genes_by_chrom = ();

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		if($lineCounter >= 1)
		{
			my @lineContents = split(/\t/, $line);
			my $chrom = $lineContents[0];
			my $tx_start = $lineContents[1];
			my $tx_stop = $lineContents[2];
			my $gene_symbol = $chrom . $lineContents[3];
			
			$chrom =~ s/chr//g;
			$tx_start++;
			$tx_stop++;
			
			if($genes_by_chrom{$chrom})
			{
				$genes_by_chrom{$chrom} .= "\n";
			}
			
			$genes_by_chrom{$chrom} .= join("\t", $tx_start, $tx_stop, $gene_symbol);
#			print join("\t", $chrom, $tx_start, $tx_stop, $gene_symbol) . "\n";
		}


	}
	
	close($input);	
	
	return(%genes_by_chrom);
}



################################################################################################
# Load Genotypes
#
################################################################################################

sub match_to_band
{
	my ($chrom, $position, $bands) = @_;

	if($bands)
	{
		my @bands = split(/\n/, $bands);
		
		foreach my $band (@bands)
		{
			my ($band_start, $band_stop, $band_name) = split(/\t/, $band);
			
			($band_name) = split(/\./, $band_name);
		
			if($position >= $band_start && $position <= $band_stop)
			{
				$band_name = substr($band_name, 0, length($band_name) - 2);
				return($band_name);
			}
		}		
	}

	
	return("chr" . $chrom);
}


sub code_to_genotype
{
	my $code = shift(@_);
	
	return("AA") if($code eq "A");
	return("CC") if($code eq "C");
	return("GG") if($code eq "G");
	return("TT") if($code eq "T");

	return("AC") if($code eq "M");
	return("AG") if($code eq "R");
	return("AT") if($code eq "W");
	return("CG") if($code eq "S");
	return("CT") if($code eq "Y");
	return("GT") if($code eq "K");

#	warn "Unrecognized ambiguity code $code!\n";
#	return("NN");
	return($code);
}


sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;

