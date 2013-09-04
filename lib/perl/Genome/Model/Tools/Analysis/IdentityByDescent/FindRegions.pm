
package Genome::Model::Tools::Analysis::IdentityByDescent::FindRegions;     # rename this when you give the module file a different name <--

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
use Genome::Model::Tools::Analysis::Helpers qw(
    code_to_genotype_returning_code
);

my $bin_size = 1000000;

my %genotypes = ();

class Genome::Model::Tools::Analysis::IdentityByDescent::FindRegions {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variant_file	=> { is => 'Text', doc => "List of variants to consider (annotation format)", is_optional => 0, is_input => 1},
		affected_files	=> { is => 'Text', doc => "Consensus files for affected individuals", is_optional => 0, is_input => 1},
		unaffected_files	=> { is => 'Text', doc => "Consensus files for unaffected individuals", is_optional => 1, is_input => 1},
		inheritance_model	=> { is => 'Text', doc => "Presumed model of mendelian inheritance [autosomal-dominant]", is_optional => 1, is_input => 1},
		output_file	=> { is => 'Text', doc => "Output file for QC result", is_optional => 1, is_input => 1},
		print_all	=> { is => 'Text', doc => "If set to 1, prints all variants regardless of mendelian pattern", is_optional => 1, is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Identifies regions that are identical by descent among members of a family pedigree"                 
}

sub help_synopsis {
    return <<EOS
This command identifies regions that are identical by descent among members of a family pedigree.
EXAMPLE:	gmt analysis identity-by-descent find-regions --variant-file informative.snps --affected-files sample1.cns,sample2.cns --output-file variants.report.tsv
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
	my $inheritance_model = $self->inheritance_model;
	
	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";		
	}

	my %stats = ();
	
	## Build an array of affected individuals' genotypes ##
	
	my @affected_array = ();
	my $num_affected = my $num_unaffected = 0;
	
	my @affected_files = split(/\,/, $affected_files);
	my @unaffected_files = split(/\,/, $unaffected_files) if($unaffected_files);
	
	warn "Loading Affected samples...\n";
	
	## Count the files of each type and print the header ##
	my $header = join("\t", "chrom", "chr_start", "chr_stop", "ref", "var", "unaffecteds_variant", "affecteds_variant", "affecteds_missing", "affecteds_ambiguous");
	foreach my $affected_file (@affected_files)
	{
		$num_affected++;
		$header .= "\t" if($header);
		$header .= "AFF:" . $affected_file . "\treads1\treads2\tfreq";
		warn "$affected_file...\n";
		load_consensus($affected_file);
	}

	foreach my $unaffected_file (@unaffected_files)
	{
		$num_unaffected++;
		$header .= "\t" if($header);
		$header .= "UNAFF:" . $unaffected_file . "\treads1\treads2\tfreq";
		load_consensus($unaffected_file);
	}

	warn "$num_affected affected samples\n";
	warn "$num_unaffected unaffected samples\n";


	my %chrom_bin_starts = my %chrom_bin_stops = my %chrom_bin_sums = my %chrom_bin_nums = ();


	## PRint the header ##
	

	if($self->output_file)
	{
		print OUTFILE "variant\tnum_affected\tnum_unaffected\t$header\n";
	}


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

			## Assign this guy to a bin ##
			
			my ($bin_start, $bin_stop) = assign_to_bin($chr_start);
			my $bin_key = join("\t", $chrom, $bin_start, $bin_stop);
						
			my $sample_genotype_string = "";

			## AFFECTED GENOTYPES ##
			
			my $affecteds_variant = my $affecteds_wildtype = my $affecteds_missing = my $affecteds_ambiguous = my $unaffecteds_variant = 0;
			my $affecteds_homozygous = 0;
			
			## Compile Affected Genotypes ##
			
			my %affected_gts = ();
			
			foreach my $affected_file (@affected_files)
			{
#				my %genotypes = load_consensus($affected_file);
				my $sample_genotype = "-\t-\t-\t-";

				my $key = "$affected_file\t$chrom\t$chr_start";
				
				if($genotypes{$key})
				{
					(my $sample_call, my $sample_reads1, my $sample_reads2, my $sample_freq) = split(/\t/, $genotypes{$key});
					my $sample_coverage = $sample_reads1 + $sample_reads2;
					$sample_call = code_to_genotype_returning_code($sample_call);
					$affected_gts{$affected_file} = $sample_call;
				}
				else
				{
					$affecteds_missing++;
				}

			}

			## Only evaluate if we have all data ##
			
			if(!$affecteds_missing)
			{
				## Go through affected genotypes, scoring for this position ##
				## Possible states: 	0 - neither allele matches
				## 			1 - one allele matches, so possibly same haplotype
				##			2 - both alleles match, so possibly same haplotype
				
				my $score_sum = my $score_num = 0;
				## only take lower states ##
				
				my $gt_string = "";
				
				foreach my $affected_file (keys %affected_gts)
				{
					my $genotype = $affected_gts{$affected_file};
					$gt_string .= "\t" if($gt_string);
					$gt_string .= $genotype;
					
					foreach my $test_affected_file (keys %affected_gts)
					{
						if($affected_file ne $test_affected_file)
						{
							my $test_genotype = $affected_gts{$test_affected_file};
							## Do the comparison ##					
							my $score = compare_genotypes($genotype, $test_genotype);
							$score_sum += $score;
							$score_num++;
						}
					}
				}
	
				## Calculate the average score ##
				
				my $avg_score = 0;
				if($score_num)
				{
					$avg_score = sprintf("%.2f", $score_sum / $score_num);
					$chrom_bin_sums{$bin_key} += $avg_score;
					$chrom_bin_nums{$bin_key}++;
				}
	
	
	
				## Print the results ##
				
				if($bin_start <= 68669356 && $bin_stop >= 68669356)
				{
					print join("\t", $chrom, $chr_start, $chr_stop, $ref, $var, $gt_string, $score_num, $avg_score) . "\n";				
				}				
			}




		}		
		
	}
	
	close($input);
	
	if($self->output_file)
	{
		close(OUTFILE);

	}
	
	warn $stats{'num_variants'} . " variants\n";


	## Go through each bin, determining its average score ##
	
	foreach my $bin_key (keys %chrom_bin_nums)
	{
		if($chrom_bin_nums{$bin_key})
		{
			my $avg_bin_score = sprintf("%.4f", $chrom_bin_sums{$bin_key} / $chrom_bin_nums{$bin_key});
			print join("\t", $avg_bin_score, $chrom_bin_nums{$bin_key}, $bin_key) . "\n";			
		}

	}

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



################################################################################################
# Load Genotypes
#
################################################################################################

sub assign_to_bin
{
	my $position = shift(@_);
	
	my $bin_start = 1;
	my $bin_stop = $bin_size;
	
	while($position > $bin_stop)
	{
		$bin_start += $bin_size;
		$bin_stop += $bin_size;

		return($bin_start, $bin_stop) if($position <= $bin_stop);
	}
	
	return($bin_start, $bin_stop);
	
}

################################################################################################
# Load Genotypes
#
################################################################################################

sub compare_genotypes
{
	my ($gt1, $gt2) = @_;

	my $allele1a = substr($gt1, 0, 1);
	my $allele1b = substr($gt1, 1, 1);

	my $allele2a = substr($gt2, 0, 1);
	my $allele2b = substr($gt2, 1, 1);
	
	if($gt1 eq $gt2)
	{
		return(2);# if($allele1a ne $allele1b);
#		return(1);
	}


	return(1) if($allele1a eq $allele2a || $allele1a eq $allele2b || $allele1b eq $allele2a || $allele1b eq $allele2b);

	return(0);
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

1;
