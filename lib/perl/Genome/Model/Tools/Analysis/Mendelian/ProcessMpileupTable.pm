
package Genome::Model::Tools::Analysis::Mendelian::ProcessMpileupTable;     # rename this when you give the module file a different name <--

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

my $num_affected = my $affecteds_missing = my $unaffecteds_variant = my $affecteds_variant = my $affecteds_ambiguous = 0;
my %genotypes = ();

class Genome::Model::Tools::Analysis::Mendelian::ProcessMpileupTable {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		sample_file	=> { is => 'Text', doc => "Tab-delimited file of family, sample, status, dir", is_optional => 0, is_input => 1},
		variant_file	=> { is => 'Text', doc => "Variants in sorted annotation format to build mpileup for", is_optional => 0, is_input => 1},
		family	=> { is => 'Text', doc => "Family for which to build mpileup", is_optional => 0, is_input => 1},
		min_cov_for_ref	=> { is => 'Text', doc => "The minimum coverage required to agree that a site is reference", is_optional => 0, is_input => 1, default => 10},
		reference	=> { is => 'Text', doc => "Path to the reference to use [defaults to build 37]", is_optional => 0, is_input => 1, default => '/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa'},
		cns_file	=> { is => 'Text', doc => "VarScan mpileup2cns results for mpileup file", is_optional => 1, is_input => 1},
		output_basename	=> { is => 'Text', doc => "Output file basename for merged results", is_optional => 1, is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Runs a Mendelian analysis of variants on a per-family basis"                 
}

sub help_synopsis {
    return <<EOS
This command processes the mpileup CNS in conjunction with the variant list
EXAMPLE:	gmt analysis mendelian process-mpileup-table --sample-file Family-Sample-Status-Dir.tsv --variant-file VCH004.variants --family VCH004 --cns-file VCH004.mpileup.cns --output-file VCH004.variants.mendelian
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

	my $sample_file = $self->sample_file;
	my $cns_file = $self->cns_file;
	my $output_basename = $self->output_basename;
	my $variant_file = $self->variant_file;
	my $reference = $self->reference;
	
	my %stats = ();

	my $sample_list = my $bam_list = "";

	## Load the mpileup results ##
	
	my %mpileup_results = load_mpileup_cns($cns_file);

	## Open the output file ##
	
	open(OUTPASS, ">$output_basename.pass") or die "Can't open outfile: $!\n";
	open(OUTFAIL, ">$output_basename.fail") or die "Can't open outfile: $!\n";

	## Print the variants ##

	my $input = new FileHandle ($sample_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($family, $sample_name, $affected_status, $dir) = split(/\t/, $line);	

		if($family eq $self->family && $affected_status && $affected_status ne "control")
		{
			$stats{'family_affecteds'}++;
			my $bam_file = `ls $dir/alignments/*.bam | head -1`;
			chomp($bam_file);

			$sample_list .= "\t" if($sample_list);
			$sample_list .= $sample_name;
			
			$bam_list .= " " if($bam_list);
			$bam_list .= $bam_file;
		}
		elsif($family eq $self->family)
		{
			$stats{'family_controls'}++;
		}
	}
	
	close($input);
	
	print "$sample_list\n$bam_list\n";


	my @samples = split(/\t/, $sample_list);
	my $num_samples = @samples;

	$input = new FileHandle ($variant_file);
	$lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);		

		## Determine variant type ##
		my $variant_type = "SNP";
		
		if($ref eq '-' || $ref eq '0' || length($var) > 1)
		{
			$variant_type = "INS";
		}
		elsif($var eq '-' || $var eq '0' || length($ref) > 1)
		{
			$variant_type = "DEL";
		}	

		my @lineContents = split(/\t/, $line);
		my $num_columns = @lineContents;

		if($lineCounter == 1)
		{
			## Print the start of header ##
			print OUTPASS "chrom\tchr_start\tchr_stop\tref\tvar\t";
			print OUTFAIL "chrom\tchr_start\tchr_stop\tref\tvar\t";
			## Print header to output files ##
			for(my $colCounter = 5; $colCounter < $num_columns; $colCounter++)
			{
				print OUTPASS "\t";
				print OUTFAIL "\t";
			}
			## Print the sample names ##
			print OUTPASS "$sample_list\n";
			print OUTFAIL "$sample_list\n";
			
		}
		
		$stats{'num_variants'}++;
		
		my $fail_flag = 0;
		my $fail_reason = "";

		print join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);

		if($variant_type eq "INS" || $variant_type eq "DEL")
		{
			my $indel_size = length($var);
			$indel_size = length($ref) if($variant_type eq "DEL");

			## Find the appropriate mpileup line ##
			my $best_mpileup_line = "";
			my $best_distance = 0;
			my $best_size_diff = 0;
			
			for(my $position = ($chr_start - $indel_size - 1); $position <= ($chr_stop + $indel_size + 1); $position++)
			{
				my $key = join("\t", $chrom, $position);
				if($mpileup_results{$key})
				{
					my ($mp_chrom, $mp_pos, $mp_ref, $mp_var) = split(/\t/, $mpileup_results{$key});

					if($mp_var ne '.' && length($mp_var) > 1)
					{
						my $mp_indel_size = length($mp_var) - 1;
						my $mp_indel_distance = abs($chr_start - $position);
						my $mp_size_diff = abs($indel_size - $mp_indel_size);
						
						## See if the indel type matches expected ##
						## If multiple possible indels, closest to expected size wins it. If that's a tie, closest to annotated position wins it ##
						if($variant_type eq "INS" && substr($mp_var, 0, 1) eq '+')
						{							
							if(!$best_mpileup_line || $mp_size_diff < $best_size_diff || ($mp_size_diff <= $best_size_diff && $mp_indel_distance < $best_distance))
							{
								$best_mpileup_line = $mpileup_results{$key};
								$best_distance = $mp_indel_distance;
							}
						}
						elsif($variant_type eq "DEL" && substr($mp_var, 0, 1) eq '-')
						{
							if(!$best_mpileup_line || $mp_size_diff < $best_size_diff || ($mp_size_diff <= $best_size_diff && $mp_indel_distance < $best_distance))
							{
								$best_mpileup_line = $mpileup_results{$key};
								$best_distance = $mp_indel_distance;
							}							
						}
					}
				}
			}
			
			
			## IF there's no best mpileup line, give it the chromosome position ##
			if(!$best_mpileup_line)
			{
				$best_mpileup_line = $mpileup_results{"$chrom\t$chr_start"};
			}
			
			## Proceeed if we have the best mpileup line ##
			if($best_mpileup_line)
			{
				my $sample_genotypes = "";
				$stats{'num_with_mpileup'}++;
	
	
				my @mpileup = split(/\t/, $best_mpileup_line);
				my @cns = split(/\s+/, $mpileup[10]);
				for(my $sampleCounter = 0; $sampleCounter < $num_samples; $sampleCounter++)
				{
					my $sample_name = $samples[$sampleCounter];
					my $sample_cns = $cns[$sampleCounter];
	
					my ($call, $cov, $reads1, $reads2, $freq, $pvalue) = split(/\:/, $sample_cns);
	
					$sample_genotypes .= "\t" if($sample_genotypes);
					$sample_genotypes .= $sample_cns;
					
					if(length($call) == 1 && $call ne "N" && $cov >= $self->min_cov_for_ref && $reads1 >= $self->min_cov_for_ref)
					{
						$fail_flag++;
						$fail_reason .= "; " if($fail_reason);
						$fail_reason .= "Affected $sample_name was wild-type";
					}
				}
				
				if($fail_flag)
				{
					$stats{'num_fail_mendelian'}++;
					print "\tFAIL\t$fail_reason\n";
					print OUTFAIL join("\t", $line, $sample_genotypes, $fail_reason) . "\n";
				}
				else
				{
					print "\tPASS\t$sample_genotypes\n";
					$stats{'num_pass_mendelian'}++;
					print OUTPASS join("\t", $line, $sample_genotypes, $fail_reason) . "\n";
				}				
			}
			else
			{
				my $sample_genotypes = "";
				$fail_flag = 1;
				$fail_reason = "No matching indels in mpileup2cns";
				foreach my $sample (@samples)
				{
					$sample_genotypes .= "\t" if($sample_genotypes);
					$sample_genotypes .= "--";
				}
				print OUTFAIL join("\t", $line, $sample_genotypes, $fail_reason) . "\n";			
	
				print "\tFAIL\t$fail_reason\n";				
			}
		}
		else
		{
			my $key = join("\t", $chrom, $chr_start);
			if($mpileup_results{$key})
			{
				my $sample_genotypes = "";
				$stats{'num_with_mpileup'}++;
	
#				print join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
	
				my @mpileup = split(/\t/, $mpileup_results{$key});
				my @cns = split(/\s+/, $mpileup[10]);
				for(my $sampleCounter = 0; $sampleCounter < $num_samples; $sampleCounter++)
				{
					my $sample_name = $samples[$sampleCounter];
					my $sample_cns = $cns[$sampleCounter];
	
					my ($call, $cov, $reads1, $reads2, $freq, $pvalue) = split(/\:/, $sample_cns);
	
					$sample_genotypes .= "\t" if($sample_genotypes);
					$sample_genotypes .= $sample_cns;
					
					if($call eq $ref && $cov >= $self->min_cov_for_ref && $reads1 >= $self->min_cov_for_ref)
					{
						$fail_flag++;
						$fail_reason .= "; " if($fail_reason);
						$fail_reason .= "Affected $sample_name was wild-type";
					}
				}
				
				if($fail_flag)
				{
					$stats{'num_fail_mendelian'}++;
					print "\tFAIL\t$fail_reason\n";
					print OUTFAIL join("\t", $line, $sample_genotypes, $fail_reason) . "\n";
				}
				else
				{
					print "\tPASS\t$sample_genotypes\n";
					$stats{'num_pass_mendelian'}++;
					print OUTPASS join("\t", $line, $sample_genotypes, $fail_reason) . "\n";
				}
			}
			else
			{
				my $sample_genotypes = "";
				$fail_flag = 1;
				$fail_reason = "No mpileup2cns results\n";
				foreach my $sample (@samples)
				{
					$sample_genotypes .= "\t" if($sample_genotypes);
					$sample_genotypes .= "--";
				}
				print OUTFAIL join("\t", $line, $sample_genotypes, $fail_reason) . "\n";			
	
				print "\tFAIL\t$fail_reason\n";
			}			
		}

		
#		my $cmd = "samtools mpileup -f $reference -q 10 -r $query $bam_list 2>/dev/null >>$output_file";
#		print join("\t", $lineCounter, $line) . "\n";
#		system($cmd);
		
#		return(0) if($lineCounter > 10);
	}
	
	close($input);

	print $self->family . "\n";
	print $stats{'family_affecteds'} . " affecteds\n";
	print $stats{'family_controls'} . " controls\n";
	print $stats{'num_variants'} . " variants\n";
	print $stats{'num_with_mpileup'} . " had mpileup2cns results\n";
	print $stats{'num_fail_mendelian'} . " failed Mendelian inheritance\n";
	print $stats{'num_pass_mendelian'} . " passed Mendelian inheritance\n";
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_mpileup_cns
{                               # replace with real execution logic.
	my $FileName = shift(@_);
	my %results = ();
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		my ($chrom, $position) = split(/\t/, $line);
		my $key = join("\t", $chrom, $position);
		$results{$key} = $line;
	}
	
	close($input);

	return(%results);
}

################################################################################################
# sorting subroutine by chromosome and then position
#
################################################################################################

sub byChrPos
{
	my ($chrom_a, $pos_a) = split(/\t/, $a);
	my ($chrom_b, $pos_b) = split(/\t/, $b);

	$chrom_a = 23 if($chrom_a eq 'X');
	$chrom_a = 24 if($chrom_a eq 'Y');
	$chrom_a = 25 if($chrom_a eq 'MT');
	
	$chrom_b = 23 if($chrom_b eq 'X');
	$chrom_b = 24 if($chrom_b eq 'Y');
	$chrom_b = 25 if($chrom_b eq 'MT');

	$chrom_a =~ s/[^0-9]//g;
	$chrom_b =~ s/[^0-9]//g;
	
	$chrom_a <=> $chrom_b
	or
	$pos_a <=> $pos_b;
}


1;


