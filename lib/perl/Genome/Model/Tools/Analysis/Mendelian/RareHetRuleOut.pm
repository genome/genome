
package Genome::Model::Tools::Analysis::Mendelian::RareHetRuleOut;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::Analysis::Mendelian::RareHetRuleOut {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		sample_file	=> { is => 'Text', doc => "Tab-delimited file of family, sample, status, dir", is_optional => 0, is_input => 1},
		family	=> { is => 'Text', doc => "Family for which to build mpileup", is_optional => 0, is_input => 1},
		min_cov_for_ref	=> { is => 'Text', doc => "The minimum coverage required to agree that a site is reference", is_optional => 0, is_input => 1, default => 10},
		reference	=> { is => 'Text', doc => "Path to the reference to use [defaults to build 37]", is_optional => 0, is_input => 1, default => '/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa'},
		chromosome	=> { is => 'Text', doc => "Name of chromosome to be analyzed", is_optional => 0, is_input => 1},
		common_vars	=> { is => 'Text', doc => "A list of common variant positions to exclude from analysis", is_optional => 0, is_input => 1},		
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
EXAMPLE:	gmt analysis mendelian rare-het-rule-out --sample-file Family-Sample-Status-Dir.tsv --family VCH004 
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
	my $output_basename = $self->output_basename;
	my $reference = $self->reference;
	
	my %stats = ();


	my %common = load_common($self->common_vars, $self->chromosome);

	my $sample_list = my $bam_list = "";
	my $case_list = "";

	my %snv_samples = my %snv_case_hets = my %snv_case_homs = ();
	my %snv_genotypes = ();
	my %control_snvs = ();

	## Print the variants ##

	my $input = new FileHandle ($sample_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($family, $sample_name, $affected_status, $dir) = split(/\t/, $line);	

		if($family eq $self->family)
		{
			my $bam_file = `ls $dir/alignments/*.bam | head -1`;
			chomp($bam_file);
			
			my $snvs_file = `ls $dir/variants/snvs.hq.bed`;
			chomp($snvs_file);


			warn "Loading SNVs for sample $sample_name...\n";
			my %sample_snvs = load_snvs($snvs_file, $self->chromosome);
			my $num_sample_snvs = 0;

			foreach my $snv (keys %sample_snvs)
			{
				my $genotype = $sample_snvs{$snv};
				$snv_samples{$snv}++;	## Count the number of samples with this variant ##
				
				my $key = join("\t", $sample_name, $snv);
				$snv_genotypes{$key} = $genotype;
				$num_sample_snvs++;
			}
			
			warn "$num_sample_snvs SNVs loaded\n";

			if($affected_status && $affected_status ne "control")
			{
				$stats{'family_affecteds'}++;
				$case_list .= "\n" if($case_list);
				$case_list .= $sample_name;				
			}
			else
			{
				$stats{'family_controls'}++;
				%control_snvs = %sample_snvs;
			}
			
			$sample_list .= "\t" if($sample_list);
			$sample_list .= $sample_name;
			
			$bam_list .= " " if($bam_list);
			$bam_list .= $bam_file;
		}

	}
	
	close($input);

	print $self->family . "\n";
	print $stats{'family_affecteds'} . " affecteds\n";
	print $stats{'family_controls'} . " controls\n";

	my @cases = split(/\n/, $case_list);

	## Now process the variants on the desired chromosome ##
	
	open(OUTFILE, ">$output_basename.variants") or die "Can't open outfile: $!\n";
	
	my %shared_rare_hets_per_bin = ();
	
	foreach my $snv (sort byPosition keys %snv_samples)
	{
		$stats{'num_snvs'}++;
		if($snv_samples{$snv} >= 2)
		{
			$stats{'num_snvs_min_samples'}++;
			
			if($common{$snv})
			{
				$stats{'num_snvs_min_samples_common'}++;
			}
			else
			{
				## Get the bin of this one ##
				
				my $bin = sprintf("%d", $snv / 100000);
				
				## Count the number of affecteds ##
				my $num_aff = my $num_aff_het = my $num_aff_hom = 0;
				foreach my $sample_name (@cases)
				{
					$num_aff++;
					
					my $key = join("\t", $sample_name, $snv);
					if($snv_genotypes{$key})
					{
						if(is_het($snv_genotypes{$key}))
						{
							$num_aff_het++;	
						}
						else
						{
						#	warn "$snv_genotypes{$key} not het\n";
						}
					}
				}
				
				if($num_aff_het == $num_aff)
				{
					$shared_rare_hets_per_bin{$bin}++;
				}
				
				
				$stats{'num_snvs_min_samples_rare'}++;
			}
		}
	
		print OUTFILE join("\t", $self->chromosome, $snv, $snv_samples{$snv}) . "\n";
	}

	close(OUTFILE);
	
	
	
	print $stats{'num_snvs'} . " unique SNV positions\n";
	print $stats{'num_snvs_min_samples'} . " seen in 2+ samples\n";
	print $stats{'num_snvs_min_samples_common'} . " common\n";
	print $stats{'num_snvs_min_samples_rare'} . " rare\n";
	sub byPosition
	{
		$a <=> $b;
	}
	
	open(OUTFILE, ">$output_basename.binSharedHets") or die "Can't open outfile: $!\n";
	print OUTFILE "chrom\tchr_start\tchr_stop\tshared_rare_hets\n";
	foreach my $bin (sort byPosition keys %shared_rare_hets_per_bin)
	{
		print OUTFILE join("\t", $self->chromosome, ($bin * 100000), (($bin * 100000) + 100000 - 1), $shared_rare_hets_per_bin{$bin}) . "\n";
	}

	close(OUTFILE);
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# is_het - returns true if variant call is heterozygous
#
################################################################################################

sub is_het
{
	my $genotype = shift(@_);
	my ($a1, $a2) = split(/\//, $genotype);
	
	if($a2 eq 'A' || $a2 eq 'C' || $a2 eq 'G' || $a2 eq 'T' || $a2 eq 'N')
	{
		return(0);
	}
	
	return(1);
}

################################################################################################
# is_hom - returns true if variant call is homozygous
#
################################################################################################

sub is_hom
{
	my $genotype = shift(@_);
	my ($a1, $a2) = split(/\//, $genotype);
	
	if($a2 eq 'A' || $a2 eq 'C' || $a2 eq 'G' || $a2 eq 'T')
	{
		if($a2 ne $a1)
		{
			return(1);
		}
	}
	
	return(0);
}





################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_snvs
{                               # replace with real execution logic.
	my ($FileName, $desired_chrom) = @_;
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	my %results = ();

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		my ($chrom, $chr_start, $chr_stop, $alleles) = split(/\t/, $line);

		if($chrom eq $desired_chrom)
		{
			## Save the variant ##
			$results{$chr_stop} = $alleles;
		}
	}
	
	close($input);

	return(%results);
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_common
{                               # replace with real execution logic.
	my ($FileName, $desired_chrom) = @_;
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	my %results = ();

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		my ($chrom, $position) = split(/\t/, $line);

		if($chrom eq $desired_chrom)
		{
			## Save the variant ##
			$results{$position} = 1;
		}
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


