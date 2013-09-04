
package Genome::Model::Tools::Analysis::Mendelian::GermlineComparison;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::Analysis::Mendelian::GermlineComparison {
    is => 'Command',                       

    has => [                                # specify the command's single-value properties (parameters) <--- 
        sample_file	=> { is => 'Text', doc => "Tab-delimited file of family, sample, status, dir", is_optional => 0, is_input => 1},
        inheritance_model	=> { is => 'Text', doc => "Mendelian inheritance model to use", is_optional => 0, is_input => 1, default => 'autosomal-dominant'},
        reference	=> { is => 'Text', doc => "Path to the reference to use [defaults to build 37]", is_optional => 0, is_input => 1, example_values=> ['/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa']},
        output_dir	=> { is => 'Text', doc => "Output directory to contain files", is_optional => 1, is_input => 1},
        tier1_only	=> { is => 'Text', doc => "If set to 1, restrict variants to tier 1 annotation types", is_optional => 1, is_input => 1},
        outside_tier1_only	=> { is => 'Text', doc => "If set to 1, only report sites outside tier 1", is_optional => 1, is_input => 1},
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Runs a Mendelian analysis of variants on a per-family basis"                 
}

sub help_synopsis {
    return <<EOS
This command runs a Mendelian analysis of variants on a per-family basis
EXAMPLE:	gmt analysis mendelian germline-comparison --sample-file Family-Sample-Status-Dir.tsv --output-dir mendelian_out
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
	my $output_dir = $self->output_dir;
	
	my %stats = ();

	my %family_affecteds = my %family_controls = my %all_controls = ();

	## Print the variants ##

	my $input = new FileHandle ($sample_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($family, $sample_name, $affected_status, $dir) = split(/\t/, $line);	

		if($affected_status && $affected_status ne "control" && lc($affected_status) ne "unaffected")
		{
			$family_affecteds{$family} .= "\n" if($family_affecteds{$family});
			$family_affecteds{$family} .= $sample_name . "\t" . $dir;
		}
		else
		{
			$family_controls{$family} .= "\n" if($family_controls{$family});
			$family_controls{$family} .= $sample_name . "\t" . $dir;
			
			$all_controls{$sample_name} = $dir;
		}
	}
	
	close($input);
	
	## Build a master list of variants in all controls ##
	my %control_variants = ();
	print "Loading control variants...\n";
	foreach my $sample (keys %all_controls)
	{
		print "$sample\n";
		my $dir = $all_controls{$sample};
		
		my %variants = load_variants($dir);

		foreach my $key (keys %variants)
		{
			$control_variants{$key}++;
		}
	}
	
	
	print "Analyzing Families...\n";
	
	## Go through each family and print a summary, then process ##
	
	foreach my $family (sort keys %family_affecteds)
	{
		open(SUMMARY, ">$output_dir/$family.comparison.out") or die "Can't open outfile: $!\n";
#		if($family eq "VCH026")
#		{
		## Reset family stats ##
		my %family_stats = ();
		$family_stats{'num_affecteds'} = $family_stats{'num_controls'} = 0;


		## Count number of affecteds and controls ##

		my @affecteds = split(/\n/, $family_affecteds{$family});
		$family_stats{'num_affecteds'} = @affecteds;

		my @family_controls = ();
		
		if($family_controls{$family})
		{
			@family_controls = split(/\n/, $family_controls{$family});
			$family_stats{'num_controls'} = @family_controls;
		}

		print join("\t", $family, "$family_stats{'num_affecteds'} affecteds, $family_stats{'num_controls'} controls") . "\n";	
		print SUMMARY join("\t", $family, "$family_stats{'num_affecteds'} affecteds, $family_stats{'num_controls'} controls") . "\n";	
		
		
		print SUMMARY "sample_name\tsnps\tindels\n";

		## Prepare to compile SNVs and indels ##

		my %family_snvs = my %family_indels = my %family_annotations = ();

		## Go through the cases, loading the variants for each one ##
		print "Loading affecteds...\n";
		foreach my $affected (@affecteds)
		{
			my ($sample_name, $dir) = split(/\t/, $affected);
			my %variants = load_variants($dir);

			print "$sample_name... ";
			my $num_snvs = my $num_indels = 0;
			
			foreach my $variant (keys %variants)
			{
				my ($chrom, $chr_start, $chr_stop, $ref, $var, $var_type) = split(/\t/, $variants{$variant});
				my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);

				my @varContents = split(/\t/, $variants{$variant});
				my $numContents = @varContents;
				my $annotation = "";
				for(my $colCounter = 6; $colCounter < $numContents; $colCounter++)
				{
					$annotation .= "\t" if($annotation);
					$annotation .= $varContents[$colCounter];
				}
				
				$family_annotations{$key} = $annotation;

				if($var_type eq "SNP")
				{
					$num_snvs++;
					$family_snvs{$key} .= "\n" if($family_snvs{$key});
					$family_snvs{$key} .= join("\t", $sample_name);
				}
				else
				{
					$num_indels++;
					$family_indels{$key} .= "\n" if($family_indels{$key});
					$family_indels{$key} .= join("\t", $sample_name);					
				}
			}

			print SUMMARY join("\t", $sample_name, $num_snvs, $num_indels) . "\n";
			print "$num_snvs SNVs, $num_indels indels\n";
		}
		
		
		## Open output files for SNVs ##
		
		open(NONSILENT, ">$output_dir/$family.snvs.nocontrol.nonsilent") or die "Can't open outfile: $!\n";
		open(NOCONTROL, ">$output_dir/$family.snvs.nocontrol") or die "Can't open outfile: $!\n";
		open(CONTROL, ">$output_dir/$family.snvs.control") or die "Can't open outfile: $!\n";
		print CONTROL "chrom\tchr_start\tchr_stop\tref\tvar\tnum_affected\tnum_control\n";
		print NOCONTROL "chrom\tchr_start\tchr_stop\tref\tvar\tnum_affected\tnum_control\n";
		print NONSILENT "chrom\tchr_start\tchr_stop\tref\tvar\tnum_affected\tnum_control\n";
		
		foreach my $snv (sort byChrPos keys %family_snvs)
		{
			my $annotation = $family_annotations{$snv};

			my @affecteds_variant = split(/\n/, $family_snvs{$snv});
			my $num_affecteds_variant = @affecteds_variant;
			my $affecteds_variant_list = $family_snvs{$snv};
			$affecteds_variant_list =~ s/\n/\,/g;
			my $num_controls_variant = 0;
			$num_controls_variant = $control_variants{$snv} if($control_variants{$snv});
			
			$family_stats{'snvs'}++;
			if($control_variants{$snv})
			{
				$family_stats{'snvs_control'}++;
				print CONTROL join("\t", $snv, $num_affecteds_variant, $num_controls_variant, $annotation) . "\n";
			}
			else
			{
				$family_stats{'snvs_nocontrol'}++;
				## Determine which affecteds have it ##

				print NOCONTROL join("\t", $snv, $num_affecteds_variant, $num_controls_variant, $annotation) . "\n";

				if(!($annotation =~ "\trna\t" || $annotation =~ "\tsilent\t"))
				{
					$family_stats{'snvs_nocontrol_nonsilent'}++;
					print NONSILENT join("\t", $snv, $num_affecteds_variant, $num_controls_variant, $annotation) . "\n";
					$family_stats{'snvs_nocontrol_nonsilent_' . $num_affecteds_variant . 'affecteds'}++;
				}
			}
		}

		close(NONSILENT);
		close(NOCONTROL);
		close(CONTROL);


		open(NONSILENT, ">$output_dir/$family.indels.nocontrol.nonsilent") or die "Can't open outfile: $!\n";
		open(NOCONTROL, ">$output_dir/$family.indels.nocontrol") or die "Can't open outfile: $!\n";
		open(CONTROL, ">$output_dir/$family.indels.control") or die "Can't open outfile: $!\n";
		print CONTROL "chrom\tchr_start\tchr_stop\tref\tvar\tnum_affected\tnum_control\n";
		print NOCONTROL "chrom\tchr_start\tchr_stop\tref\tvar\tnum_affected\tnum_control\n";
		print NONSILENT "chrom\tchr_start\tchr_stop\tref\tvar\tnum_affected\tnum_control\n";


		foreach my $indel (sort byChrPos keys %family_indels)
		{
			my $annotation = $family_annotations{$indel};
			my @affecteds_variant = split(/\n/, $family_indels{$indel});
			my $num_affecteds_variant = @affecteds_variant;
			my $affecteds_variant_list = $family_indels{$indel};
			$affecteds_variant_list =~ s/\n/\,/g;
			my $num_controls_variant = 0;
			$num_controls_variant = $control_variants{$indel} if($control_variants{$indel});

			$family_stats{'indels'}++;
			if($control_variants{$indel})
			{
				$family_stats{'indels_control'}++;
				print CONTROL join("\t", $indel, $num_affecteds_variant, $num_controls_variant, $annotation) . "\n";				
			}
			else
			{
				$family_stats{'indels_nocontrol'}++;
				## Determine which affecteds have it ##

				print NOCONTROL join("\t", $indel, $num_affecteds_variant, $num_controls_variant, $annotation) . "\n";				

				if(!($annotation =~ "\trna\t" || $annotation =~ "\tsilent\t"))
				{
					$family_stats{'indels_nocontrol_nonsilent'}++;
					print NONSILENT join("\t", $indel, $num_affecteds_variant, $num_controls_variant, $annotation) . "\n";				
					$family_stats{'indels_nocontrol_nonsilent_' . $num_affecteds_variant . 'affecteds'}++;
				}
			}
		}		

		close(NONSILENT);
		close(NOCONTROL);
		close(CONTROL);

		
		print SUMMARY $family_stats{'snvs'} . " SNVs\n";
		print SUMMARY $family_stats{'snvs_control'} . " present in control(s)\n";
		print SUMMARY $family_stats{'snvs_nocontrol'} . " not seen in control(s)\n";
		print SUMMARY $family_stats{'snvs_nocontrol_nonsilent'} . " non-silent\n";
		
		foreach my $key (sort keys %family_stats)
		{
			print SUMMARY "$family_stats{$key} $key\n" if($key =~ 'snvs_nocontrol_');
		}
		
		print SUMMARY $family_stats{'indels'} . " INDELs\n";
		print SUMMARY $family_stats{'indels_control'} . " present in control(s)\n";
		print SUMMARY $family_stats{'indels_nocontrol'} . " not seen in control(s)\n";
		print SUMMARY $family_stats{'indels_nocontrol_nonsilent'} . " non-silent\n";

		foreach my $key (sort keys %family_stats)
		{
			print SUMMARY "$family_stats{$key} $key\n" if($key =~ 'indels_nocontrol_nonsilent_');
		}
		
		print join("\t", $family, $family_stats{'snvs'}, $family_stats{'snvs_nocontrol'}, $family_stats{'snvs_nocontrol_nonsilent'}, $family_stats{'indels'}, $family_stats{'indels_nocontrol'}, $family_stats{'indels_nocontrol_nonsilent'}) . "\n";


#	}		
	}
	
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
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

################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_variants
{                               # replace with real execution logic.
	my $dir = shift(@_);
	my $self = shift(@_);
	my %variants = ();
	my $post_annotation_file = "$dir/variants/filtered.variants.post_annotation";

	my $input = new FileHandle ($post_annotation_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);		
		my @lineContents = split(/\t/, $line);
		my $variant_type = $lineContents[5];
		my $trv_type = $lineContents[13];
		my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
		
		if($self->tier1_only)
		{
			if(is_tier1($trv_type))
			{
				$variants{$key} = $line;				
			}
		}
		elsif($self->outside_tier1_only)
		{
			if(!is_tier1($trv_type))
			{
				$variants{$key} = $line;				
			}
		}
		else
		{
			$variants{$key} = $line;
		}
	}
	
	close($input);

	return(%variants);
}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub is_tier1
{
	my $trv_type = shift(@_);
	return(1) if($trv_type eq 'frame_shift_del');
	return(1) if($trv_type eq 'frame_shift_ins');
	return(1) if($trv_type eq 'in_frame_del');
	return(1) if($trv_type eq 'in_frame_ins');
	return(1) if($trv_type eq 'missense');
	return(1) if($trv_type eq 'nonsense');
	return(1) if($trv_type eq 'nonstop');
	return(1) if($trv_type eq 'rna');
	return(1) if($trv_type eq 'silent');
#	return(1) if($trv_type eq 'splice_region');
#	return(1) if($trv_type eq 'splice_region_del');
#	return(1) if($trv_type eq 'splice_region_ins');
	return(1) if($trv_type eq 'splice_site');
	return(1) if($trv_type eq 'splice_site_del');
	return(1) if($trv_type eq 'splice_site_ins');

	return(0);
}

1;
