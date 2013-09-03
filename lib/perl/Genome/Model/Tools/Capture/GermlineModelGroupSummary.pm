
package Genome::Model::Tools::Capture::GermlineModelGroupSummary;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ModelGroup - Build Genome Models for Capture Datasets
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	12/09/2009 by D.K.
#	MODIFIED:	12/09/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it
use Genome::Model::Tools::Capture::Helpers 'iupac_to_base';

## Declare global statistics hash ##

my %stats = ();

my %snps = my %snp_hets = my %snp_homs = ();

class Genome::Model::Tools::Capture::GermlineModelGroupSummary {
	is => 'Genome::Model::Tools::Capture',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		group_id		=> { is => 'Text', doc => "ID of model group" , is_optional => 0},
		dbsnp_file		=> { is => 'Text', doc => "A list of dbSNPs in annotation format for concordance calculation" , is_optional => 1},
		roi_file		=> { is => 'Text', doc => "A BED file of ROI targets to restrict the analysis" , is_optional => 1},
		output_file		=> { is => 'Text', doc => "An output file to contain the results" , is_optional => 1},
		output_variants		=> { is => 'Text', doc => "A directory to receive per-sample variant calls" , is_optional => 1},
		output_snp_list		=> { is => 'Text', doc => "A file containing the list of all SNPs called and their MAFs" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Operate on germline model groups"                 
}

sub help_synopsis {
    return <<EOS
Operate on germline model groups
EXAMPLE:	gmt capture somatic-model-group --group-id 3328
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

	my $group_id = $self->group_id;

	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
		print OUTFILE "sample\tbuild_dir\tsnps\tdbsnp\tdbsnp\%\tsnp_hets\tsnp_homs\tnum_indels\tnum_ins\tnum_del\n";
	}

	## Get the models in each model group ##

	my $model_group = Genome::ModelGroup->get($group_id);
	my @models = $model_group->models;
	
	foreach my $model (@models)
	{
		$stats{'models_in_group'}++;
		
		my $model_id = $model->genome_model_id;
		my $model_name = $model->name;
		my $subject_name = $model->subject_name;

		if($model->last_succeeded_build_directory)# && $stats{'models_in_group'} <= 20)
		{
			$stats{'samples_included'}++;
			my $build_dir = $model->last_succeeded_build_directory;
			my $variant_summary = $self->summarize_variants($subject_name, $build_dir, $self);

			if($self->output_file)
			{
				print OUTFILE join("\t", $subject_name, $build_dir, $variant_summary) . "\n";
			}

			print join("\t", $stats{'models_in_group'}, $subject_name, $build_dir, $variant_summary) . "\n";
		}

	}	

	## List All of the Variants ##
	
	if($self->output_snp_list)
	{
		open(OUTFILE, ">" . $self->output_snp_list) or die "Can't open outfile: $!\n";		
	}
	
	foreach my $key (sort byChrPos keys %snps)
	{
		$snp_hets{$key} = 0 if(!$snp_hets{$key});
		$snp_homs{$key} = 0 if(!$snp_homs{$key});
		my $maf_estimate = ($snp_hets{$key} + 2 * $snp_homs{$key}) / ($stats{'samples_included'} * 2);
		$maf_estimate = sprintf("%.3f", $maf_estimate);
		print join("\t", $key, $snp_hets{$key}, $snp_homs{$key}, $maf_estimate) . "\n";
		if($self->output_snp_list)
		{
			print OUTFILE join("\t", $snps{$key}, $snp_hets{$key}, $snp_homs{$key}, $maf_estimate) . "\n";			
		}

	}
	
	if($self->output_snp_list)
	{
		close(OUTFILE);
	}

	return 1;
}



################################################################################################
# Summarize variants - count variants in the ROI
#
################################################################################################

sub summarize_variants
{
    my $self = shift;
	my ($sample_name, $build_dir, $self) = @_;

	my %var_stats = ();
	$var_stats{'num_snps'} = $var_stats{'num_snps_dbsnp'} = $var_stats{'dbsnp_concordance'} = $var_stats{'num_insertions'} = $var_stats{'num_deletions'} = 0;
	$var_stats{'num_snp_hets'} = $var_stats{'num_snp_homs'} = 0;


	my $snp_bed_file = "$build_dir/variants/snvs.hq.bed";
	my $variants_file = "$build_dir/variants/filtered.variants.post_annotation";

	## IF ROI file provided, do a limit operation ##
	
	## Build temp file for positions where readcounts are needed ##
	if($self->roi_file)
	{
		my ($tfh,$temp_path) = Genome::Sys->create_temp_file;
		
		unless($tfh)
		{
			$self->error_message("Unable to create temporary file.");
			die;
		}

		close($tfh);
		
		my $cmd = "java -jar $ENV{GENOME_SW_LEGACY_JAVA}/VarScan/VarScan.jar limit $variants_file --regions-file " . $self->roi_file . " --output-file $temp_path";
		system($cmd);
		
		$variants_file = $temp_path;
	}


	## Parse the variants file ##
	if(-e $variants_file)
	{
		## Get a temp file for storing SNP positions ##
		my ($snpfh,$snp_path) = Genome::Sys->create_temp_file;
		unless($snpfh)
		{
			$self->error_message("Unable to create temporary file.");
			die;
		}
	
		if($self->output_variants)
		{
			open(OUTSNP, ">" . $self->output_variants . "/$sample_name.snp") or die "Can't open outfile: $!\n";
			open(OUTINDEL, ">" . $self->output_variants . "/$sample_name.indel") or die "Can't open outfile: $!\n";
		}


		my $input = new FileHandle ($variants_file);
		my $lineCounter = 0;
		
		while (<$input>)
		{
			chomp;
			my $line = $_;
			
			my @lineContents = split(/\t/, $line);
			my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);
			my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);

			my $var_type = $lineContents[5];
			my $trv_type = $lineContents[13];
	
			if($var_type eq "SNP")
			{
				$var_stats{'num_snps'}++;
				## Save annotation for this SNP ##
				$snps{$key} = $line;
				print $snpfh "$line\n";
				print OUTSNP "$line\n";
			}
			elsif($var_type eq "INS")
			{
				$var_stats{'num_insertions'}++;
				print OUTINDEL "$line\n";
			}
			elsif($var_type eq "DEL")
			{
				$var_stats{'num_deletions'}++;
				print OUTINDEL "$line\n";
			}
		}
	
		close($input);

		close($snpfh);

		
		## Calculate dbSNP concordance ##
		if($self->dbsnp_file)
		{
			my ($dbsnpfh,$dbsnp_path) = Genome::Sys->create_temp_file;
			unless($dbsnpfh)
			{
				$self->error_message("Unable to create temporary file.");
				die;
			}
			close($dbsnpfh);
			
			my $cmd = "java -jar $ENV{GENOME_SW_LEGACY_JAVA}/VarScan/VarScan.jar limit $snp_path --regions-file " . $self->dbsnp_file . " --output-file $dbsnp_path";
			system($cmd);
	
			## Calculate dbSNP concordance ##
			
			$var_stats{'num_snps_dbsnp'} = `cat $dbsnp_path | wc -l`;
			chomp($var_stats{'num_snps_dbsnp'});
			
			## calculate dbSNP concordance ##
			$var_stats{'dbsnp_concordance'} = sprintf("%.2f", $var_stats{'num_snps_dbsnp'} / $var_stats{'num_snps'} * 100) . '%';

			if($self->output_variants)
			{
				my $output_file = $self->output_variants . "/" . $sample_name . ".snp.known";
				system("cp -r $dbsnp_path $output_file");
			}

		}

		## Total Indels ##
		
		$var_stats{'num_indels'} = $var_stats{'num_insertions'} + $var_stats{'num_deletions'};

		## Parse the SNP BED file ##
		($var_stats{'num_snp_hets'}, $var_stats{'num_snp_homs'}) = $self->parse_snp_bed($snp_bed_file);
	}


	my $snp_summary = join("\t", $var_stats{'num_snps'}, $var_stats{'num_snps_dbsnp'}, $var_stats{'dbsnp_concordance'}, $var_stats{'num_snp_hets'}, $var_stats{'num_snp_homs'});
	my $indel_summary = join("\t", $var_stats{'num_indels'}, $var_stats{'num_insertions'}, $var_stats{'num_deletions'});
	return(join("\t", $snp_summary, $indel_summary));
	
}



#############################################################
# parse_file - parses the file
#
#############################################################

sub parse_snp_bed
{
    my $self = shift;
	my $FileName = shift(@_);

	my $num_hets = my $num_homs = 0;
	
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		my ($chrom, $chr_start, $chr_stop, $alleles) = split(/\t/, $line);
		my ($ref, $cns) = split(/\//, $alleles);
		$chr_start++ if($chr_start < $chr_stop);
		
		my $var = iupac_to_base($ref, $cns);
		my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
		
		if(is_homozygous($cns))
		{
			$snp_homs{$key}++;
			$num_homs++;
		}
		else
		{
			$snp_hets{$key}++;
			$num_hets++;
		}
	}

	close($input);

	return($num_hets, $num_homs);
}



#############################################################
# parse_file - parses the file
#
#############################################################

sub parse_file
{
	my $FileName = shift(@_);

	my @lines = ();

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
	}

	close($input);
}


#############################################################
# byChrPos - a sorting subroutine by chrom then pos
#
#############################################################

sub byChrPos
{
    (my $chrom_a, my $pos_a) = split(/\t/, $a);
    (my $chrom_b, my $pos_b) = split(/\t/, $b);

	$chrom_a =~ s/X/23/;
	$chrom_a =~ s/Y/24/;
	$chrom_a =~ s/MT/25/;
	$chrom_a =~ s/M/25/;
	$chrom_a =~ s/[^0-9]//g;

	$chrom_b =~ s/X/23/;
	$chrom_b =~ s/Y/24/;
	$chrom_b =~ s/MT/25/;
	$chrom_b =~ s/M/25/;
	$chrom_b =~ s/[^0-9]//g;

    $chrom_a <=> $chrom_b
    or
    $pos_a <=> $pos_b;
    
#    $chrom_a = 23 if($chrom_a =~ 'X');
#    $chrom_a = 24 if($chrom_a =~ 'Y');
    
}



#############################################################
# is_homozygous - returns 1 if cns code is homozygous
#
#############################################################

sub is_homozygous
{
	my $cns = shift(@_);
	if($cns eq "A" || $cns eq "C" || $cns eq "G" || $cns eq "T")
	{
		return(1);
	}
	return(0);
}

1;
