package Genome::Model::Tools::Capture::MutationsFromGroup;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MutationsFromGroup - Retrieve, reformat, and annotate mutations from a somatic-variation group
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	02/09/2012 by D.K.
#	MODIFIED:	03/15/2012 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it
use Genome::Model::Tools::Capture::Helpers qw/iupac_to_base sortByChrPos/;

## Declare global statistics hash ##

my %stats = ();


class Genome::Model::Tools::Capture::MutationsFromGroup {
	is => 'Genome::Model::Tools::Capture',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		group_id		=> { is => 'Text', doc => "ID of somatic-variation model group" , is_optional => 0},
		output_dir	=> { is => 'Text', doc => "An output directory to hold copy number results" , is_optional => 0},
		output_report	=> { is => 'Text', doc => "An output file for a summary report" , is_optional => 1},
		varscan_params	=> { is => 'Text', doc => "Parameters to pass to VarScan copynumber" , is_optional => 1, default => '--min-coverage 20 --min-segment-size 25 --max-segment-size 100'},
		reference	=> 	{ is => 'Text', doc => "Reference to use for bam-readcounts-based filters; defaults to build 37" , is_optional => 1, example_values => ['/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa']},
		test_only	=> { is => 'Text', doc => "Set to 1 to see how the command would run without doing anything" , is_optional => 1},
		run_limit	=> { is => 'Text', doc => "If specified, will only launch this many copy number in parallel jobs" , is_optional => 1},		
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Retrieves, reformat, and annotate mutations from a somatic-variation group"                 
}

sub help_synopsis {
    return <<EOS
This command will retrieve, reformat, and annotate mutations from a model group of somatic-variation models
EXAMPLE:	gmt capture mutations-from-group --group-id 3328 --output-dir somatic_mutations
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
	my $output_dir = $self->output_dir;

	if($self->output_report)
	{
		open(REPORT, ">" . $self->output_report) or die "Can't open outfile: $!\n";
		print REPORT "tumor_sample\tnormal_sample\tsnvs_novel\tsnvs_tier1\tindels_novel\tindels_tier1\tbuild_dir\n";
	}

	## Get the models in each model group ##

	my $model_group = Genome::ModelGroup->get($group_id);
	my @models = $model_group->models;
	my $debug_counter = 0;

	foreach my $model (@models)
	{
		$stats{'models_in_group'}++;
		
		my $model_id = $model->genome_model_id;
		my $model_name = $model->name;
		my $subject_name = $model->subject_name;
		$subject_name = "Model" . $model_id if(!$subject_name);
		my $build_dir = $model->last_succeeded_build_directory;

		## Get normal and tumor model ##
		
		my $normal_model = $model->normal_model;
		my $tumor_model = $model->tumor_model;

		## Get Model IDs ##
		
		my $normal_model_id = $normal_model->id;
		my $tumor_model_id = $tumor_model->id;

		## Get TCGA-Compliant Subject Names ##

		my $normal_subject_name = $normal_model->subject_name;
		my $tumor_subject_name = $tumor_model->subject_name;

		my ($sample_prefix) = split(/\-/, $normal_subject_name);
		$normal_subject_name =~ s/$sample_prefix/TCGA/;
		
		($sample_prefix) = split(/\-/, $tumor_subject_name);
		$tumor_subject_name =~ s/$sample_prefix/TCGA/;

		my @temp = split(/\-/, $tumor_subject_name);
		my $patient_id = join("-", $temp[0], $temp[1], $temp[2]);

		my $tumor_sample = $tumor_model->subject_name;
		my $normal_sample = $normal_model->subject_name;
		my $tumor_model_dir = $tumor_model->last_succeeded_build_directory;
		my $normal_model_dir = $normal_model->last_succeeded_build_directory;
		my $tumor_bam = `ls $tumor_model_dir/alignments/*.bam | head -1`; chomp($tumor_bam);
		my $normal_bam = `ls $normal_model_dir/alignments/*.bam | head -1`; chomp($normal_bam);

		my $output_dir = $self->output_dir . "/" . $tumor_sample . "-" . $normal_sample;
	
		print join("\t", $tumor_sample, $normal_sample) . "\n";
		
		if($self->test_only)
		{

		}
		else
		{
			if($build_dir)
			{
				mkdir($output_dir) if(!(-d $output_dir));
				
				## Get Novel SNVs file ##
				my $novel_snvs_file = "$build_dir/novel/snvs.hq.novel.v2.bed";
				my $tier1_snvs_file = "$build_dir/effects/snvs.hq.novel.tier1.v2.bed";				
				my $tier1_snvs_annotation_file = "$build_dir/effects/snvs.hq.tier1.v1.annotated.top";
				## Load Tier 1 ##
				my %tier1_snvs = $self->load_variants_from_bed($tier1_snvs_file);
				
				## Get Novel Indels file ##
				my $novel_indels_file = "$build_dir/novel/indels.hq.novel.v2.bed";
				my $tier1_indels_file = "$build_dir/effects/indels.hq.novel.tier1.v2.bed";
				my $tier1_indels_annotation_file = "$build_dir/effects/indels.hq.tier1.v1.annotated.top";

				## Load Tier1 ##				
				my %tier1_indels = $self->load_variants_from_bed($tier1_indels_file);
				
				if(-e $novel_snvs_file)
				{
					my $mpileup_file = "$output_dir/snvs.novel.mpileup";

					## Copy over annotation files ##
					system("cp -r $tier1_snvs_annotation_file $output_dir/snvs.novel.tier1.transcript-annotation");

					open(NOVELBED, ">$output_dir/snvs.novel.bed") or die "Can't open outfile: $!\n";
					open(TIER1BED, ">$output_dir/snvs.novel.tier1.bed") or die "Can't open outfile: $!\n";
					open(NOVEL, ">$output_dir/snvs.novel.formatted") or die "Can't open outfile: $!\n";
					open(SCRIPT, ">$output_dir/script-mpileup-snvs.sh") or die "Can't open outfile: $!\n";
					print SCRIPT "#!/gsc/bin/sh\n";
					my $novel_snvs = load_file_to_string($novel_snvs_file);
					my @novel_snvs = split(/\n/, $novel_snvs);
					my $num_snvs = @novel_snvs;
					my $num_snvs_tier1 = 0;

					foreach my $snv (@novel_snvs)
					{
						print NOVELBED "$snv\n";
						my ($chrom, $chr_start, $chr_stop, $alleles) = split(/\t/, $snv);

						## Reformat ##
						$chr_start++ if($chr_start < $chr_stop);
						my ($ref, $cns) = split(/\//, $alleles);
						my $var = iupac_to_base($ref, $cns);
						my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
						
						print NOVEL join("\t", $chrom, $chr_start, $chr_stop, $ref, $var) . "\n";
						my $query = $chrom . ":" . $chr_start . "-" . $chr_stop;
						
						if($tier1_snvs{$key})
						{
							print TIER1BED "$snv\n";
							$num_snvs_tier1++;
						}

						print SCRIPT "samtools mpileup -f " . $self->reference . " -q 1 -r $query $normal_bam $tumor_bam >>$mpileup_file\n";
					}
					
					close(NOVEL);
					close(NOVELBED);
					close(TIER1BED);
					close(SCRIPT);


					
					## PROCESS TIER 1 INDELS ONLY ##

					my $num_indels = my $num_indels_tier1 = 0;



					if(-e $novel_indels_file)
					{
						my %novel_indels = $self->load_variants_from_bed($novel_indels_file);
						
						
						open(NOVELBED, ">$output_dir/indels.novel.bed") or die "Can't open outfile: $!\n";
						open(TIER1BED, ">$output_dir/indels.novel.tier1.bed") or die "Can't open outfile: $!\n";
						open(NOVEL, ">$output_dir/indels.novel.formatted") or die "Can't open outfile: $!\n";
						open(SCRIPT, ">$output_dir/script-mpileup-indels.sh") or die "Can't open outfile: $!\n";
						print SCRIPT "#!/gsc/bin/sh\n";
	
						my %printed = ();
	
						foreach my $indel (sortByChrPos(keys %novel_indels))
						{
							$num_indels++;
							
							if(!$printed{$indel})
							{
								print NOVEL "$indel\n";
								print NOVELBED $novel_indels{$indel} . "\n";
	
								my ($chrom, $chr_start, $chr_stop) = split(/\t/, $indel);
								$chr_start--;
								$chr_stop++;
								my $query = $chrom . ":" . $chr_start . "-" . $chr_stop;
								
								if($tier1_indels{$indel})
								{
									print TIER1BED $novel_indels{$indel} . "\n";
									$num_indels_tier1++;
									print SCRIPT "samtools mpileup -f " . $self->reference . " -q 1 -r $query $normal_bam $tumor_bam >>$mpileup_file\n";
								}

								$printed{$indel} = 1;
							}

						}
						
						close(NOVEL);
						close(NOVELBED);
						close(TIER1BED);
						close(SCRIPT);

						## Get Annotation ##
						system("sort -u $tier1_indels_annotation_file >$output_dir/indels.novel.tier1.transcript-annotation");
						system("gmt capture sort-by-chr-pos --input $output_dir/indels.novel.tier1.transcript-annotation --output $output_dir/indels.novel.tier1.transcript-annotation 2>/dev/null");						
					}
					
					if($self->output_report)
					{
						print REPORT join("\t", $tumor_sample, $normal_sample, $num_snvs, $num_snvs_tier1, $num_indels, $num_indels_tier1, $build_dir) . "\n";
					}


					print join("\t", $num_snvs, $num_snvs_tier1, $num_indels, $num_indels_tier1, $build_dir, $output_dir, $novel_snvs_file) . "\n";
				

#					exit(0);

				}
			}				
		}

#		exit(0);

	}	
	
	foreach my $key (sort keys %stats)
	{
		print $stats{$key} . " " . $key . "\n";
	}
	
	close(REPORT) if($self->output_report);
	
	return 1;
}



#############################################################
# parse_file - parses the file
#
#############################################################

sub load_file_to_string
{
	my $FileName = shift(@_);

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	my $lines = "";
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		$lines .= "\n" if($lines);
		$lines .= $line;
		$lineCounter++;
		
		
	}
	
	close($input);
	
	return($lines);
}



#############################################################
# parse_file - parses the file
#
#############################################################

sub load_variants_from_bed
{
    my $self = shift;
	my $FileName = shift;

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	my %variants = ();
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my ($chrom, $chr_start, $chr_stop, $alleles) = split(/\t/, $line);
		my ($ref, $cns) = split(/\//, $alleles);
		my $var = "";

		if(length($ref) > 1 || length($cns) > 1 || $ref eq "0" || $cns eq "0" || $ref eq "-" || $cns eq "-" || $ref eq '*' || $cns eq '*')
		{
			## INDEL ##
			
			if($cns eq "0" || $cns eq "-" || $cns eq '*')
			{
				## DELETION ##
				$var = "-";
				$chr_start++ if($chr_start < $chr_stop);
			}
			elsif($ref eq "0" || $ref eq "-" || $ref eq '*')
			{
				## INSERTION ##
				$ref = "-";
				$var = $cns;
				$chr_start-- if($chr_start == $chr_stop);
			}
		}
		else
		{
			## SNV ##
			$var = iupac_to_base($ref, $cns);
			$chr_start++ if($chr_start < $chr_stop);
		}
		
		my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
		$variants{$key} = $line;
		
	}
	
	close($input);
	
	return(%variants);
}

1;
