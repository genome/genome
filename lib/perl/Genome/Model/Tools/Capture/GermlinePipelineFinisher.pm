
package Genome::Model::Tools::Capture::GermlinePipelineFinisher;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# GermlinePipelineFinisher - Merge a set of maf files into one giant one -- for GERMLINE events, also serves as a testing ground for new ideas to be placed into the germline pipeline
#					
#	AUTHOR:		Will Schierding (wschierd@genome.wustl.edu)
#
#	CREATED:	09/29/2010 by W.S.
#	MODIFIED:	09/29/2010 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

## Declare global statistics hash ##

my %stats = ();

class Genome::Model::Tools::Capture::GermlinePipelineFinisher {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		data_dir	=> { is => 'Text', doc => "Base Data Directory i.e. /gscmnt/sata424/info/medseq/Freimer-Boehnke/Analysis-1033Samples/ " , is_optional => 0},
		project_name	=> { is => 'Text', doc => "Name of the project i.e. ASMS" , is_optional => 0},
		sample_list	=> { is => 'Text', doc => "File of sample names to include, 1 name per line, no headers" , is_optional => 0},
		model_list	=> { is => 'Text', doc => "Same as input to germline pipeline, no headers, (space or tab delim) model_id, sample_name, build_id, build_status, build_dir" , is_optional => 0},
		center 		=> { is => 'Text', doc => "Genome center name" , is_optional => 1, default => "genome.wustl.edu"},
		build 		=> { is => 'Text', doc => "Reference genome build = either '36' or '37'" , is_optional => 0, default => "36"},
		lsf_user_id	=> { is => 'Text', doc => "gsc username to supply to lsf for e-mails" , is_optional => 1, default => "wschierd"},
		sequence_phase	=> { is => 'Text', doc => "Sequencing phase" , is_optional => 1, default => "4"},
		sequence_source	=> { is => 'Text', doc => "Sequence source" , is_optional => 1, default => "Capture"},
		sequencer	=> { is => 'Text', doc => "Sequencing platform name" , is_optional => 1, default => "IlluminaGAIIx"},
		make_vcf	=> { is => 'Text', doc => "Make a vcf file for each sample" , is_optional => 1, default => "0"},
		make_maf	=> { is => 'Text', doc => "Make a single maf file for each sample, named by output-file" , is_optional => 1, default => "0"},
		output_file	=> { is => 'Text', doc => "Name of MAF File" , is_optional => 1},
	        skip_if_output_present => { is => 'Boolean', is_optional => 1, is_input => 1, default => 0, doc => 'enable this flag to shortcut through vcf file making if the output_file is already present.',},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Merge a set of maf files into one giant one -- for GERMLINE projects, also serves as a testing ground for new ideas to be placed into the germline pipeline"                 
}

sub help_synopsis {
    return <<EOS
Merge a set of maf files into one giant one -- for GERMLINE events, also serves as a testing ground for new ideas to be placed into the germline pipeline
EXAMPLE:	gmt capture germline-pipeline-finisher
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

	my $data_dir = $self->data_dir;
	my $project_name = $self->project_name;
	my $sample_list_file = $self->sample_list;
	my $model_list_file = $self->model_list;	

	my $center = $self->center;
	my $build = $self->build;
	my $sequence_phase = $self->sequence_phase;
	my $sequence_source = $self->sequence_source;
	my $sequencer = $self->sequencer;

	my $user = $self->lsf_user_id;
	my $skipifoutputpresent = $self->skip_if_output_present;

	unless ($self->make_maf || $self->make_vcf) {
		die "Must specify maf or vcf (or can call both) to make";
	}

	my $return = 1;

	my %sample_list;
	my $sample_input = new FileHandle ($sample_list_file);
	while (my $sample = <$sample_input>) {
		chomp($sample);
		$sample_list{$sample}++;
	}
	my $sample_count = 0;
	foreach (sort keys %sample_list) {
		$sample_count++;
	}
	print "Sample List Loaded, $sample_count Samples in List\n";

	my %model_hash;
	my $model_input = new FileHandle ($model_list_file);
	while (my $line = <$model_input>) {
		chomp($line);
		$line =~ s/\s+/\t/g;
		my ($model_id, $sample_name, $build_id, $build_status, $build_dir, $bam_file) = split(/\t/, $line);
		if ($bam_file) {
			if (exists $sample_list{$sample_name}) {
				$model_hash{$sample_name} = "$model_id\t$build_id\t$build_dir\t$bam_file";
			}
			else {
				print "Sample in Model list but not in sample list, skipped sample: $sample_name\n";
			}
		}
		else {
			print "No bam file specified for sample: $sample_name\n";
			if (exists $sample_list{$sample_name}) {
				$model_hash{$sample_name} = "$model_id\t$build_id\t$build_dir";
			}
			else {
				print "Sample in Model list but not in sample list, skipped sample: $sample_name\n";
			}
		}
	}
	my $model_count = 0;
	foreach (sort keys %model_hash) {
		$model_count++;
	}
	print "Model List Loaded, $model_count Models in List\n";

	if ($self->make_maf) {
		unless ($self->output_file) { die "Must supply output-file name for maf file";}
		my $output_file = $self->output_file;
		## Open the outfile ##
		my $outfile = $data_dir . $output_file;
		open(OUTFILE, ">$outfile") or die "Can't open output file: $!\n";
		my $firstfile = 1;
		foreach my $sample_name (sort keys %model_hash) {
			my $sample_output_dir = $data_dir . "/" . $sample_name;
			my ($model_id, $build_id, $build_dir, $bam_file) = split(/\t/, $model_hash{$sample_name});
			my $maf_file = $sample_output_dir . "/merged.germline.ROI.tier1.out.maf";
			if (-e $maf_file) {
				my $input = new FileHandle ($maf_file);
				my $header = <$input>;
				if ($firstfile) {
					print OUTFILE "$header";
					$firstfile = 0;
				}
				while (<$input>) {
					my $line = $_;
					print OUTFILE "$line";
				}
				close($input);
			}
			else {
				print "Sample $sample_name does not have maf file: $maf_file\n";
			}
		}
	}

	if ($self->make_vcf) {
		foreach my $sample_name (sort keys %model_hash) {
			my $sample_output_dir = $data_dir . "/" . $sample_name . "/";
			my $raw_merged_vcf = "$sample_output_dir/merged.variants.vcf";
			my $outfile = $sample_output_dir . "merged.germline.ROI.tier1.vcf";
			my $annotation_infile = $sample_output_dir . "merged.germline.ROI.tier1.vcf.to_annotation";
			my $annotation_outfile = $sample_output_dir . "merged.germline.ROI.tier1.vcf.from_annotation";
			my ($model_id, $build_id, $build_dir, $bam_file) = split(/\t/, $model_hash{$sample_name});
			my $snp_file = $build_dir . "/snp_related_metrics/snps_all_sequences.filtered";
			my $indel_file = $build_dir . "/snp_related_metrics/indels_all_sequences.filtered";
			unless ($bam_file) {
				my $model = Genome::Model->get($model_id);
				my $build = $model->last_succeeded_build;
				$bam_file = $build->whole_rmdup_bam_file;
			}
			my $varscan_snv_file = $sample_output_dir.'/varScan.output.snp';
			my $samtools_snv_file = $sample_output_dir.'/samtools.output.snp.adaptor';
			my $varscan_indel_file = $sample_output_dir.'/varScan.output.indel.formatted';
			my $samtools_indel_file = $sample_output_dir.'/samtools.output.indel.formatted';
			my $gatk_indel_file = $sample_output_dir.'/GATK.output.indel.vcf';
			my $gatk_dindel_indel_file = $sample_output_dir.'/GATK.ug.output.indel';

			my $pre_roi_file = $sample_output_dir.'/merged.germline.snp';
			my $roi_file = $sample_output_dir.'/merged.germline.snp.ROI';
			my $pre_roi_file_indel = $sample_output_dir.'/merged.germline.indel';
			my $roi_file_indel = $sample_output_dir.'/merged.germline.indel.ROI';

			my $snv_filter_fail;
			my $indel_filter;
			my $indel_filter_fail;
			my $snv_annotation;
			my $indel_annotation;
			my $indel;
			my $varfile;

			my $snv_filter = $sample_output_dir.'/merged.germline.snp.ROI.tier1.out.strandfilter';
			my $dbsnp = $sample_output_dir.'/merged.germline.snp.ROI.tier1.out.dbsnp';
			if (-e $snv_filter) {
				$snv_filter = $sample_output_dir.'/merged.germline.snp.ROI.tier1.out.strandfilter';
				$snv_filter_fail = $sample_output_dir.'/merged.germline.snp.ROI.tier1.out.strandfilter_filtered';
				$indel_filter = $sample_output_dir.'/merged.germline.indel.ROI.tier1.out.strandfilter';
				$indel_filter_fail = $sample_output_dir.'/merged.germline.indel.ROI.tier1.out.strandfilter_filtered';
				$indel_annotation = $sample_output_dir.'/annotation.germline.indel.transcript';
				$snv_annotation = $sample_output_dir.'/annotation.germline.snp.transcript';
				if (-e $indel_filter) {
					$indel = $indel_filter;
				}
				else {
					$indel = $sample_output_dir.'/merged.germline.indel.ROI.tier1.out';
				}
				$varfile = $snv_filter;#$sample_output_dir.'/merged.germline.snp.ROI.tier1.out';
			}
			else {
				$snv_filter = $sample_output_dir.'/merged.germline.snp.ROI.strandfilter';
				$snv_filter_fail = $sample_output_dir.'/merged.germline.snp.ROI.strandfilter_filtered';
				$indel_filter = $sample_output_dir.'/merged.germline.indel.ROI.strandfilter';
				$indel_filter_fail = $sample_output_dir.'/merged.germline.indel.ROI.strandfilter_filtered';
				$indel_annotation = $sample_output_dir.'/annotation.germline.indel.strandfilter.transcript';
				$snv_annotation = $sample_output_dir.'/annotation.germline.snp.strandfilter.transcript';
				$indel = $sample_output_dir.'/merged.germline.indel.ROI.strandfilter.tier1.out';
				$varfile = $sample_output_dir.'/merged.germline.snp.ROI.strandfilter.tier1.out';
			}

			#create VCF file from varscan snps (unformatted file)
			my $cmd_varsnp = "gmt vcf vcf-maker-varscan --output-file $varscan_snv_file.vcf --varscan-file $varscan_snv_file --type snv --sample-id $sample_name --dbsnp-file ~/gscmnt/sata921/info/medseq/cmiller/annotations/snp130.txt";

			#create VCF file from varscan indels (formatted file)
			my $cmd_varindel = "gmt vcf vcf-maker-varscan --output-file $varscan_indel_file.vcf --varscan-file $varscan_indel_file --type indel --sample-id $sample_name --dbsnp-file /gscmnt/sata921/info/medseq/cmiller/annotations/snp130.txt";

			#merge the two
			my $cmd_varmerge = "gmt vcf vcf-merge --output-file $sample_output_dir/varScan.snp.indel.vcf --vcf-files $varscan_snv_file.vcf,$varscan_indel_file.vcf";

			#create VCF file from samtools snps
			my $cmd_samsnp = "gmt vcf vcf-maker-samtools --output-file $samtools_snv_file.vcf --samtools-file $samtools_snv_file --type snv --sample-id $sample_name --dbsnp-file /gscmnt/sata921/info/medseq/cmiller/annotations/snp130.txt";

			#create VCF file from samtools indels
			my $cmd_samindel = "gmt vcf vcf-maker-samtools --output-file $samtools_indel_file.vcf --samtools-file $samtools_indel_file --type indel --sample-id $sample_name --dbsnp-file /gscmnt/sata921/info/medseq/cmiller/annotations/snp130.txt";

			#merge the two
			my $cmd_sammerge = "gmt vcf vcf-merge --output-file $sample_output_dir/samtools.snp.indel.vcf --vcf-files $samtools_snv_file.vcf,$samtools_indel_file.vcf";

			my $filelist;
			my $source_ids;
			#now, do a three-way merge and label the sources
			if (-e "$sample_output_dir/varScan.snp.indel.vcf") {
				$filelist .= "$sample_output_dir/varScan.snp.indel.vcf,";
				$source_ids .= "varscan,";
			}
			if (-e "$sample_output_dir/samtools.snp.indel.vcf") {
				$filelist .= "$sample_output_dir/samtools.snp.indel.vcf,";
				$source_ids .= "samtools,";
			}
			if (-e "$gatk_indel_file") {
				$filelist .= "$gatk_indel_file,";
				$source_ids .= "GATK_igv2,";
			}
			if (-e "$gatk_dindel_indel_file") {
				$filelist .= "$gatk_dindel_indel_file,";
				$source_ids .= "GATK_ug,";
			}
			$filelist =~ s/,$//;
			$source_ids =~ s/,$//;
			my $cmd_allthree = "gmt vcf vcf-merge --output-file $raw_merged_vcf --vcf-files $filelist --source-ids \"$source_ids\"";

#			my $cmd = "bsub -q apipe -o /gscuser/wschierd/Deleteme/$sample_name.out -e /gscuser/wschierd/Deleteme/$sample_name.err -R\"select[type==LINUX64 && model != Opteron250 && mem>4000] rusage[mem=4000]\" -M 4000000 \"$base_cmd\"";

#			my $wc = `wc -l $raw_merged_vcf`;
#			(my $wc2) = $wc =~ m/(\d+)\s+/;
			unless ($skipifoutputpresent && -s $raw_merged_vcf && "1" && "1") { #extra &&&& so that coloring works in gedit
				my $bsub = 'bsub -q apipe -R "select[model!=Opteron250 && type==LINUX64 && mem>8000 && tmp>10000] rusage[mem=8000, tmp=10000]" -M 8000000 '."-u $user ";
				my $jobid1 = `$bsub -J varsnpvcf \'$cmd_varsnp\'`;
				   $jobid1=~/<(\d+)>/;
				   $jobid1= $1;
				   print "$jobid1\n";
				my $jobid2 = `$bsub -J varindelvcf \'$cmd_varindel\'`;
				   $jobid2=~/<(\d+)>/;
				   $jobid2= $1;
				   print "$jobid2\n";
				my $jobid3 = `$bsub -J samsnpvcf \'$cmd_samsnp\'`;
				   $jobid3=~/<(\d+)>/;
				   $jobid3= $1;
				   print "$jobid3\n";
				my $jobid4 = `$bsub -J samindelvcf \'$cmd_samindel\'`;
				   $jobid4=~/<(\d+)>/;
				   $jobid4= $1;
				   print "$jobid4\n";
				my $jobid5 = `$bsub -J varmergevcf -w \'ended($jobid1) && ended($jobid2)\' \'$cmd_varmerge\'`;
				   $jobid5=~/<(\d+)>/;
				   $jobid5= $1;
				   print "$jobid5\n";
				my $jobid6 = `$bsub -J sammergevcf -w \'ended($jobid3) && ended($jobid4)\' \'$cmd_sammerge\'`;
				   $jobid6=~/<(\d+)>/;
				   $jobid6= $1;
				   print "$jobid6\n";
				my $jobid7 = `$bsub -J finalvcf -w \'ended($jobid5) && ended($jobid6)\' \'$cmd_allthree\'`;
				   $jobid7=~/<(\d+)>/;
				   $jobid7= $1;
				   print "$jobid7\n";
			}
			if (-s $raw_merged_vcf && "1" && "1") {
				unless ($skipifoutputpresent && -s $outfile  && 1 && 1) { #extra &&&& so that coloring works in gedit
					open(OUTFILE, ">$outfile") or die "Can't open output file: $!\n";
					my %final_snvs;
					my $snv_input = new FileHandle ($varfile);
					while (my $line = <$snv_input>) {
						chomp($line);
						my ($chr, $start, $stop, $ref, $var, @everything_else) = split(/\t/, $line);
						my $matcher = "$chr\t$start";
						$final_snvs{$matcher}++;
					}
					my %final_indels;
					my $indel_input = new FileHandle ($indel);
					while (my $line = <$indel_input>) {
						chomp($line);
						my ($chr, $start, $stop, $ref, $var, @everything_else) = split(/\t/, $line);
						my $matcher = "$chr\t$start";
						$final_indels{$matcher}++;
					}

					my %filter_reason_snvs;
					if (-s ($snv_filter_fail) && 1) {
						$snv_input = new FileHandle ($snv_filter_fail);
						while (my $line = <$snv_input>) {
							chomp($line);
							my ($chr, $start, $stop, $ref, $var, @everything_else) = split(/\t/, $line);
							my $length = @everything_else;
							my $filter_reason = $everything_else[$length-1];
							my $matcher = "$chr\t$start";
							$filter_reason_snvs{$matcher} = $filter_reason;
						}
					}
					my %filter_reason_indels;
					if (-s ($indel_filter_fail) && 1) {
						$indel_input = new FileHandle ($indel_filter_fail);
						while (my $line = <$indel_input>) {
							chomp($line);
							my ($chr, $start, $stop, $ref, $var, @everything_else) = split(/\t/, $line);
							my $length = @everything_else;
							my $filter_reason = $everything_else[$length-1];
							my $matcher = "$chr\t$start";
							$filter_reason_indels{$matcher} = $filter_reason;
						}
					}

					my %roi_filter;
					$snv_input = new FileHandle ($pre_roi_file);
					while (my $line = <$snv_input>) {
						chomp($line);
						my ($chr, $start, $stop, $ref, $var, @everything_else) = split(/\t/, $line);
						my $matcher = "$chr\t$start";
						$roi_filter{$matcher}++;
					}
					$snv_input = new FileHandle ($roi_file);
					while (my $line = <$snv_input>) {
						chomp($line);
						my ($chr, $start, $stop, $ref, $var, @everything_else) = split(/\t/, $line);
						my $matcher = "$chr\t$start";
						delete $roi_filter{$matcher};
					}

					my %roi_filter_indel;
					$indel_input = new FileHandle ($pre_roi_file_indel);
					while (my $line = <$indel_input>) {
						chomp($line);
						my ($chr, $start, $stop, $ref, $var, @everything_else) = split(/\t/, $line);
						my $matcher = "$chr\t$start";
						$roi_filter_indel{$matcher}++;
					}
					$indel_input = new FileHandle ($roi_file_indel);
					while (my $line = <$indel_input>) {
						chomp($line);
						my ($chr, $start, $stop, $ref, $var, @everything_else) = split(/\t/, $line);
						my $matcher = "$chr\t$start";
						delete $roi_filter_indel{$matcher};
					}

					my $raw_vcf_input = new FileHandle ($raw_merged_vcf);
					while (my $line = <$raw_vcf_input>) {
						if ($line =~ m/#/) {
							print OUTFILE "$line";
							next;
						}
						chomp($line);
						my (@vcf_line) = split(/\t/, $line);
#						my ($chr, $pos, $dbsnp, $ref, $var, $score, $filter, $info, @everything_else) = split(/\t/, $line);
						my $matcher = $vcf_line[0]."\t".$vcf_line[1];
						if ($vcf_line[7] =~ m/SNP/i || $vcf_line[7] =~ m/SNV/i) {
							my $newline;
							if (defined $filter_reason_snvs{$matcher}) {
								$vcf_line[6] = $filter_reason_snvs{$matcher};
								$newline = join("\t",@vcf_line);
							}
							elsif (defined $roi_filter{$matcher}) {
								$vcf_line[6] = "OUTSIDE_ROI";
								$newline = join("\t",@vcf_line);
							}
							elsif (defined $final_snvs{$matcher}) {
								$vcf_line[6] = "PASS";
								$newline = join("\t",@vcf_line);
							}
							else {
								$vcf_line[6] = "FILTER_LIKELY_OUTSIDE_TIER1";
								$newline = join("\t",@vcf_line);
							}
							print OUTFILE "$newline\n";
						}
						elsif ($vcf_line[7] =~ m/INDEL/i) {
							my $newline;
							if (defined $filter_reason_indels{$matcher}) {
								$vcf_line[6] = $filter_reason_indels{$matcher};
								$newline = join("\t",@vcf_line);
							}
							elsif (defined $roi_filter_indel{$matcher}) {
								$vcf_line[6] = "OUTSIDE_ROI";
								$newline = join("\t",@vcf_line);
							}
							elsif (defined $final_indels{$matcher}) {
								$vcf_line[6] = "PASS";
								$newline = join("\t",@vcf_line);
							}
							else {
								$vcf_line[6] = "FILTER_LIKELY_OUTSIDE_TIER1";
								$newline = join("\t",@vcf_line);
							}
							print OUTFILE "$newline\n";
						}
					}
				}
			}
			if (-s $outfile && "1" && "1") {
				unless ($skipifoutputpresent && -s $annotation_outfile  && "1" && "0") { #extra &&&& so that coloring works in gedit
					open(ANNOT, ">$annotation_infile") or die "Can't open output file: $!\n";
					my $vcf_input = new FileHandle ($outfile);
					while (my $line = <$vcf_input>) {
						if ($line =~ m/^#/) {next;}
						chomp($line);
						my ($chrom, $position, $dbsnp, $ref, $var, $confidence, $filter, @stuff) = split(/\t/, $line);
						unless ($filter) { print "This line will fail to have a filter status:\nsample:$sample_name\nline:$line\n";}
						if ($filter eq 'PASS') {
							if ($line =~ m/SNP/) {
								my (@snps) = split(/,/, $var);
								foreach my $variant (@snps) {
									print ANNOT "$chrom\t$position\t$position\t$ref\t$variant\n";
								}
							}
###INDELS ARENT ANNOTATED CORRECTLY BECAUSE OF BAD POSITION
							elsif ($line =~ m/INDEL/) {
								$ref =~ s/^.//;
								$var =~ s/^.//;
								if ($ref eq '') {
									$ref = 0;
								}
								if ($var eq '') {
									$var = 0;
								}
								print ANNOT "$chrom\t$position\t$position\t$ref\t$var\n";
							}
						}
					}
					close(ANNOT);
					
					my $bsub = 'bsub -u wschierd@genome.wustl.edu -J vcf_annot -R "select[type==LINUX64 && mem>8000] rusage[mem=8000]" -M 8000000 ';
					my $annot_cmd_b36 = $bsub."\'perl -I ~/genome-stable/ `which gmt` annotate transcript-variants --variant-file $annotation_infile --output-file $annotation_outfile --annotation-filter top --reference-transcripts NCBI-human.combined-annotation/54_36p_v3\'";
					my $annot_cmd_b37 = $bsub."\'perl -I ~/genome-stable/ `which gmt` annotate transcript-variants --variant-file $annotation_infile --output-file $annotation_outfile --annotation-filter top --reference-transcripts NCBI-human.combined-annotation/58_37c_v2\'";
					if ($build =~ m/36/) {
						system($annot_cmd_b36);
					}
					else {
						system($annot_cmd_b37);
					}
				}
			}



#			$return = Genome::Sys->shellcmd(
#	                           cmd => "$cmd",
#	                           output_files => [$outfile],
#	                           skip_if_output_is_present => 0,
#	                       );
#			unless($return) { 
#				$self->error_message("Failed to execute Vcf Maker: Returned $return");
#				die $self->error_message;
#			}
		}
	}

	return $return;

}




################################################################################################
# SUBS
#
################################################################################################


#############################################################
# IUPAC to base - convert IUPAC code to variant base
#
#############################################################

sub iupac_to_base
{
	(my $allele1, my $allele2) = @_;
	
	return($allele2) if($allele2 eq "A" || $allele2 eq "C" || $allele2 eq "G" || $allele2 eq "T");
	
	if($allele2 eq "M")
	{
		return("C") if($allele1 eq "A");
		return("A") if($allele1 eq "C");
	}
	elsif($allele2 eq "R")
	{
		return("G") if($allele1 eq "A");
		return("A") if($allele1 eq "G");		
	}
	elsif($allele2 eq "W")
	{
		return("T") if($allele1 eq "A");
		return("A") if($allele1 eq "T");		
	}
	elsif($allele2 eq "S")
	{
		return("C") if($allele1 eq "G");
		return("G") if($allele1 eq "C");		
	}
	elsif($allele2 eq "Y")
	{
		return("C") if($allele1 eq "T");
		return("T") if($allele1 eq "C");		
	}
	elsif($allele2 eq "K")
	{
		return("G") if($allele1 eq "T");
		return("T") if($allele1 eq "G");				
	}	
	
	return($allele2);
}

#############################################################
# ParseBlocks - takes input file and parses it
#
#############################################################

sub trv_to_mutation_type
{
	my $trv_type = shift(@_);
	
	return("Missense_Mutation") if($trv_type eq "missense");	
	return("Nonsense_Mutation") if($trv_type eq "nonsense" || $trv_type eq "nonstop");	
	return("Silent") if($trv_type eq "silent");		
	return("Splice_Site_SNP") if($trv_type eq "splice_site");
	return("Splice_Site_Indel") if($trv_type eq "splice_site_del");		
	return("Splice_Site_Indel") if($trv_type eq "splice_site_ins");		
	return("Frame_Shift_Del") if($trv_type eq "frame_shift_del");		
	return("Frame_Shift_Ins") if($trv_type eq "frame_shift_ins");		
	return("In_Frame_Del") if($trv_type eq "in_frame_del");		
	return("In_Frame_Ins") if($trv_type eq "in_frame_ins");		
	return("RNA") if($trv_type eq "rna");		

	warn "Unknown mutation type $trv_type\n";
	return("Unknown");
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub byChrPos
{
	my ($chrom_a, $pos_a) = split(/\t/, $a);
	my ($chrom_b, $pos_b) = split(/\t/, $b);
	
	$chrom_a =~ s/X/23/;
	$chrom_a =~ s/Y/24/;
	$chrom_a =~ s/MT/25/;
	$chrom_a =~ s/[^0-9]//g;

	$chrom_b =~ s/X/23/;
	$chrom_b =~ s/Y/24/;
	$chrom_b =~ s/MT/25/;
	$chrom_b =~ s/[^0-9]//g;

	$chrom_a <=> $chrom_a
	or
	$pos_a <=> $pos_b;
}

1;
