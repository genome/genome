
package Genome::Model::Tools::Capture::ProcessModels;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ProcessModels - Compare germline reference models to find germline events
#					
#	AUTHOR:		Will Schierding (wschierd@genome.wustl.edu)
#
#	CREATED:	6/19/2009 by W.S.
#	MODIFIED:	6/19/2009 by W.S.
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

class Genome::Model::Tools::Capture::ProcessModels {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		output_dir	=> { is => 'Text', doc => "Output directory for comparison files" , is_optional => 0},
		model_list	=> { is => 'Text', doc => "Text file id,subject_name,build_ids,build_statuses,last_succeeded_build_directory, one per line - space delim" , is_optional => 0},
		regions_file	=> { is => 'Text', doc => "Optional limit to regions file" , is_optional => 1},
		skip_if_output_present => { is => 'Text', doc => "Do not attempt to run pipeline if output present" , is_optional => 1},
		verbose => { is => 'Text', doc => "Display Lots of Output" , is_optional => 1, default => 1},
		use_stable => { is => 'Text', doc => "1 if you want to submit genome-stable version, 0 if you want to use your present dir version" , is_optional => 1, default => 0},
		skip_roi => { is => 'Text', doc => "1 if you want to skip roi filtering (exome), 0 if you want to use the filter" , is_optional => 1, default => 0},
		only_tier_1 => { is => 'Text', doc => "1 if you want to have only tier 1 results (skip ucsc), 0 if you want full tiering" , is_optional => 1, default => 0},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Perform downstream analysis on a list of genome models"                 
}

sub help_synopsis {
    return <<EOS
Perform downstream analysis on a list of genome models.  The list should be tab-delimited with model_id and sample_name as the first two columns.
EXAMPLE:	gmt capture process-models ...
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

	## Get required parameters ##
	my $model_list = $self->model_list;
	my $verbose = $self->verbose;
	my $stable = $self->use_stable;
	my $skip_roi = $self->skip_roi;
	my $only_tier_1 = $self->only_tier_1;
	my $output_dir = "./";
	$output_dir = $self->output_dir if($self->output_dir);
	my $regions_file;
	if($self->regions_file) {
		$regions_file = $self->regions_file ;
	}
	else {
		if ($skip_roi) {
			$regions_file = '/gscmnt/sata424/info/medseq/Freimer-Boehnke/ExomeComparison/WuSpace_2514360.bed'; # dummy roi file = wuspace bed file
		}
		else {
			die "No regions file submitted. Please specify roi file or enable --roi-skip option";
		}
	}

	my $input = new FileHandle ($model_list);
	my $lineCounter = 0;
	my $i = 0;
	while (<$input>)
	{
		$i++;
		chomp;
		my $line = $_;
		$stats{'num_pairs'}++;
		$lineCounter++;
		$line =~ s/\s+/\t/g;
		my ($model_id, $sample_name, $build_id, $build_status, $build_dir, $bam_file) = split(/\t/, $line);
		## get the bam file ##
		unless ($bam_file) {
			my $model = Genome::Model->get($model_id);
			my $build = $model->last_succeeded_build;
			$bam_file = $build->whole_rmdup_bam_file;
		}

		## Establish sample output dir ##
		my $sample_output_dir = $output_dir . "/" . $sample_name;
		mkdir($sample_output_dir) if(!(-d $sample_output_dir));
		if ($verbose) {
			print "$model_id\t$sample_name\t$build_status\t$build_dir\t$bam_file\n";
		}

		my $snp_file = $build_dir . "/snp_related_metrics/snps_all_sequences.filtered";
		my $indel_file = $build_dir . "/snp_related_metrics/indels_all_sequences.filtered";

		if(-e $bam_file && -e $snp_file && -e $indel_file) {
			my $varscan_snps = "";
			$varscan_snps = `cat $sample_output_dir/varScan.output.snp | wc -l` if(-e "$sample_output_dir/varScan.output.snp");
			chomp($varscan_snps) if($varscan_snps);

			my $final_snp_file = "$sample_output_dir/merged.germline.snp.ROI.strandfilter.tier1.out";
			my $final_indel_file = "$sample_output_dir/merged.germline.indel.ROI.strandfilter.tier1.out";

			my $final_snp_file2 = "$sample_output_dir/merged.germline.snp.ROI.strandfilter.tier2.out";
			my $final_indel_file2 = "$sample_output_dir/merged.germline.indel.ROI.strandfilter.tier2.out";
			my $final_snp_file3 = "$sample_output_dir/merged.germline.snp.ROI.strandfilter.tier3.out";
			my $final_indel_file3 = "$sample_output_dir/merged.germline.indel.ROI.strandfilter.tier3.out";
			my $final_snp_file4 = "$sample_output_dir/merged.germline.snp.ROI.strandfilter.tier4.out";
			my $final_indel_file4 = "$sample_output_dir/merged.germline.indel.ROI.strandfilter.tier4.out";

			my $snpexists = 0;
			my $indelexists = 0;
			if (-s $final_snp_file) {
				$snpexists = 1;
			}
			elsif (!$only_tier_1 && (-s $final_snp_file2 || -s $final_snp_file3 || -s $final_snp_file4)) {
				print "No tier 1 SNP File: $sample_name at: $final_snp_file\n";
				$snpexists = 1;
			}
			else {
				print "No tiered SNP File: $sample_name\n";
			}

			if (-s $final_indel_file) {
				$indelexists = 1;
			}
			elsif (!$only_tier_1 && (-s $final_indel_file2 || -s $final_indel_file3 || -s $final_indel_file4)){
				print "No tier 1 Indel File: $sample_name at: $final_indel_file\n";
				$indelexists = 1;
			}
			else {
				print "No tiered Indel File: $sample_name\n";
			}

			if ($snpexists == 1) {
				my $output_completed = `grep -i "Successfully completed" $sample_output_dir/$sample_name.output`;
				unless($output_completed) {
					$snpexists = 0;
					print "LSF Didn't Report Success: $sample_name";
				}
			}

			if($self->skip_if_output_present && $snpexists && $indelexists)
			{
				## Skip because valid output ##
				print "skipped $sample_name for already having valid output\n";
			}
			else
			{
				if($verbose) {
					print "$model_id\t$sample_name\t$build_status\t$build_dir\n";
				}
				my @outfile_list = qw(annotation.germline.indel.ucsc merged.germline.indel merged.germline.indel.ROI.tier4.out merged.germline.snp.ROI samtools.output.indel.formatted varScan.output.snp annotation.germline.indel.unannot-ucsc merged.germline.indel.ROI merged.germline.indel.shared merged.germline.snp.ROI.tier1.out samtools.output.snp.adaptor varScan.output.snp.filter annotation.germline.snp.transcript merged.germline.indel.ROI.tier1.out merged.germline.indel.sniper-only merged.germline.snp.ROI.tier2.out varScan.output.indel varScan.output.snp.formatted annotation.germline.snp.ucsc merged.germline.indel.ROI.tier2.out merged.germline.indel.varscan-only merged.germline.snp.ROI.tier3.out varScan.output.indel.filter varScan.output.snp.variants annotation.germline.indel.transcript annotation.germline.snp.unannot-ucsc merged.germline.indel.ROI.tier3.out merged.germline.snp merged.germline.snp.ROI.tier4.out varScan.output.indel.formatted $sample_name.out $sample_name.err merged.germline.snp.ROI.tier1.out.strandfilter.readcounts merged.germline.snp.ROI.tier1.out.strandfilter_filtered merged.germline.snp.ROI.tier1.out.strandfilter GATK.output.indel GATK.output.indel.vcf GATK.output.indel.bed GATK.output.indel.formatted GATK.output.indel.adaptor merged.germline.snp.ROI.tier1.out.dbsnp merged.germline.indel.ROI.tier1.out.strandfilter merged.germline.indel.ROI.tier1.out.strandfilter.readcounts merged.germline.indel.ROI.tier1.out.strandfilter_filtered merged.germline.ROI.tier1.out.maf merged.germline.snp.ROI.strandfilter.tier1.out merged.germline.indel.ROI.strandfilter.tier1.out annotation.germline.indel.strandfilter.unannot-ucsc annotation.germline.snp.strandfilter.unannot-ucsc merged.germline.snp.ROI.strandfilter.readcounts merged.germline.snp.ROI.strandfilter merged.germline.snp.ROI.strandfilter_filtered annotation.germline.snp.strandfilter.ucsc annotation.germline.indel.strandfilter.ucsc annotation.germline.snp.strandfilter.transcript GATK.ug.output.indel GATK.ug.output.indel.idx GATK.ug.output.indel.adaptor merged.germline.indel.gatk-only merged.germline.indel.ROI.strandfilter.readcounts merged.germline.indel.ROI.strandfilter_filtered merged.germline.indel.ROI.strandfilter annotation.germline.indel.strandfilter.transcript);

				foreach my $file (@outfile_list) {
					my $del_file = "$sample_output_dir/$file";
					if (-e $del_file) {
						unlink("$del_file");
					}
				}
				my $cmd;
				if ($stable) {
					$cmd = "perl -I /gsc/scripts/opt/genome/current/pipeline/lib/perl/ `which gmt` germline capture-bams --build-id $build_id --germline-bam-file $bam_file --filtered-indelpe-snps $snp_file --indels-all-sequences-filtered $indel_file --data-directory $sample_output_dir --skip-roi $skip_roi --only-tier-1 $only_tier_1 --only-tier-1-indel $only_tier_1 --regions-file $regions_file";

				}
				else {
					$cmd = "gmt germline capture-bams --build-id $build_id --germline-bam-file $bam_file --filtered-indelpe-snps $snp_file --indels-all-sequences-filtered $indel_file --data-directory $sample_output_dir --skip-roi $skip_roi --only-tier-1 $only_tier_1 --only-tier-1-indel $only_tier_1 --regions-file $regions_file";
				}
				if($verbose) {
					print "$cmd\n";
				}
				my $job_name = "$sample_output_dir/";
				my $output_name = "$sample_output_dir/$sample_name.output";
				my $error_name = "$sample_output_dir/$sample_name.err";
				my $output_name_bak = "$sample_output_dir/$sample_name.output.bak";
				my $error_name_bak = "$sample_output_dir/$sample_name.err.bak";
				if (-e $output_name) {
					my $cmd_mv = "mv $output_name $output_name_bak";
					system($cmd_mv);
				}
				if (-e $error_name) {;
					my $cmd_mv = "mv $error_name $error_name_bak";
					system($cmd_mv);
				}
				system("bsub -u wschierd\@genome.wustl.edu -q apipe -R\"select[type==LINUX64 && model != Opteron250 && mem>4000] rusage[mem=4000]\" -M 4000000 -J $job_name -o $output_name -e $error_name \"$cmd\"");
				sleep(1);
			}
		}
		else {
			print "-e bam_file && -e snp_file && -e indel_file failed";
			exit;
		}
		my $longqueue_pending=`bjobs -q long | grep PEND | wc -l`; chomp $longqueue_pending;
		my $apipequeue_pending=`bjobs -q apipe | grep PEND | wc -l`; chomp $apipequeue_pending;
		if ($longqueue_pending >= 75) {
			sleep(600);
		}
		elsif ($apipequeue_pending >= 50) {
			sleep(600);
		}
	}

	close($input);
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


1;

