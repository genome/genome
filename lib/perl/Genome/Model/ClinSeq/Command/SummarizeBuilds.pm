package Genome::Model::ClinSeq::Command::SummarizeBuilds;

#Written by Malachi Griffith

use strict;
use warnings;
use Genome;

use Genome::Utility::List;
use Genome::Utility::Text;

use Cwd;

class Genome::Model::ClinSeq::Command::SummarizeBuilds {
    is => 'Command::V2',
    has_input => [
        builds => {
              is => 'Genome::Model::Build::ClinSeq',
              is_many => 1,
              shell_args_position => 1,
              require_user_verify => 0,
              doc => 'clinseq build to summarize',
        },
        outdir => {
              is => 'FilesystemPath',
              doc => 'Directory where output files will be written',
        },
        skip_lims_reports => {
              is => 'Number',
              is_optional => 1,
              doc => 'Use this option to skip LIMS use of illumina_info run/lane/library report tool',
        },
        log_file => {
              is => 'FilesystemPath',
              is_optional => 1,
              doc => 'All status messages go to the specified log file',
        },
        normal_sample_common_names => {
              #TODO: Is there a better way to determine which samples are 'normal'?
              is => 'Text',
              default => 'normal',
              doc => 'The possible sample common names used in the database to specify a Normal sample',
        },
        tumor_sample_common_names => {
              #TODO: Is there a better way to determine which samples are 'tumor'?
              is => 'Text',
              default => 'tumor|met|post treatment|recurrence met|pre-treatment met|pin lesion|relapse',
              doc => 'The possible sample common names used in the database to specify a Tumor sample',
        },
    ],
    doc => 'summarize the inputs of a clinseq build (models/builds, processing profiles, etc.)',
};

sub help_synopsis {
    return <<EOS

genome model clin-seq summarize-builds --outdir=/tmp/  119971814

genome model clin-seq summarize-builds --outdir=/tmp/  --skip-lims-reports=1  119971814

genome model clin-seq summarize-builds --outdir=/tmp/  id=119971814

genome model clin-seq summarize-builds --outdir=/tmp/  model.id=2882726707

genome model clin-seq summarize-builds --outdir=/tmp/  "model.name='ClinSeq - ALL1 - (Nov. 2011 PP) - v2'"

genome model clin-seq summarize-builds --outdir=/tmp/  'id in [119971814,119971932]'

EOS
}

sub help_detail {
    return <<EOS
Summarize the status and key metrics for one or more clinseq models.

(put more content here)
EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);

  unless (-e $self->outdir && -d $self->outdir) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['outdir'],
	                                          desc => "Outdir: " . $self->outdir . " not found or not a directory",
                                          );
  }
  return @errors;
}

sub execute {
  my $self = shift;
  my @builds = $self->builds;
  my $outdir = $self->outdir;
  my $log_file = $self->log_file;

  unless (-e $outdir) {
    Genome::Sys->create_directory($outdir);
  }

  if ($log_file) {
    $self->queue_status_messages(1);
  }

  unless ($outdir =~ /\/$/){
    $outdir .= "/";
  }

  my $clinseq_build_count = scalar(@builds);
  for my $clinseq_build (@builds) {

    #If there is more than one clinseq build supplied.  Create sub-directories for each
    my $build_outdir;
    if ($clinseq_build_count > 1){
      $build_outdir = $outdir . $clinseq_build->id . "/";
      Genome::Sys->create_directory($build_outdir);
    }else{
      $build_outdir = $outdir;
    }

    $self->summarize_clinseq_build($clinseq_build, $build_outdir);
  }

  $self->status_message("\n\n");

  my @output = $self->status_messages();
  if ($log_file) {
      my $log = IO::File->new(">$log_file");
      $log->print(join("\n", @output));
      $log->close;

      $self->queue_status_messages(0);
      $self->status_message("Log file written to $log_file\n");
  }

  return 1;
}

sub summarize_clinseq_build {
    my ($self, $clinseq_build, $build_outdir) = @_;

    #Store grand summary values that apply to the entire ClinSeq build
    my $summary_dir = $build_outdir . "summary/";
    Genome::Sys->create_directory($summary_dir);

    my $stats_file = $summary_dir . "Stats.tsv";
    open (my $stats_fh, ">$stats_file") || die "\n\nCould not open output file: $stats_file\n\n";
    print $stats_fh "Question\tAnswer\tData_Type\tAnalysis_Type\tStatistic_Type\tExtra_Description\n";

    my $model = $clinseq_build->model;

    $self->status_message("\n***** " . $clinseq_build->__display_name__ . " ****");

    #Get build objects for all builds that might underly a clinseq model (Ref Align, Somatic variation, RNA-seq)
    my ($wgs_somvar_build, $exome_somvar_build, $tumor_rnaseq_build, $normal_rnaseq_build, $wgs_normal_refalign_build, $wgs_tumor_refalign_build, $exome_normal_refalign_build, $exome_tumor_refalign_build);

    $wgs_somvar_build = $clinseq_build->wgs_build;
    $exome_somvar_build = $clinseq_build->exome_build;
    $tumor_rnaseq_build = $clinseq_build->tumor_rnaseq_build;
    $normal_rnaseq_build = $clinseq_build->normal_rnaseq_build;
    $wgs_normal_refalign_build = $wgs_somvar_build->normal_build if ($wgs_somvar_build);
    $wgs_tumor_refalign_build = $wgs_somvar_build->tumor_build if ($wgs_somvar_build);
    $exome_normal_refalign_build = $exome_somvar_build->normal_build if ($exome_somvar_build);
    $exome_tumor_refalign_build = $exome_somvar_build->tumor_build if ($exome_somvar_build);

    #set the build types
    my %data_types;
    $self->_set_data_types($clinseq_build, \%data_types);

    #Gather all builds into a single array
    my @builds = ($wgs_normal_refalign_build, $wgs_tumor_refalign_build, $wgs_somvar_build, $exome_normal_refalign_build, $exome_tumor_refalign_build, $exome_somvar_build, $tumor_rnaseq_build, $normal_rnaseq_build, $clinseq_build);
    @builds = grep { $_ } @builds;

    #Get a list of sample names for samples associated with this clinseq build
    my %model_samples;
    my @model_samples;
    for my $build (@builds){
      next unless ($build->model->subject->class eq "Genome::Sample");
      $model_samples{$build->model->subject->id}{sample} = $build->model->subject;
    }
    foreach my $sample_id (keys %model_samples){
      push (@model_samples, $model_samples{$sample_id}{sample});
    }

    my $individual = $model->subject;
    $self->summarize_individual($individual);

    #Display instrument data counts for each sample/build actually associated with the clinseq model
    $self->status_message("\n\nSamples and instrument data counts (for samples associated with this clinseq model only - not all samples of the individual)");
    my %instdata_counts;
    for my $tn ('tumor', 'normal') { for my $type ('dna', 'rna') { $instdata_counts{$tn.$type} = 'n/a'; }}

    for my $sample (@model_samples) {
     my ($key, $value) = $self->summarize_sample($sample);
     $instdata_counts{$key} = $value;
    }

    print $stats_fh "Tumor Genomic DNA Instrument Data Count\t$instdata_counts{tumordna}\tlims\tClinseq Build Summary\tCount\tNumber of lanes of instrument data generated for tumor genomic DNA\n";
    print $stats_fh "Normal Genomic DNA Instrument Data Count\t$instdata_counts{normaldna}\tlims\tClinseq Build Summary\tCount\tNumber of lanes of instrument data generated for normal genomic DNA\n";
    print $stats_fh "Tumor RNA Instrument Data Count\t$instdata_counts{tumorrna}\tlims\tClinseq Build Summary\tCount\tNumber of lanes of instrument data generated for tumor RNA\n";
    print $stats_fh "Normal RNA Instrument Data Count\t$instdata_counts{normalrna}\tlims\tClinseq Build Summary\tCount\tNumber of lanes of instrument data generated for normal RNA\n";

    #Locations of useful methods need to do the following:
    #... /Genome/lib/perl/Genome/IntrumentData.pm
    #... /Genome/lib/perl/Genome/Model/ReferenceAlignment.pm
    #... /Genome/lib/perl/Genome/Model/Build/ReferenceAlignment.pm
    #... /Genome/lib/perl/Genome/Model/SomaticVariation.pm
    #... /Genome/lib/perl/Genome/Model/Build/SomaticVariation.pm
    #... /Genome/lib/perl/Genome/Model/RnaSeq.pm
    #... /Genome/lib/perl/Genome/Model/Build/RnaSeq.pm


    #Summarize the instrument data used by each model
    $self->status_message("\n\nInstrument data actually used by each build");
    for my $build (@builds){
      my $build_type = $build->type_name;
      my @instdata = $build->instrument_data;
      my $instdata_count = scalar(@instdata);

      if ($instdata_count > 0){
        $self->status_message("\nbuild " . $build->__display_name__ . " ($build_type)" . " uses " . $instdata_count . " instrument data");
        foreach my $instdata (@instdata){
          my $run_name = $instdata->run_name || "[UNDEF run_name]";
          $self->status_message("\t" . $instdata->id . "\t" . $run_name . "\t" . $instdata->library_name . "\t" . $instdata->sample_name);
        }
      }
    }

    #Summarize the build IDs and status of each build comprising the ClinSeq model
    $self->status_message("\n\nBuilds and status of each");
    for my $build (@builds) {
      $self->status_message("build '" . $build->__display_name__ . "' has status " . $build->status);
    }

    #List the build dirs for each build
    $self->status_message("\n\nBuild dir associated with each model/build");
    for my $build (@builds){
      my $build_type = $build->type_name;
      my $build_dir = $build->data_directory || "[UNDEF data_directory]";
      my $build_dn = $build->__display_name__;
      $self->status_message("$build_dir\t$build_type\t$build_dn");
    }

    #Summarize the processing profiles associated with each model
    $self->status_message("\n\nProcessing profiles associated with each model");
    for my $build (@builds){
      my $m = $build->model;
      my $pp = $m->processing_profile;
      my $pp_type = $pp->type_name;
      $self->status_message("model '" . $m->id . "' used processing profile '" . $pp->__display_name__ . "' ($pp_type)");
    }

    #Summarize the reference sequence build associated with each model
    $self->status_message("\n\nReference sequence build associated with each model");
    for my $build (@builds){
      my $m = $build->model;
      my $build_type = $build->type_name;
      if ($m->can("reference_sequence_build")){
        my $rb = $m->reference_sequence_build;
        $self->status_message("model '" . $m->__display_name__ . " ($build_type)" . "' used reference sequence build " . $rb->__display_name__);
      }
    }

    #Summarize the annotation reference build associated with each model
    $self->status_message("\n\nAnnotation reference build associated with each model");
    for my $build (@builds){
      my $m = $build->model;
      my $build_type = $build->type_name;
      if ($m->can("annotation_reference_build")){
        my $ab = $m->annotation_reference_build;
        if ($ab){
          my $ab_name = $ab->name || "[UNDEF annotation_name]";
          my $ab_dname = $ab->__display_name__ || "[UNDEF annotation_reference_build]";
          $self->status_message("model '" . $m->__display_name__ . " ($build_type)" . "' used annotation reference build " . $ab_dname . " ($ab_name)");
        }else{
          $self->status_message("model '" . $m->__display_name__ . " ($build_type)" . "' did NOT have an annotation reference build defined!");
        }
      }elsif ($m->can("annotation_build")){
        my $ab = $m->annotation_build;
        my $ab_name = $ab->name || "[UNDEF annotation_name]";
        my $ab_dname = $ab->__display_name__ || "[UNDEF annotation_reference_build]";
        $self->status_message("model '" . $m->__display_name__ . " ($build_type)" . "' used annotation reference build " . $ab_dname . " ($ab_name)");
      }
    }

    #Summarize the genotype microarray build used with each model
    $self->status_message("\n\nGenotype microarray build associated with each model");
    for my $build (@builds){
      my $m = $build->model;
      my $build_type = $build->type_name;
      if ($build->can("genotype_microarray_build")){
        my $gb = $build->genotype_microarray_build;
        if ($gb){
          my $gb_name = $gb->__display_name__ || "[UNDEF genotype_microarray_build]";
          $self->status_message("model '" . $m->__display_name__ . " ($build_type)" . "' used genotype microarray build " . $gb_name);
        }else{
          $self->status_message("model '" . $m->__display_name__ . " ($build_type)" . "' did NOT have an associated microarray build ");
        }
      }
    }

    #Summarize the dbsnp build used with each model
    $self->status_message("\n\ndbSNP build associated with each model");
    for my $build (@builds){
      my $m = $build->model;
      my $build_type = $build->type_name;
      if ($m->can("dbsnp_build")){
        my $db = $build->dbsnp_build;
        my $dm = $m->dbsnp_model;
        my $db_name = $db->__display_name__ || "[UNDEF dbsnp_build]";
        my $dm_name = $dm->__display_name__ || "[UNDEF dbsnp_model]";
        my $dm_id = $dm->id || "[UNDEF dbsnp_model]";

        $self->status_message("model '" . $m->__display_name__ . " ($build_type)" . "' used dbSNP build " . $db_name . " ($dm_id)");
      }
    }

    #Obtain the haploid coverage of each WGS/Exome ref alignment model
    #Obtain the array SNP concordance metrics for each ref alignment model (these are the 'Gold' SNP concordance values where SNPs called by WGS are compared to SNPs called by microarray of the same sample)
    #Tool that gives basic coverage for a ref align build ID
    #gmt analysis report-coverage --build-ids 120559668
    #Grab code from here: Genome/Model/Tools/Analysis/ReportCoverage.pm
    $self->status_message("\n\nHaploid coverage of each Exome/WGS ref alignment model");
    $self->status_message(join("\t",
        qw(subject_name pp_type sequence_type lane_count common_name
        tissue_desc extraction_type sequence_amount_gbp haploid_coverage array_het_snp_count
        array_het_snp_depth array_het_snp_concordance total_snp_positions_found_unfiltered
        total_snp_positions_found_filtered snp_positions_in_dbsnp snp_positions_not_in_dbsnp overall_dbsnp_concordance
        build_id)));
    for my $build (@builds) {
        $self->summarize_haploid_coverage_for_build($build, $stats_fh);
    }

    #Obtain the read duplication rate for each WGS BAM file from the Reference alignment models - note that duplication rate is defined on a per library basis (as it should be)
    #Also obtain other basic mapping information from the alignments files
    #Get per library info from the .metrics files in the alignments subdir of the refalign build dir
    #Get per BAM (all libraries combined) info from the .flagstat file in the alignments subdir of the refalign build dir
    $self->status_message("\n\nSample and library metrics from the reference alignment build directory");
    $self->status_message(join("\t",
        qw(subject_name pp_type sequence_type lane_count common_name tissue_desc extraction_type
        sequence_amount_gbp haploid_coverage sample_total_single_read_count sample_mapped_read_percent
        sample_properly_paired_read_percent sample_duplicate_read_percent)));
    for my $build (@builds) {
        $self->summarize_sample_and_library_metrics_for_build($build, $stats_fh);
    }

    $self->generate_APIPE_reports(\@builds, $build_outdir, $stats_fh);
    $self->generate_LIMS_reports(\@builds, $build_outdir);

    #Get BAMQC results for all ref-align builds
    for my $build (@builds) {
        next unless $self->_is_reference_alignment_build($build);
        $self->get_bamqc_results($build, $build_outdir, \%data_types);
        $self->get_perlane_bamqc_results($build, $build_outdir, \%data_types);
    }

    #Summarize exome coverage statistics for each WGS/Exome reference alignment model
    # cd /gscmnt/gc8001/info/model_data/2882774248/build120412367/reference_coverage/wingspan_0
    # cat *_STATS.tsv | perl -ne '@line=split("\t", $_); if ($line[12]==20){print "$_"}' | cut -f 8 | perl -ne 'chomp($_); $n++; $c+=$_; $a=$c/$n; print "$a\n”’
    my %exome_builds_with_coverage;
    $self->status_message("\n\nExome coverage values for each WGS/Exome reference alignment build");
    for my $build (@builds) {
      next unless $self->_is_reference_alignment_build($build);

      my $build_dir = $build->data_directory;
      my $tissue_desc = "[UNDEF tissue_desc]";
      my $extraction_type = "[UNDEF extraction_type]";
      my $subject = $build->subject;
      my $subject_name = $subject->name;
      my $common_name = $self->_get_subject_common_name($subject);
      if ($subject->can("tissue_desc")){
        $tissue_desc = $subject->tissue_desc;
      }
      if ($subject->can("extraction_type")){
        $extraction_type = $subject->extraction_type;
      }
      my @lanes = $build->instrument_data;
      my $sequence_type = $self->_determine_wgs_or_exome_for_instrument_data(@lanes);

      my $wingspan_0_dir = "$build_dir/reference_coverage/wingspan_0/";
      #my $wingspan_0_path = $wingspan_0_dir . "*_STATS.tsv";
      next unless (-e $wingspan_0_dir && -d $wingspan_0_dir);
      my $wingspan_0_file;

      opendir (my $dh, $wingspan_0_dir);
      my @files = readdir($dh);
      closedir($dh);
      foreach my $file (@files){
        if ($file =~ /\_STATS\.tsv/){
          $wingspan_0_file = $wingspan_0_dir . $file;
        }
      }

      next unless ($wingspan_0_file);

      my $build_id = $build->id;
      $exome_builds_with_coverage{$build_id}=1;

      $self->status_message("\nSample: $subject_name ($common_name | $tissue_desc | $extraction_type | $sequence_type)");

      #Summarize the average coverage of regions of interest (e.g. exons) at each minimum coverage cutoff value used when the reference coverage report was generated
      #Output file format(stats_file):
      #[0] Region Name (column 4 of BED file)
      #[1] Percent of Reference Bases Covered
      #[2] Total Number of Reference Bases
      #[3] Total Number of Covered Bases
      #[4] Number of Missing Bases
      #[5] Average Coverage Depth
      #[6] Standard Deviation Average Coverage Depth
      #[7] Median Coverage Depth
      #[8] Number of Gaps
      #[9] Average Gap Length
      #[10] Standard Deviation Average Gap Length
      #[11] Median Gap Length
      #[12] Min. Depth Filter
      #[13] Discarded Bases (Min. Depth Filter)
      #[14] Percent Discarded Bases (Min. Depth Filter)

      #Calculate:
      #total number of ROIs
      #average of all median ROI coverage levels
      #percent of all ROIs covered at X depth or greater across 80% or greater of breadth of ROIs
      my %covs;
      my $min_breadth = 80;
      open (REFCOV, "$wingspan_0_file");
      while(<REFCOV>){
        # cat *_STATS.tsv | perl -ne '@line=split("\t", $_); if ($line[12]==20){print "$_"}' | cut -f 8 | perl -ne 'chomp($_); $n++; $c+=$_; $a=$c/$n; print "$a\n”’
        chomp($_);
        my @line = split("\t", $_);
        my $min_cov = $line[12];
        my $percent_bases_covered = $line[1];
        my $median_coverage_depth = $line[7];
        if (defined($covs{$min_cov})){
          $covs{$min_cov}{count}++;
          $covs{$min_cov}{cum_median_coverage} += $median_coverage_depth;
          if ($percent_bases_covered > $min_breadth){
            $covs{$min_cov}{min_breadth_count}++;
          }
          $covs{$min_cov}{avg_median_coverage} = sprintf("%.2f", ($covs{$min_cov}{cum_median_coverage} / $covs{$min_cov}{count}));
          $covs{$min_cov}{min_breadth_count_percent} = sprintf("%.2f", (($covs{$min_cov}{min_breadth_count}/$covs{$min_cov}{count})*100));
        }else{
          $covs{$min_cov}{count} = 1;
          $covs{$min_cov}{cum_median_coverage} = $median_coverage_depth;
          $covs{$min_cov}{min_breadth_count} = 0;
          if ($percent_bases_covered > $min_breadth){
            $covs{$min_cov}{min_breadth_count}++;
          }
          $covs{$min_cov}{avg_median_coverage} = sprintf("%.2f", $median_coverage_depth);
          $covs{$min_cov}{min_breadth_count_percent} = sprintf("%.2f", (($covs{$min_cov}{min_breadth_count}/$covs{$min_cov}{count})*100));
        }
      }
      close (REFCOV);

      $self->status_message("min_coverage\troi_count\tmin_breadth_count_"."$min_breadth\tmin_breadth_count_percent_"."$min_breadth\taverage_median_coverage");
      foreach my $min_cov (sort {$a <=> $b} keys %covs){
        $self->status_message("$min_cov\t$covs{$min_cov}{count}\t$covs{$min_cov}{min_breadth_count}\t$covs{$min_cov}{min_breadth_count_percent}\t$covs{$min_cov}{avg_median_coverage}");
        print $stats_fh "Median ROI Coverage at >= $min_cov X\t$covs{$min_cov}{avg_median_coverage}\t$common_name\tClinseq Build Summary\tAverage\tAverage of Median Exon Coverage Values at a Min Coverage of $min_cov for $common_name $extraction_type data\n";
        print $stats_fh "Percent ROI Coverage at >= $min_cov X and >= $min_breadth % breadth\t$covs{$min_cov}{min_breadth_count_percent}\t$common_name\tClinseq Build Summary\tAverage\tPercent of Exons Covered at a Min Coverage of $min_cov and Min Breadth of $min_breadth for $common_name $extraction_type data\n";
      }
    }

    #Create ROI (i.e. exome) coverage view for a series of model IDs
    #   $ENV{GENOME_SYS_SERVICES_SEARCH_URL}/view/genome/model/set/coverage.html?&id=2880223426&id=2880225204&id=2880225668
    #Create coverage view for a series of build IDs
    #   $ENV{GENOME_SYS_SERVICES_SEARCH_URL}/view/genome/model/build/set/coverage.html?&id=115088238&id=115092892&id=115667927
    $self->status_message("\n\nCreate exome coverage reports for exome reference alignment builds");
    my $id_string = '';
    foreach my $exome_build (keys %exome_builds_with_coverage){
      $id_string .= "&id=$exome_build";
    }
    if ($id_string =~ /\d+/){
      my $cov_report_url = Genome::Utility::List::join_with_single_slash($ENV{GENOME_SYS_SERVICES_WEB_VIEW_URL}, "/genome/model/build/set/coverage.html?"."$id_string");
      $self->status_message("\nView here:\n$cov_report_url");
      $self->status_message("\n");
    }


    #Summarize basic RNA-seq alignment stats for each RNA-seq model
    #Files to grab or summarize if present
    #$build_dir/alignments/alignment_stats.txt (OLD)
    #$build_dir/expression/cufflinks.out (OLD)
    #$build_dir/metrics/PicardRnaSeqMetrics.txt (OLD)
    #$build_dir/metrics/PicardRnaSeqMetrics.png (OLD)
    #$build_dir/metrics/PicardRnaSeqChart.pdf (OLD)

    #$build_dir/alignment_stats/alignment_stats.txt (NEW)
    #$build_dir/junctions/summary/PercentGeneJunctionsCovered_BreadthvsDepth_BoxPlot.pdf (NEW)
    #$build_dir/junctions/summary/ObservedJunctions_SpliceSiteAnchorTypes_Pie.pdf
    #$build_dir/junctions/summary/ObservedJunctions_SpliceSiteUsage_Pie.pdf
    #$build_dir/junctions/summary/TranscriptJunctionReadCounts_Log2_Hist.pdf
    #$build_dir/bam-qc/*.pdf (NEW)
    #$build_dir/bam-qc/*.html (NEW)

    #$ENV{GENOME_SYS_SERVICES_FILES_URL}/$build_dir/bam-qc/
    #$ENV{GENOME_SYS_SERVICES_FILES_URL}/$build_dir/junctions/summary/

    
    for my $build (@builds) {
      next unless $self->_is_rna_seq_build($build);
      #Summarize RNA-seq metrics for each build
      $self->get_rnaseq_metrics($build, $build_outdir, $stats_fh);
    }

    #Summarize basic stats for the WGS and Exome somatic variation models - try using metrics object of somatic variation build?
    #e.g. number of tier 1,2,3,4 SNVs and InDels (both 'novel' and 'previously known')
    #e.g. number of SVs
    $self->status_message("\n\nGet basic somatic variation stats");
    $self->status_message("pp_name\ttier1_snv_count\ttier2_snv_count\ttier3_snv_count\ttier4_snv_count\ttier1_indel_count\ttier2_indel_count\ttier3_indel_count\ttier4_indel_count\tsv_count\tbuild_id");
    for my $build (@builds) {
      my $m = $build->model;
      my $pp = $m->processing_profile;
      my $pp_name = $pp->name;
      my $data_type = $self->_determine_wgs_or_exome_for_build($build, \%data_types);

      next unless $self->_is_somatic_variation_build($build);

      my $build_id = $build->id;
      my $tissue_desc = "[UNDEF tissue_desc]";
      my $extraction_type = "[UNDEF extraction_type]";
      my $subject = $build->subject;
      my $subject_name = $subject->name;
      my $common_name = $self->_get_subject_common_name($subject);
      if ($subject->can("tissue_desc")){
        if ($subject->tissue_desc){
          $tissue_desc = $subject->tissue_desc;
        }
      }
      if ($subject->can("extraction_type")){
        if ($subject->extraction_type){
          $extraction_type = $subject->extraction_type;
        }
      }
      #grab metrics and build a hash from them
      my %metrics = map {$_->name => $_->value} $build->metrics;
      my $tier1_snv_count = "n/a";
      my $tier2_snv_count = "n/a";
      my $tier3_snv_count = "n/a";
      my $tier4_snv_count = "n/a";
      my $tier1_indel_count = "n/a";
      my $tier2_indel_count = "n/a";
      my $tier3_indel_count = "n/a";
      my $tier4_indel_count = "n/a";
      my $sv_count = "n/a";

      if (defined($metrics{'tier1_snv_count'})){
        $tier1_snv_count = $metrics{'tier1_snv_count'};
      }
      if (defined($metrics{'tier2_snv_count'})){
        $tier2_snv_count = $metrics{'tier2_snv_count'};
      }
      if (defined($metrics{'tier3_snv_count'})){
        $tier3_snv_count = $metrics{'tier3_snv_count'};
      }
      if (defined($metrics{'tier4_snv_count'})){
        $tier4_snv_count = $metrics{'tier4_snv_count'};
      }
      if (defined($metrics{'tier1_indel_count'})){
        $tier1_indel_count = $metrics{'tier1_indel_count'};
      }
      if (defined($metrics{'tier2_indel_count'})){
        $tier2_indel_count = $metrics{'tier2_indel_count'};
      }
      if (defined($metrics{'tier3_indel_count'})){
        $tier3_indel_count = $metrics{'tier3_indel_count'};
      }
      if (defined($metrics{'tier4_indel_count'})){
        $tier4_indel_count = $metrics{'tier4_indel_count'};
      }
      if (defined($metrics{'sv_count'})){
        $sv_count = $metrics{'sv_count'};
      }

      $self->status_message("$pp_name\t$tier1_snv_count\t$tier2_snv_count\t$tier3_snv_count\t$tier4_snv_count\t$tier1_indel_count\t$tier2_indel_count\t$tier3_indel_count\t$tier4_indel_count\t$sv_count\t$build_id");

      print $stats_fh "SomVar Tier1 SNV Count\t$tier1_snv_count\t$data_type\tClinseq Build Summary\tCount\tSomatic variation tier 1 SNV count for $common_name $extraction_type data\n";
      print $stats_fh "SomVar Tier2 SNV Count\t$tier2_snv_count\t$data_type\tClinseq Build Summary\tCount\tSomatic variation tier 2 SNV count for $common_name $extraction_type data\n";
      print $stats_fh "SomVar Tier3 SNV Count\t$tier3_snv_count\t$data_type\tClinseq Build Summary\tCount\tSomatic variation tier 3 SNV count for $common_name $extraction_type data\n";
      print $stats_fh "SomVar Tier4 SNV Count\t$tier4_snv_count\t$data_type\tClinseq Build Summary\tCount\tSomatic variation tier 4 SNV count for $common_name $extraction_type data\n";
      print $stats_fh "SomVar Tier1 INDEL Count\t$tier1_indel_count\t$data_type\tClinseq Build Summary\tCount\tSomatic variation tier 1 INDEL count for $common_name $extraction_type data\n";
      print $stats_fh "SomVar Tier2 INDEL Count\t$tier2_indel_count\t$data_type\tClinseq Build Summary\tCount\tSomatic variation tier 2 INDEL count for $common_name $extraction_type data\n";
      print $stats_fh "SomVar Tier3 INDEL Count\t$tier3_indel_count\t$data_type\tClinseq Build Summary\tCount\tSomatic variation tier 3 INDEL count for $common_name $extraction_type data\n";
      print $stats_fh "SomVar Tier4 INDEL Count\t$tier4_indel_count\t$data_type\tClinseq Build Summary\tCount\tSomatic variation tier 4 INDEL count for $common_name $extraction_type data\n";
      print $stats_fh "SomVar SV Count\t$sv_count\t$data_type\tClinseq Build Summary\tCount\tSomatic variation SV count for $common_name $extraction_type data\n";
    }

    #Summarize SV annotation file from somatic variation results
    $self->status_message("\n\nGet more detailed merged somatic SV stats");
    $self->status_message("pp_name\tsv_count\tctx_count\tdel_count\tinv_count\titx_count\tbuild_id");
    for my $build (@builds) {
      my $m = $build->model;
      my $pp = $m->processing_profile;
      my $pp_name = $pp->name;
      my $data_type = $self->_determine_wgs_or_exome_for_build($build, \%data_types);
      next unless ($self->_is_somatic_variation_build($build) && $data_type eq "WGS");

      my $build_id = $build->id;
      my $build_dir = $build->data_directory;


      my $sv_annot_search1 = $build_dir . "/variants/sv/union-union*/svs.hq.merge.annot.somatic";
      my $sv_annot_search2 = $build_dir . "/variants/sv/union-sv*/svs.hq.merge.annot.somatic";

      my $sv_annot_file = 'NULL';

      my $sv_annot_file1 = `ls $sv_annot_search1 2>/dev/null` || "NULL";
      chomp($sv_annot_file1);

      my $sv_annot_file2 = `ls $sv_annot_search2 2>/dev/null` || "NULL";
      chomp($sv_annot_file2);

      if (-e $sv_annot_file1){
        $sv_annot_file = $sv_annot_file1;
      }elsif (-e $sv_annot_file2){
        $sv_annot_file = $sv_annot_file2;
      }

      chomp($sv_annot_file);
      my $sv_count = 0;
      my $ctx_count = 0;
      my $del_count = 0;
      my $inv_count = 0;
      my $itx_count = 0;
      if (-e $sv_annot_file){
        open (SV, "$sv_annot_file") || die "\n\nCould not open SV annot file: $sv_annot_file\n\n";
        while(<SV>){
          chomp($_);
          next if ($_ =~ /^\#/);
          $sv_count++;
          my @line = split("\t", $_);
          if ($_ =~ /\s+CTX\s+/){$ctx_count++;}
          if ($_ =~ /\s+DEL\s+/){$del_count++;}
          if ($_ =~ /\s+INV\s+/){$inv_count++;}
          if ($_ =~ /\s+ITX\s+/){$itx_count++;}
        }
        close(SV);
        $self->status_message("$pp_name\t$sv_count\t$ctx_count\t$del_count\t$inv_count\t$itx_count\t$build_id");
      }else{
        $self->status_message("Could not find SV file with search:\nls $sv_annot_search1\nls $sv_annot_search2");
      }
    }

    #Summarize SV annotation file from somatic variation results
    $self->status_message("\n\nGet basic expression count from RNA-seq");
    $self->status_message("pp_name\tgenes_fpkm_greater_1\ttranscripts_fpkm_greater_1\tbuild_id");
    for my $build (@builds) {
      my $m = $build->model;
      my $pp = $m->processing_profile;
      my $pp_name = $pp->name;
      my $subject = $build->subject;
      my $subject_name = $subject->name;
      my $common_name = $self->_get_subject_common_name($subject);

      next unless $self->_is_rna_seq_build($build);

      my $build_id = $build->id;
      my $build_dir = $build->data_directory;
      my $genes_fpkm_greater_1 = 0;
      my $transcripts_fpkm_greater_1 = 0;

      my $ge_file = $build_dir . "/expression/genes.fpkm_tracking";
      if (-e $ge_file){
        open (GE, "$ge_file") || die "\n\nCould not open gene expression file: $ge_file\n\n";
        while(<GE>){
          chomp($_);
          if ($_ =~ /tracking/){
            next();
          }
          my @line = split("\t", $_);
          if ($line[9] > 1){
            $genes_fpkm_greater_1++;
          }
        }
        close(GE);
      }else{
        $self->status_message("Could not find gene expression file in RNA-seq build dir: $ge_file");
        $genes_fpkm_greater_1 = "n/a";
      }

      my $te_file = $build_dir . "/expression/isoforms.fpkm_tracking";
      if (-e $te_file){
        open (TE, "$te_file") || die "\n\nCould not open transcript expression file: $te_file\n\n";
        while(<TE>){
          chomp($_);
          if ($_ =~ /tracking/){
            next();
          }
          my @line = split("\t", $_);
          if ($line[9] > 1){
            $transcripts_fpkm_greater_1++;
          }
        }
        close(TE);
      }else{
        $self->status_message("Could not find transcript expression file in RNA-seq build dir: $te_file");
        $transcripts_fpkm_greater_1 = "n/a";
      }
      $self->status_message("$pp_name\t$genes_fpkm_greater_1\t$transcripts_fpkm_greater_1\t$build_id");
    }


    #Print BAMs for all reference alignment and RNA-seq builds
    $self->status_message("\n\nGet all BAM file locations");
    for my $build (@builds) {
      next unless($self->_is_reference_alignment_build($build) or $self->_is_rna_seq_build($build));

      my $build_id = $build->id;
      my $build_dir = $build->data_directory;
      my $tissue_desc = "[UNDEF tissue_desc]";
      my $extraction_type = "[UNDEF extraction_type]";
      my $subject = $build->subject;
      my $subject_name = $subject->name;
      my $common_name = $self->_get_subject_common_name($subject);
      if ($subject->can("tissue_desc")){
        if ($subject->tissue_desc){
          $tissue_desc = $subject->tissue_desc;
        }
      }
      if ($subject->can("extraction_type")){
        if ($subject->extraction_type){
          $extraction_type = $subject->extraction_type;
        }
      }

      my $bam_file = "NOT FOUND";
      if ($build->can("whole_rmdup_bam_file")){
        $bam_file = $build->whole_rmdup_bam_file;
      }elsif($build->can("alignment_result")){
        my $alignment_result = $build->alignment_result;
        if ($alignment_result->can("bam_file")){
          $bam_file = $alignment_result->bam_file;
        }
      }
      $self->status_message("$subject_name ($common_name | $tissue_desc | $extraction_type)\t$bam_file");
    }
    close($stats_fh);

    return 1;
}

sub summarize_individual {
    my $self = shift;
    my $individual = shift;

    #Summarize Individual (id, name, gender, etc.). Check if all models hit the same individual and warn if not
    #... /Genome/lib/perl/Genome/Individual.pm
    $self->status_message("\n\nIndividual (subject/patient) information");
    my $individual_id = $individual->id || "[UNDEF individual id]";
    my $individual_type = $individual->subject_type || "[UNDEF individual subject type]";
    my $individual_species_name = $individual->species_name || "[UNDEF individual species name]";
    my $individual_name = $individual->name || "[UNDEF individual name]";
    my $individual_gender = $individual->gender || "[UNDEF individual gender]";
    my $individual_common_name = $individual->common_name || "[UNDEF individual common name]";
    my $individual_upn =  $individual->upn || "[UNDEF individual upn]";
    my $individual_ethnicity =  $individual->ethnicity || "[UNDEF individual ethnicity]";
    my $individual_race =  $individual->race || "[UNDEF individual race]";
    $self->status_message("individual_id: $individual_id");
    $self->status_message("individual_type: $individual_type");
    $self->status_message("individual_species_name: $individual_species_name");
    $self->status_message("individual_name: $individual_name");
    $self->status_message("individual_gender: $individual_gender");
    $self->status_message("individual_common_name: $individual_common_name");
    $self->status_message("individual_upn: $individual_upn");
    $self->status_message("individual_ethnicity: $individual_ethnicity");
    $self->status_message("individual_race: $individual_race");

    return 1;
}

sub summarize_sample {
    my $self = shift;
    my $sample = shift;

    my $sample_id = $sample->id;
    my @instdata = $sample->instrument_data;
    my $scn = $sample->common_name || "[UNDEF sample common_name]";
    my $tissue_desc = $sample->tissue_desc || "[UNDEF sample tissue_desc]";
    my $extraction_type = $sample->extraction_type || "[UNDEF sample extraction_type]";
    $self->status_message("sample " . $sample->__display_name__ . " ($tissue_desc - $extraction_type) has " . scalar(@instdata) . " instrument data");
    my $normal_sample_common_names = $self->normal_sample_common_names;
    my $tumor_sample_common_names = $self->tumor_sample_common_names;

    my ($tumor_or_normal, $dna_or_rna) = ('n/a') x 2;
    if ($scn =~ /$tumor_sample_common_names/i) {
       $tumor_or_normal = 'tumor';
    } elsif ($scn =~ /$normal_sample_common_names/i) {
        $tumor_or_normal = 'normal';
    }

    if($extraction_type eq 'genomic dna') {
        $dna_or_rna = 'dna';
    } elsif ($extraction_type eq 'rna') {
        $dna_or_rna = 'rna';
    }

    my $count = scalar(@instdata);

    return ($tumor_or_normal . $dna_or_rna) => $count;
}

sub summarize_haploid_coverage_for_build {
    my $self = shift;
    my $build = shift;
    my $stats_fh = shift;

    #Only perform the following for reference alignment builds!
    return unless $self->_is_reference_alignment_build($build);

    my $build_type = $build->type_name;
    my $build_dir = $build->data_directory;
    my $tissue_desc = "[UNDEF tissue_desc]";
    my $extraction_type = "[UNDEF extraction_type]";
    my $subject = $build->subject;
    my $subject_name = $subject->name;
    my $common_name = $self->_get_subject_common_name($subject);
    if ($subject->can("tissue_desc")){
        if ($subject->tissue_desc){
            $tissue_desc = $subject->tissue_desc;
        }
    }
    if ($subject->can("extraction_type")){
        if ($subject->extraction_type){
            $extraction_type = $subject->extraction_type;
        }
    }
    #grab metrics and build a hash from them
    my %metrics = map {$_->name => $_->value} $build->metrics;
    my $build_id = $build->id;
    my $gbp = "n/a";
    my $haploid_coverage = "n/a";

    #Get total amount of instrument data from the 'metrics' object
    if (defined($metrics{"instrument data total kb"})){
        $gbp = sprintf("%0.01f", $metrics{"instrument data total kb"} / 10**6);
    }

    #Get haploid coverage from the 'metrics' object
    if (defined($metrics{'haploid_coverage'})){
        $haploid_coverage = sprintf("%0.01f", $metrics{'haploid_coverage'});
    }

    my @lanes = $build->instrument_data;
    my $sequence_type = $self->_determine_wgs_or_exome_for_instrument_data(@lanes);
    my $lane_count = scalar(@lanes);

    #Get array vs. sequence concordance of SNP positions (only those on the Illumina bead array).  These are referred to in the system as 'Gold SNP Concordance' values
    my $gold_filtered_het_snp_count = "n/a";                 # 'filtered:heterozygous snv:hits:match:heterozygous (1 alleles) snv:count'
    my $gold_filtered_total = "n/a";                         # 'filtered:heterozygous snv:total'
    my $gold_filtered_het_snp_depth = "n/a";                 # 'filtered:heterozygous snv:hits:match:heterozygous (1 alleles) snv:qual'
    my $gold_filtered_het_snp_percent_concordance = "n/a";   # ($gold_filtered_het_snp_count / $gold_filtered_total) * 100

    #Note, do not calculate for Exome as the SNP concordance values are not being calculated correctly... not accounting for the target space when determining percent concordance!
    if (($sequence_type eq "wgs" || $sequence_type eq "mixed") && defined($metrics{'filtered:heterozygous snv:hits:match:heterozygous (1 alleles) snv:count'}) && defined($metrics{'filtered:heterozygous snv:total'}) && defined($metrics{'filtered:heterozygous snv:hits:match:heterozygous (1 alleles) snv:qual'})){
        $gold_filtered_het_snp_count = $metrics{'filtered:heterozygous snv:hits:match:heterozygous (1 alleles) snv:count'};
        $gold_filtered_total = $metrics{'filtered:heterozygous snv:total'};
        $gold_filtered_het_snp_depth = $metrics{'filtered:heterozygous snv:hits:match:heterozygous (1 alleles) snv:qual'};
        if ($gold_filtered_total){
            $gold_filtered_het_snp_percent_concordance = sprintf("%.2f", (($gold_filtered_het_snp_count / $gold_filtered_total)*100));
        }
    }

    #Obtain the dbSNP concordance metrics for each ref alignment model (these are the SNPs called by WGS compared to publicly known dbSNPS)
    #NOTE: These metrics do not appear to be stored in a metrics object...  Will have to parse from a file routinely created for reference alignment builds
    #e.g. /gscmnt/gc7001/info/model_data/2879615516/build114445127/reports/dbsnp_concordance.filtered.txt
    #"Overall concordance: 93.6585"  "Total Snvs: 3656996"  "Hits: 3425089"  "Misses: 231907"
    my $total_snp_positions_found_filtered = "n/a";
    my $snp_positions_in_dbsnp = "n/a";
    my $snp_positions_not_in_dbsnp = "n/a";
    my $overall_dbsnp_concordance = "n/a";
    my $dbsnp_concordance_file = $build_dir . "/reports/dbsnp_concordance.filtered.txt";
    if (-e $dbsnp_concordance_file){
        open (DBSNPC, "$dbsnp_concordance_file");
        while(<DBSNPC>){
            chomp($_);
            if ($_ =~ /Overall\s+concordance\:\s+(\S+)/){
                $overall_dbsnp_concordance = $1;
                if ($overall_dbsnp_concordance =~ /\d+/){
                    $overall_dbsnp_concordance = sprintf("%.2f", $overall_dbsnp_concordance)
                }
            }
            if ($_ =~ /Hits\:\s+(\S+)/){
                $snp_positions_in_dbsnp = $1;
            }
            if ($_ =~ /Misses\:\s+(\S+)/){
                $snp_positions_not_in_dbsnp = $1;
            }
            if ($_ =~ /Total\s+Snvs\:\s+(\S+)/){
                $total_snp_positions_found_filtered = $1;
            }
        }
        close (DBSNPC);
    }

    #Get the total unfiltered SNP calls for each reference alignment build
    my $total_snp_positions_found_unfiltered = "n/a";
    $dbsnp_concordance_file = $build_dir . "/reports/dbsnp_concordance.txt";
    if (-e $dbsnp_concordance_file){
        open (DBSNPC, "$dbsnp_concordance_file");
        while(<DBSNPC>){
            chomp($_);
            if ($_ =~ /Total\s+Snvs\:\s+(\S+)/){
                $total_snp_positions_found_unfiltered = $1;
            }
        }
        close (DBSNPC);
    }

    #Display gathered info for both Exome and WGS reference alignment models
    $self->status_message(join("\t",
        $subject_name, $build_type, $sequence_type, $lane_count, $common_name,
        $tissue_desc, $extraction_type, $gbp, $haploid_coverage, $gold_filtered_het_snp_count,
        $gold_filtered_het_snp_depth, $gold_filtered_het_snp_percent_concordance,
        $total_snp_positions_found_unfiltered, $total_snp_positions_found_filtered,
        $snp_positions_in_dbsnp, $snp_positions_not_in_dbsnp, $overall_dbsnp_concordance,
        $build_id
    ));

    #Resolve data type
    my $data_type;
    if ($subject->can("common_name")){
        if ($subject->common_name){
            $data_type = $sequence_type . "_" . $common_name;
        }else{
            $data_type = $sequence_type . "_" . $subject_name;
        }
    }else{
        $data_type = $sequence_type . "_" . $subject_name;
    }

    print $stats_fh "Data amount (Gbp)\t$gbp\t$data_type\tClinseq Build Summary\tFloat\tData amount (Gbp) for $sequence_type $common_name data\n";
    print $stats_fh "Haploid coverage\t$haploid_coverage\t$data_type\tClinseq Build Summary\tFloat\tHaploid coverage for $sequence_type $common_name data\n";
    print $stats_fh "Gold SNP Concordance\t$gold_filtered_het_snp_percent_concordance\t$data_type\tClinseq Build Summary\tPercent\tSNP array vs. sequencing SNP concordance (gold, filtered, het, snps) for $sequence_type $common_name data\n";
    print $stats_fh "dbSNP SNP Concordance\t$overall_dbsnp_concordance\t$data_type\tClinseq Build Summary\tPercent\tdbSNP vs. sequencing SNP concordance for $sequence_type $common_name data\n";

    return 1;
}

sub summarize_sample_and_library_metrics_for_build {
    my $self = shift;
    my $build = shift;
    my $stats_fh = shift;

    return unless $self->_is_reference_alignment_build($build);

    my $build_type = $build->type_name;
    my $build_dir = $build->data_directory;
    my $tissue_desc = "[UNDEF tissue_desc]";
    my $extraction_type = "[UNDEF extraction_type]";
    my $subject = $build->subject;
    my $subject_name = $subject->name;
    my $common_name = $self->_get_subject_common_name($subject);
    if ($subject->can("tissue_desc")){
        if ($subject->tissue_desc){
            $tissue_desc = $subject->tissue_desc;
        }
    }
    if ($subject->can("extraction_type")){
        if ($subject->extraction_type){
            $extraction_type = $subject->extraction_type;
        }
    }
    #grab metrics and build a hash from them
    my %metrics = map {$_->name => $_->value} $build->metrics;
    my $build_id = $build->id;
    my $gbp = "n/a";
    my $haploid_coverage = "n/a";

    #Get total amount of instrument data from the 'metrics' object
    if (defined($metrics{"instrument data total kb"})){
        $gbp = sprintf("%0.01f", $metrics{"instrument data total kb"} / 10**6);
    }

    #Get haploid coverage from the 'metrics' object
    if (defined($metrics{'haploid_coverage'})){
        $haploid_coverage = sprintf("%0.01f", $metrics{'haploid_coverage'});
    }

    my @lanes = $build->instrument_data;
    my $sequence_type = $self->_determine_wgs_or_exome_for_instrument_data(@lanes);
    my $lane_count = scalar(@lanes);

    #Parse the flagstat file for sample metrics of the BAM file
    my $sample_total_single_read_count = "n/a";
    my $sample_mapped_read_percent = "n/a";
    my $sample_properly_paired_read_percent = "n/a";
    my $sample_duplicate_read_percent = "n/a";

    my $flagstat_file = $build->merged_alignment_result->merged_alignment_bam_flagstat;
    if (-e $flagstat_file){
        my $flagstat_data = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_file);
        $sample_total_single_read_count = $flagstat_data->{total_reads};

        $sample_mapped_read_percent = $flagstat_data->{reads_mapped_percentage}
            if $flagstat_data->{reads_mapped_percentage};

        $sample_properly_paired_read_percent = $flagstat_data->{reads_mapped_in_proper_pairs_percentage}
            if $flagstat_data->{reads_mapped_in_proper_pairs_percentage};

        my $duplicate_count = $flagstat_data->{reads_marked_duplicates};
        $sample_duplicate_read_percent = sprintf("%.2f", (($duplicate_count/$sample_total_single_read_count)*100));
    }else{
        $self->status_message("Warning: Could not find flagstat file: $flagstat_file");
    }

    $self->status_message(join("\t",
            $subject_name, $build_type, $sequence_type, $lane_count, $common_name, $tissue_desc, $extraction_type,
            $gbp, $haploid_coverage, $sample_total_single_read_count, $sample_mapped_read_percent,
            $sample_properly_paired_read_percent, $sample_duplicate_read_percent
        ));

    #Resolve data type
    my $data_type;
    if ($subject->can("common_name")){
        if ($subject->common_name){
            $data_type = $sequence_type . "_" . $common_name;
        }else{
            $data_type = $sequence_type . "_" . $subject_name;
        }
    }else{
        $data_type = $sequence_type . "_" . $subject_name;
    }

    print $stats_fh "Total Single Read Count\t$sample_total_single_read_count\t$data_type\tClinseq Build Summary\tCount\tTotal single read count for $sequence_type $common_name data\n";
    print $stats_fh "Percent Reads Mapped\t$sample_mapped_read_percent\t$data_type\tClinseq Build Summary\tPercent\tRead mapping percent for $sequence_type $common_name data\n";
    print $stats_fh "Percent Reads Properly Paired\t$sample_properly_paired_read_percent\t$data_type\tClinseq Build Summary\tPercent\tPercent of reads that are properly paired for $sequence_type $common_name data\n";
    print $stats_fh "Read Duplication Rate (sample level)\t$sample_duplicate_read_percent\t$data_type\tClinseq Build Summary\tPercent\tPercent read duplication at sample level (all libraries combined) for $sequence_type $common_name data\n";

    my $alignments_dir = $build_dir . "/alignments/";

    #Search for metrics files.  Hard to predict which build ID will be used for the name of this file because of short cutting
    my @metrics_files;
    opendir (my $dh, $alignments_dir);
    my @files = readdir($dh);
    closedir($dh);
    foreach my $file (@files){
        if ($file =~ /\.metrics/){
            my $metrics_file = $alignments_dir . $file;
            push(@metrics_files, $metrics_file);
        }
    }

    $self->status_message("\tlibrary_name\tlibrary_duplication_rate");
    foreach my $file (@metrics_files){
        my $library_name = "n/a";
        my $library_duplication_rate = "n/a";
        open (METRICS, "$file") || die "\n\nCould not open file: $file\n\n";
        while(<METRICS>){
            chomp($_);
            my @line = split("\t", $_);
            my $col_count = scalar(@line);
            if ($col_count == 9 && $line[7] =~ /\d+/){
                $library_name = $line[0];
                $library_duplication_rate = sprintf("%.2f", ($line[7]*100));
                $self->status_message("\t$library_name\t$library_duplication_rate");
            }
        }
        close (METRICS);
    }
    return 1;
}


sub get_rnaseq_metrics {
    my $self = shift;
    my $rnaseq_build = shift;
    my $build_outdir = shift;
    my $stats_fh = shift;
    my $rnaseq_build_id = $rnaseq_build->id;
    my $rnaseq_build_dir = $rnaseq_build->data_directory;
    my $tissue_desc = "[UNDEF tissue_desc]";
    my $extraction_type = "[UNDEF extraction_type]";
    my $subject = $rnaseq_build->subject;
    my $subject_name = $subject->name;
    my $common_name = $self->_get_subject_common_name($subject);
    if ($subject->can("tissue_desc")){
      $tissue_desc = $subject->tissue_desc;
    }
    if ($subject->can("extraction_type")){
      $extraction_type = $subject->extraction_type;
    }
    $self->status_message("\n\nGet basic RNA-seq alignment stats");
    $self->status_message("\nsample");
    $self->status_message("$subject_name ($common_name | $tissue_desc | $extraction_type)");
    $self->get_rnaseq_alignment_stats($common_name, $extraction_type, $rnaseq_build, 
          $rnaseq_build_dir, $rnaseq_build_id, $stats_fh);
    $self->get_cufflinks_metrics($common_name, $extraction_type, $rnaseq_build_dir,
          $rnaseq_build_id, $stats_fh);
    $self->get_picard_metrics($common_name, $extraction_type, $rnaseq_build_dir, 
          $rnaseq_build_id, $stats_fh);
    $self->copy_from_rnaseq_build($rnaseq_build_dir, $build_outdir, $common_name);
    $self->display_handy_rnaseq_urls($rnaseq_build_dir);
}

sub display_handy_rnaseq_urls {
    my $self = shift;
    my $build_dir = shift;
    #display urls to some handy locations in the rna-seq build
    if (-e "$build_dir/bam-qc/"){
      $self->status_message("\nRNA-seq BAM-QC results:\n"
        . Genome::Utility::List::join_with_single_slash($ENV{GENOME_SYS_SERVICES_FILES_URL},
          "$build_dir/bam-qc/"));
    }
    if (-e "$build_dir/junctions/summary/"){
      $self->status_message("\nRNA-seq junction summary results:\n"
        . Genome::Utility::List::join_with_single_slash($ENV{GENOME_SYS_SERVICES_FILES_URL},
          "$build_dir/junctions/summary/"));
    }
}

sub get_rnaseq_alignment_stats {
    my $self = shift;
    my $common_name = shift;
    my $extraction_type = shift;
    my $rnaseq_build = shift;
    my $build_dir = shift;
    my $build_id = shift;
    my $stats_fh = shift;
    
    my $alignment_stats_file = $rnaseq_build->alignment_stats_file();
    my $total_top_alignments = 0;
    my $total_top_spliced_alignments = 0;
    my $mt_top_alignments = 0;

    my $total_reads = "n/a";
    my $unmapped_reads_p = "n/a";
    my $total_reads_mapped_p = "n/a";
    if (-e $alignment_stats_file){
      my $unmapped_reads = "n/a";
      my $total_reads_mapped = "n/a";
      open (ALIGN, "$alignment_stats_file");
      while(<ALIGN>){
        chomp($_);
        my @line = split("\t", $_);
        if ($_ =~ /Total\s+Reads\:\s+(\d+)/){
          $total_reads = $1;
        }
        if ($_ =~ /Unmapped\s+Reads\:\s+(\d+)/){
          $unmapped_reads = $1;
        }
        if ($_ =~ /Total\s+Reads\s+Mapped\:\s+(\d+)/){
          $total_reads_mapped = $1;
        }
        if ($total_reads =~ /\d+/ && $unmapped_reads =~ /\d+/ && $total_reads_mapped =~ /\d+/){
          if ($total_reads > 0){
            $unmapped_reads_p = sprintf("%.2f", ($unmapped_reads/$total_reads)*100);
            $total_reads_mapped_p = sprintf("%.2f", ($total_reads_mapped/$total_reads)*100);
          }
        }
        if (scalar(@line) == 9){
          next if ($line[0] =~ /^chr$/);
          $total_top_alignments += $line[1];
          $total_top_spliced_alignments += $line[2];
          if ($line[0] eq "MT"){
            $mt_top_alignments = $line[1];
          }
        }
      }
      close (ALIGN);
    }else{
      $self->status_message("Could not find alignment_stats.txt file for build: $build_id");
    }
    my $spliced_alignments_p = "n/a";
    my $mt_alignments_p = "n/a";
    if ($total_top_alignments){
      $spliced_alignments_p = sprintf("%.2f", (($total_top_spliced_alignments/$total_top_alignments)*100));
      $mt_alignments_p = sprintf("%.2f", (($mt_top_alignments/$total_top_alignments)*100));
    }
    print $stats_fh "RNA-seq Total Reads\t$total_reads\t$common_name\tClinseq Build Summary\tCount\tTotal RNA-seq reads for $common_name $extraction_type data\n";
    print $stats_fh "RNA-seq Percent Spliced Alignments\t$spliced_alignments_p\t$common_name\tClinseq Build Summary\tPercent\tPercent of RNA-seq reads mapped across splice junctions for $common_name $extraction_type data\n";
    print $stats_fh "RNA-seq Percent MT Alignments\t$mt_alignments_p\t$common_name\tClinseq Build Summary\tPercent\tPercent of RNA-seq reads mapped to the MT chromosome for $common_name $extraction_type data\n";
    print $stats_fh "RNA-seq Percent Reads Mapped\t$total_reads_mapped_p\t$common_name\tClinseq Build Summary\tPercent\tPercent of RNA-seq reads mapped for $common_name $extraction_type data\n";
    print $stats_fh "RNA-seq Percent Reads UnMapped\t$unmapped_reads_p\t$common_name\tClinseq Build Summary\tPercent\tPercent of RNA-seq reads unmapped for $common_name $extraction_type data\n";
    $self->status_message("total_reads\ttotal_reads_mapped_percent\tunmapped_reads_percent");
    $self->status_message("$total_reads\t$total_reads_mapped_p\t$unmapped_reads_p");
}

sub get_cufflinks_metrics {
    my $self = shift;
    my $common_name = shift;
    my $extraction_type = shift;
    my $build_dir = shift;
    my $build_id = shift;
    my $stats_fh = shift;
    
    #/gscmnt/gc2016/info/model_data/2880794613/build115909698/expression/cufflinks.out
    my $cufflinks_out_file = $build_dir . "/expression/cufflinks.out";
    my $frag_size_mean = "n/a";
    my $frag_size_std = "n/a";
    if (-e $cufflinks_out_file){
      open (CUFF, "$cufflinks_out_file");
      while(<CUFF>){
        chomp($_);
        if ($_ =~ /Estimated\s+Mean\:\s+(\S+)/){
          $frag_size_mean = $1;
        }
        if ($_ =~ /Estimated\s+Std\s+Dev\:\s+(\S+)/){
          $frag_size_std = $1;
        }
      }
      close(CUFF);
    }else{
      $self->status_message("Could not find cufflinks.out file for build: $build_id");
    }
    print $stats_fh "RNA-seq Mean Fragment Size\t$frag_size_mean\t$common_name\tClinseq Build Summary\tFloat\tMean cDNA fragment size inferred from RNA-seq reads for $common_name $extraction_type data\n";
    print $stats_fh "RNA-seq StDev of Fragment Size\t$frag_size_std\t$common_name\tClinseq Build Summary\tFloat\tStandard deviation of cDNA fragment size inferred from RNA-seq reads for $common_name $extraction_type data\n";
    $self->status_message("fragment_size_mean\tfragment_size_std");
    $self->status_message("$frag_size_mean\t$frag_size_std");
}

sub get_picard_metrics {
    my $self = shift;
    my $common_name = shift;
    my $extraction_type = shift;
    my $build_dir = shift;
    my $build_id = shift;
    my $stats_fh = shift;
    #/gscmnt/gc2016/info/model_data/2880794613/build115909698/metrics/PicardRnaSeqMetrics.txt
    my $picard_metrics_file = $build_dir . "/metrics/PicardRnaSeqMetrics.txt";
    my $pct_ribosomal_bases = "n/a";
    my $pct_coding_bases = "n/a";
    my $pct_utr_bases = "n/a";
    my $pct_intronic_bases = "n/a";
    my $pct_intergenic_bases = "n/a";
    if (-e $picard_metrics_file){
      open (PIC, "$picard_metrics_file");
      while(<PIC>){
        chomp($_);
        next if ($_ =~ /^\#/);
        my @line = split("\t", $_);
        if (scalar(@line) == 22){
          if ($_ =~ /^\d+/){
            $pct_ribosomal_bases = sprintf("%.2f", $line[10]*100);
            $pct_coding_bases = sprintf("%.2f", $line[11]*100);
            $pct_utr_bases = sprintf("%.2f", $line[12]*100);
            $pct_intronic_bases = sprintf("%.2f", $line[13]*100);
            $pct_intergenic_bases = sprintf("%.2f", $line[14]*100);
          }
        }
      }
      close(PIC);
    }else{
      $self->status_message("Could not find PicardRnaSeqMetrics.txt file for build: $build_id");
    }
    print $stats_fh "RNA-seq Percent Coding Bases\t$pct_coding_bases\t$common_name\tClinseq Build Summary\tPercent\tPercent of all mapped RNA-seq reads corresponding to coding regions for $common_name $extraction_type data\n";
    print $stats_fh "RNA-seq Percent UTR Bases\t$pct_utr_bases\t$common_name\tClinseq Build Summary\tPercent\tPercent of all mapped RNA-seq reads corresponding to UTR regions for $common_name $extraction_type data\n";
    print $stats_fh "RNA-seq Percent Intronic Bases\t$pct_intronic_bases\t$common_name\tClinseq Build Summary\tPercent\tPercent of all mapped RNA-seq reads corresponding to intronic regions for $common_name $extraction_type data\n";
    print $stats_fh "RNA-seq Percent Intergenic Bases\t$pct_intergenic_bases\t$common_name\tClinseq Build Summary\tPercent\tPercent of all mapped RNA-seq reads corresponding to intergenic regions for $common_name $extraction_type data\n";
    print $stats_fh "RNA-seq Percent Ribosomal Bases\t$pct_ribosomal_bases\t$common_name\tClinseq Build Summary\tPercent\tPercent of all mapped RNA-seq reads corresponding to coding ribosomal for $common_name $extraction_type data\n";
    $self->status_message("percent_coding_bases\tpercent_utr_bases\tpercent_intronic_bases\tpercent_intergenic_bases\tpercent_ribosomal_bases");
    $self->status_message("$pct_coding_bases\t$pct_utr_bases\t$pct_intronic_bases\t$pct_intergenic_bases\t$pct_ribosomal_bases");
    $self->status_message("build_id\n$build_id");
}

sub copy_from_rnaseq_build {
    my $self = shift;
    my $build_dir = shift;
    my $build_outdir = shift;
    my $common_name = shift;

    my @rnaseq_files_to_copy;
    push (@rnaseq_files_to_copy, "$build_dir/metrics/PicardRnaSeqMetrics.png");
    push (@rnaseq_files_to_copy, "$build_dir/metrics/PicardRnaSeqChart.pdf");
    push (@rnaseq_files_to_copy, "$build_dir/junctions/summary/PercentGeneJunctionsCovered_BreadthvsDepth_BoxPlot.pdf");
    push (@rnaseq_files_to_copy, "$build_dir/junctions/summary/ObservedJunctions_SpliceSiteAnchorTypes_Pie.pdf");
    push (@rnaseq_files_to_copy, "$build_dir/junctions/summary/ObservedJunctions_SpliceSiteUsage_Pie.pdf");
    push (@rnaseq_files_to_copy, "$build_dir/junctions/summary/TranscriptJunctionReadCounts_Log2_Hist.pdf");
    push (@rnaseq_files_to_copy, "$build_dir/bam-qc/*.pdf");
    push (@rnaseq_files_to_copy, "$build_dir/bam-qc/*.html");
    my $rnaseq_metrics_dir =  $build_outdir . "/rnaseq/$common_name/";
    Genome::Sys->shellcmd(cmd => "mkdir -p $rnaseq_metrics_dir");
    
    #Make copies of read locations .png and end bias plots for convenience
    foreach my $file (@rnaseq_files_to_copy){
      my $cp_cmd = "cp $file $rnaseq_metrics_dir";
      Genome::Sys->shellcmd(cmd => $cp_cmd, allow_failed_exit_code => 1);
    }
}

#Generate LIMS library quality reports (including alignment and quality metrics) for each flowcell associated with each sample
#e.g. illumina_info --sample H_KA-306905-S.4294 --report library --format tsvsub 
sub generate_LIMS_reports {
    my $self = shift;
    my $builds = shift;
    my $build_outdir = shift;
    my %samples_processed = ();
    $self->status_message("\n\nSample sequencing metrics from LIMS");
    $self->status_message("See results files in: $build_outdir\n");
    for my $build (@$builds) {
        next unless ($build->model->subject->class eq "Genome::Sample");
        my $subject = $build->subject;
        my $subject_name = $subject->name;
        my $common_name = $self->_get_subject_common_name($subject);
        #Only process each sample once
        unless ($samples_processed{$subject_name}){
            $samples_processed{$subject_name} = 1;
            my $lims_sample_outdir = $build_outdir . "/LIMS_reports/$common_name/";
            Genome::Sys->shellcmd(cmd=> "mkdir -p $lims_sample_outdir");
            $self->summarize_library_quality_reports_for_build($build, $lims_sample_outdir);
        }
    }
}

#Generate APIPE instrument data reports (including quality metrics) for each sample
#e.g. genome instrument-data list solexa --filter sample_name='H_LF-10-0372-09-131-1135122'  --show='id,flow_cell_id,lane,sample_name,library_name,read_length,is_paired_end,clusters,median_insert_size,sd_above_insert_size,target_region_set_name,fwd_filt_error_rate_avg,rev_filt_error_rate_avg' --style=csv
sub generate_APIPE_reports {
    my $self = shift;
    my $builds = shift;
    my $build_outdir = shift;
    my $stats_fh = shift;
    
    $self->status_message("\n\nSample sequencing metrics from APIPE");
    my %samples_processed;
    for my $build (@$builds) {
        next unless ($build->model->subject->class eq "Genome::Sample");
        my $subject = $build->subject;
        my $subject_name = $subject->name;
        my $common_name = $self->_get_subject_common_name($subject);
        #Only process each sample once
        unless ($samples_processed{$subject_name}){
            $samples_processed{$subject_name} = 1;
            my $apipe_sample_outdir = $build_outdir . "/APIPE_reports/$common_name/";
            Genome::Sys->shellcmd(cmd=> "mkdir -p $apipe_sample_outdir");
            $self->summarize_apipe_instrument_data_reports($build, $apipe_sample_outdir, $stats_fh);
        }
    }
}

sub summarize_apipe_instrument_data_reports {
    my $self = shift;
    my $build = shift;
    my $build_outdir = shift;
    my $stats_fh = shift;

    return unless $self->_is_reference_alignment_build($build);

    my $build_dir = $build->data_directory;
    my $tissue_desc = "[UNDEF tissue_desc]";
    my $extraction_type = "[UNDEF extraction_type]";
    my $subject = $build->subject;
    my $subject_name = $subject->name;
    my $common_name = $self->_get_subject_common_name($subject);
    if ($subject->can("tissue_desc")){
        if ($subject->tissue_desc){
            $tissue_desc = $subject->tissue_desc;
        }
    }
    if ($subject->can("extraction_type")){
        if ($subject->extraction_type){
            $extraction_type = $subject->extraction_type;
        }
    }

    #Occasionally the sample name will contain problem characters like "(" and ")" which need to be escaped or replaced
    my $subject_name_escaped = Genome::Utility::Text::sanitize_string_for_filesystem($subject_name);
    my $id_sample_summary_file_csv = $build_outdir . $subject_name_escaped . "_APIPE_Sample_Sequence_QC.csv";
    my $id_sample_summary_file_html = $build_outdir . $subject_name_escaped . "_APIPE_Sample_Sequence_QC.html";

    #Produce the sample sequencing summary in csv format
    $self->status_message("\n");
    $self->_run_solexa_lister($subject_name, 'csv', $id_sample_summary_file_csv);
    $self->_run_solexa_lister($subject_name, 'html', $id_sample_summary_file_html);
    $self->status_message("\nSample: $subject_name ($common_name | $tissue_desc | $extraction_type)");
    #Parse the csv file and store as tsv in the main output
    my ($insert_size_sum, $insert_size_count, $avg_insert_size, $fwd_error_rate_sum, $fwd_error_rate_count, $avg_fwd_error_rate, $rev_error_rate_sum, $rev_error_rate_count, $avg_rev_error_rate) = (0,0,0,0,0,0,0,0,0);
    if (-e $id_sample_summary_file_csv){
        open (SID, "$id_sample_summary_file_csv");
        while(<SID>){
            $_ =~ s/\,/\t/g;
            $self->status_message("$_");
            chomp($_);
            my @line = split("\t", $_);
            next unless (scalar(@line) == 13);
            if ($line[8] =~ /\d+/){
                $insert_size_sum += $line[8]; $insert_size_count++;
            }
            if ($line[11] =~ /\d+/){
                $fwd_error_rate_sum += $line[11]; $fwd_error_rate_count++;
            }
            if ($line[12] =~ /\d+/){
                $rev_error_rate_sum += $line[12]; $rev_error_rate_count++;
            }
        }
        close (SID);
    }
    if ($insert_size_sum){
        $avg_insert_size = sprintf("%.0f", $insert_size_sum/$insert_size_count);
    }else{
        $avg_insert_size = "n/a";
    }
    if ($fwd_error_rate_sum){
        $avg_fwd_error_rate = sprintf("%.2f", ($fwd_error_rate_sum/$fwd_error_rate_count));
    }else{
        $avg_fwd_error_rate = "n/a";
    }
    if ($rev_error_rate_sum){
        $avg_rev_error_rate = sprintf("%.2f", ($rev_error_rate_sum/$rev_error_rate_count));
    }else{
        $avg_rev_error_rate = "n/a";
    }
    print $stats_fh "Average Library Insert Size\t$avg_insert_size\t$common_name\tClinseq Build Summary\tAverage\tAverage of library-by-library median insert sizes for $common_name $extraction_type data\n";
    print $stats_fh "Forward Read Average Error Rate\t$avg_fwd_error_rate\t$common_name\tClinseq Build Summary\tAverage\tAverage of library-by-library sequencing forward read error rates for $common_name $extraction_type data\n";
    print $stats_fh "Reverse Read Average Error Rate\t$avg_rev_error_rate\t$common_name\tClinseq Build Summary\tAverage\tAverage of library-by-library sequencing reverse read error rates for $common_name $extraction_type data\n";
}

#Get the per-lane bam QC results for refalign builds.
sub get_perlane_bamqc_results {
    my $self = shift;
    my $build = shift;
    my $outdir = shift;
    my $data_types = shift;
    my $build_type = $self->_determine_wgs_or_exome_for_build($build, $data_types);
    $build_type = lc $build_type;
    my $subject = $build->subject;
    my $common_name = $self->_get_subject_common_name($subject);
    my $qc_dir = $outdir . "/$build_type/per_lane_bam_qc/$common_name";
    Genome::Sys->shellcmd(cmd => "mkdir -p $qc_dir");
    my $bam_qc_metrics = Genome::Model::ReferenceAlignment::Command::InstrumentDataAlignmentBams->create(
            build_id => $build->id, outdir => $qc_dir);
    my %lane_bamqc_path;
    $bam_qc_metrics->get_lane_bamqc_path($build, \%lane_bamqc_path);
    foreach my $lane  (keys(%lane_bamqc_path)) {
        if($lane_bamqc_path{$lane} ne "-") {
            my $perlane_bamqc_results_dir = $lane_bamqc_path{$lane};
            my $perlane_bamqc_op_dir = $qc_dir . "/lane" . $lane . "/";
            Genome::Sys->shellcmd(cmd => "mkdir -p $perlane_bamqc_op_dir");
            Genome::Sys->shellcmd(cmd => "cp -rf $perlane_bamqc_results_dir/* $perlane_bamqc_op_dir");
        }
    }
}

#Get the bam QC results for refalign builds.
sub get_bamqc_results {
    my $self = shift;
    my $build = shift;
    my $outdir = shift;
    my $data_types = shift;
    my $build_type = $self->_determine_wgs_or_exome_for_build($build, $data_types);
    $build_type = lc $build_type;
    my $subject = $build->subject;
    my $common_name = $self->_get_subject_common_name($subject);
    my $qc_dir = $outdir . "/$build_type/summary_bam_qc/$common_name";
    Genome::Sys->shellcmd(cmd => "mkdir -p $qc_dir");
    my $bam_qc_metrics = Genome::Model::ReferenceAlignment::Command::BamQcMetrics->create(
            build_id => $build->id, output_directory => $qc_dir);
    $bam_qc_metrics->execute();
}

sub summarize_library_quality_reports_for_build {
    my $self = shift;
    my $build = shift;
    my $outdir = shift;

    return unless $self->_is_reference_alignment_build($build);
    return if $self->skip_lims_reports;

    my @formats = qw (csv tsv html);
    my @reports = qw (lane library library_index_summary run);
    my %rf;
    $rf{lane}{file} = "lane";
    $rf{library}{file} = "library_index_summary";
    $rf{library_index_summary}{file} = "library_index_summary";
    $rf{run}{file} = "run";

    my $build_dir = $build->data_directory;
    my $tissue_desc = "[UNDEF tissue_desc]";
    my $extraction_type = "[UNDEF extraction_type]";
    my $subject = $build->subject;
    my $subject_name = $subject->name;
    my $common_name = $self->_get_subject_common_name($subject);
    if ($subject->can("tissue_desc")){
        $tissue_desc = $subject->tissue_desc;
    }
    if ($subject->can("extraction_type")){
        $extraction_type = $subject->extraction_type;
    }

    my $previous_cwd = Cwd::getcwd();
    chdir($outdir); #illumina_info writes to cwd

    foreach my $format (@formats){
        foreach my $report (@reports){
            my $rf_file = $rf{$report}{file};
            my $id_summary_file = $outdir . $subject_name . "_LIMS_Sample_Sequence_QC_" . $report . "." . $format;
            next if (-e $id_summary_file); #Shortcut for testing purposes
            my $tmp_file = $outdir . $rf_file . "." . $format;
            my $id_list_cmd = "PATH=/gsc/bin/:\$PATH illumina_info --sample $subject_name --report $report --format $format";
            #my $id_list_cmd = "PATH=/gsc/bin/:\$PATH illumina_info --sample $subject_name --report $report --format $format 1>/dev/null 2>/dev/null";
            unless ($self->skip_lims_reports){
                Genome::Sys->shellcmd(cmd => $id_list_cmd, allow_failed_exit_code => 1);
                if (-e $tmp_file){
                    Genome::Sys->move_file($tmp_file, $id_summary_file);
                }
            }
        }
    }
    chdir($previous_cwd);
    return 1;
}

sub _get_subject_common_name {
    my $self = shift;
    my $subject = shift;
    my $common_name = "[UNDEF common_name]";
    if ($subject->can("common_name")){
        if ($subject->common_name){
            $common_name = $subject->common_name;
        }
    }
    return $common_name;
}

sub _run_solexa_lister {
    my $self = shift;
    my $sample_name = shift;
    my $style = shift;
    my $output_file = shift;

    my $output_fh = Genome::Sys->open_file_for_writing($output_file);
    my $id_list_cmd1 = Genome::InstrumentData::Command::List::Solexa->create(
        filter => "sample_name='$sample_name'",
        show => 'id,flow_cell_id,lane,sample_name,library_name,read_length,is_paired_end,clusters,median_insert_size,sd_above_insert_size,target_region_set_name,fwd_filt_error_rate_avg,rev_filt_error_rate_avg',
        style => $style,
        output => $output_fh,
    );
    $id_list_cmd1->execute() or die 'failed to run lister to generate ' . $output_file;
    $output_fh->close();

    return 1;
}

sub _set_data_types {
    my $self = shift;
    my $clinseq_build = shift;
    my $data_types = shift;
    my $wgs_somvar_build = $clinseq_build->wgs_build;
    my $exome_somvar_build = $clinseq_build->exome_build;
    my $tumor_rnaseq_build = $clinseq_build->tumor_rnaseq_build;
    my $normal_rnaseq_build = $clinseq_build->normal_rnaseq_build;
    my $wgs_normal_refalign_build = $wgs_somvar_build->normal_build if ($wgs_somvar_build);
    my $wgs_tumor_refalign_build = $wgs_somvar_build->tumor_build if ($wgs_somvar_build);
    my $exome_normal_refalign_build = $exome_somvar_build->normal_build if ($exome_somvar_build);
    my $exome_tumor_refalign_build = $exome_somvar_build->tumor_build if ($exome_somvar_build);
    if($wgs_tumor_refalign_build) {
        $data_types->{$wgs_tumor_refalign_build->id} = "WGS";
    }
    if($wgs_normal_refalign_build) {
        $data_types->{$wgs_normal_refalign_build->id} = "WGS";
    }
    if($exome_tumor_refalign_build) {
        $data_types->{$exome_tumor_refalign_build->id} = "Exome";
    }
    if($exome_normal_refalign_build) {
        $data_types->{$exome_normal_refalign_build->id} = "Exome";
    }
    if($tumor_rnaseq_build) {
        $data_types->{$tumor_rnaseq_build->id} = "RNAseq";
    }
    if($normal_rnaseq_build) {
        $data_types->{$normal_rnaseq_build->id} = "RNAseq";
    }
}

sub _determine_wgs_or_exome_for_build {
    my $self = shift;
    my $build = shift;
    my $data_types = shift;
    my $m = $build->model;
    my $id = $build->id;
    if(defined $data_types->{$id}) {
        return $data_types->{$id};
    }
    
    my $pp = $m->processing_profile;
    my $pp_name = $pp->name;
    my $data_type = "Unknown";
    if ($pp_name =~ /wgs/i){
        $data_type = "WGS";
    }elsif($pp_name =~ /exome/i){
        $data_type = "Exome";
    }

    return $data_type;
}

sub _determine_wgs_or_exome_for_instrument_data {
    my $self = shift;
    my @lanes = shift;

    my $lane_count = scalar(@lanes);
    #Infer whether the data is wgs or exome based on the presence of target region set names
    my $trsn_count = 0;
    for my $id (@lanes) {
        if ($id->target_region_set_name){
            $trsn_count++;
        }
    }
    my $sequence_type = "unknown";
    if ($trsn_count == $lane_count){
        $sequence_type = "exome";
    }elsif($trsn_count == 0){
        $sequence_type = "wgs";
    }else{
        $sequence_type = "mixed";
    }

    return $sequence_type;
}

sub _is_reference_alignment_build {
    my $self = shift;
    my $build = shift;

    return $build->type_name eq 'reference alignment';
}

sub _is_rna_seq_build {
    my $self = shift;
    my $build = shift;

    return $build->type_name eq 'rna seq';
}

sub _is_somatic_variation_build {
    my $self = shift;
    my $build = shift;

    return $build->type_name eq 'somatic variation';
}

1;

