package Genome::Model::ClinSeq;

use strict;
use warnings;
use Genome;
use List::MoreUtils;

# these are used below, and are also used in the documentation on commands in the tree
# to provide the most useful examples possible
our $DEFAULT_CANCER_ANNOTATION_DB_ID = 'tgi/cancer-annotation/human/build37-20150205.1';
our $DEFAULT_MISC_ANNOTATION_DB_ID   = 'tgi/misc-annotation/human/build37-20130113.1';
our $DEFAULT_COSMIC_ANNOTATION_DB_ID = 'cosmic/65.3';

class Genome::Model::ClinSeq {
    is                 => ['Genome::Model', 'Genome::Model::ClinSeq::Util',],
    has_optional_input => [
        wgs_model => {
            is  => 'Genome::Model::SomaticVariation',
            doc => 'somatic variation model for wgs data'
        },
        exome_model => {
            is  => 'Genome::Model::SomaticVariation',
            doc => 'somatic variation model for exome data'
        },
        tumor_rnaseq_model => {
            is  => 'Genome::Model::RnaSeq',
            doc => 'rnaseq model for tumor rna-seq data'
        },
        normal_rnaseq_model => {
            is  => 'Genome::Model::RnaSeq',
            doc => 'rnaseq model for normal rna-seq data'
        },
        de_model => {
            is  => 'Genome::Model::DifferentialExpression',
            doc => 'differential-expression for tumor vs normal rna-seq data'
        },

        cancer_annotation_db => {
            is            => 'Genome::Db::Tgi::CancerAnnotation',
            default_value => $DEFAULT_CANCER_ANNOTATION_DB_ID
        },
        misc_annotation_db => {
            is            => 'Genome::Db::Tgi::MiscAnnotation',
            default_value => $DEFAULT_MISC_ANNOTATION_DB_ID
        },
        cosmic_annotation_db => {
            is            => 'Genome::Db::Cosmic',
            default_value => $DEFAULT_COSMIC_ANNOTATION_DB_ID
        },

        force => {
            is  => 'Boolean',
            doc => 'skip sanity checks on input models'
        },

        #processing_profile      => { is => 'Genome::ProcessingProfile::ClinSeq', id_by => 'processing_profile_id', default_value => { } },
    ],
    #Processing profile parameters
    #If you add/modify params here make sure to add/modify params in
    # Genome/Test/Factory/ProcessingProfile/ClinSeq.pm
    has_param => [
        bam_readcount_version => {
            is  => 'Text',
            doc => 'The bam readcount version to use during clonality analysis'
        },
        sireport_min_tumor_vaf => {
            is  => 'Number',
            doc => 'Variants with a tumor VAF less than this (in any tumor sample) will be filtered out.'
        },
        sireport_max_normal_vaf => {
            is  => 'Number',
            doc => 'Variants with a normal VAF greater than this (in any normal sample) will be filtered out.'
        },
        sireport_min_coverage => {
            is  => 'Number',
            doc => 'Variants with coverage less than this (in any sample) will be filtered out.'
        },
        sireport_min_mq_bq => {
            is  => 'Text',
            doc => 'Comma separated increasing-list of minimum mapping qualities for bam-readcounting.'
        },
        exome_cnv => {
            is  => 'Boolean',
            doc => 'Should we run the Exome-cnv step or not.'
        },
    ],
    has_optional_metric => [
        common_name => {
            is  => 'Text',
            doc => 'the name chosen for the root directory in the build'
        },
    ],
    has_calculated => [
        expected_common_name => {
            is             => 'Text',
            calculate_from => [qw /wgs_model exome_model tumor_rnaseq_model normal_rnaseq_model/],
            calculate      => q|
              my ($wgs_common_name, $exome_common_name, $tumor_rnaseq_common_name, $normal_rnaseq_common_name, $wgs_name, $exome_name, $tumor_rnaseq_name, $normal_rnaseq_name);
              if ($wgs_model) {
                  $wgs_common_name = $wgs_model->subject->individual->common_name;
                  $wgs_name = $wgs_model->subject->individual->name;
              }
              if ($exome_model) {
                  $exome_common_name = $exome_model->subject->individual->common_name;
                  $exome_name = $exome_model->subject->individual->name;
              }
              if ($tumor_rnaseq_model) {
                  $tumor_rnaseq_common_name = $tumor_rnaseq_model->subject->individual->common_name;
                  $tumor_rnaseq_name = $tumor_rnaseq_model->subject->individual->name;
              }
              if ($normal_rnaseq_model) {
                  $normal_rnaseq_common_name = $normal_rnaseq_model->subject->individual->common_name;
                  $normal_rnaseq_name = $normal_rnaseq_model->subject->individual->name;
              }
              my @names = ($wgs_common_name, $exome_common_name, $tumor_rnaseq_common_name, $normal_rnaseq_common_name, $wgs_name, $exome_name, $tumor_rnaseq_name, $normal_rnaseq_name);
              my $final_name = "UnknownName";
              foreach my $name (@names){
                if ($name){
                  $final_name = $name;
                  last();
                }
              }
              return $final_name;
            |
        },
    ],
    has => [
        individual_common_name => {
            via => '__self__',
            to  => 'expected_common_name',
            doc => 'alias for benefit of Solr indexing',
        },
    ],
    doc => 'clinical and discovery sequencing data analysis and convergence of RNASeq, WGS and exome capture data',
};

sub define_by {'Genome::Model::Command::Define::BaseMinimal'}

sub _help_synopsis {
    my $self = shift;
    return <<"EOS"
    genome processing-profile create clin-seq  --name 'November 2011 Clinical Sequencing'
    genome model define clin-seq  --processing-profile='November 2011 Clinical Sequencing'  --wgs-model=2882504846 --exome-model=2882505032 --tumor-rnaseq-model=2880794613
    # Automatically builds if/when the models have a complete underlying build
EOS
}

sub help_detail_for_create_profile {
    return <<EOS
The ClinSeq pipeline has the following processing-profile parameters.
  bam_readcount_version => version of bamreadcounts to use for read-counting steps,
  sireport_min_tumor_vaf => minimum tumor VAF cutoff to filter out variants in SnvIndelReport,
  sireport_max_normal_vaf   => maximum normal VAF cutoff to filter out variants in SnvIndelReport,
  sireport_min_coverage => minimum coverage cutoff to filter out variants in SnvIndelReport,
  sireport_min_mq_bq => semicolon delimited list of comma separated mapping quality and base quality scores, these
    are used as cutoffs for read-counting in SnvIndelReport.
EOS
}

sub _help_detail_for_model_define {
    return <<EOS

The ClinSeq pipeline takes four models, each of which is optional, and produces data sets potentially useful in a clinical or discovery project.

There are several primary goals:

1. Generate results that help to identify clinically actionable events (gain of function mutations, amplifications, etc.)

2. Converge results across data types for a single case (e.g. variant read counts from wgs + exome + rnaseq)

3. Automatically generate statistics, tables and figures to be used in manuscripts

4. Test novel analysis methods in a pipeline that runs faster and is more agile than say somatic-variation

EOS
}

sub _resolve_subject {
    my $self     = shift;
    my @subjects = $self->_infer_candidate_subjects_from_input_models();
    if (@subjects > 1) {
        if ($self->force) {
            @subjects = ($subjects[0]);
        }
        else {
            die $self->error_message(
                "Conflicting subjects on input models!:\n\t" . join("\n\t", map {$_->__display_name__} @subjects));
        }
    }
    elsif (@subjects == 0) {
        die $self->error_message("No subjects on input models?");
    }
    return $subjects[0];
}

sub _resolve_annotation {
    my $self        = shift;
    my @annotations = $self->_infer_annotations_from_input_models();
    if (@annotations > 1) {
        if ($self->force) {
            @annotations = ($annotations[0]);
        }
        else {
            die $self->error_message("Conflicting annotations on input models!:\n\t"
                    . join("\n\t", map {$_->__display_name__} @annotations));
        }
    }
    elsif (@annotations == 0) {
        die $self->error_message("No annotation builds on input models?");
    }
    return $annotations[0];
}

sub _resolve_reference {
    my $self       = shift;
    my @references = $self->_infer_references_from_input_models();
    if (@references > 1) {
        if ($self->force) {
            @references = ($references[0]);
        }
        else {
            die $self->error_message("Conflicting reference sequence builds on input models!:\n\t"
                    . join("\n\t", map {$_->__display_name__} @references));
        }
    }
    elsif (@references == 0) {
        die $self->error_message("No reference builds on input models?");
    }
    return $references[0];
}

#Implement specific error checking here, any error that is added to the @errors array will prevent the model from being commited to the database
#Could also implement a --force input above to allow over-riding of errors
sub __errors__ {
    my $self   = shift;
    my @errors = $self->SUPER::__errors__;
    return @errors;
}

sub _initialize_build {
    my $self  = shift;
    my $build = shift;

    # this is currently tracked as a metric on the build but it should really be an output/metric
    my $common_name = $self->expected_common_name;
    $build->common_name($common_name);

    return 1;
}

#MAP WORKFLOW INPUTS
sub map_workflow_inputs {
    my $self  = shift;
    my $build = shift;

    my $model          = $self;
    my $data_directory = $build->data_directory;

    # inputs
    my $wgs_build           = $build->wgs_build;
    my $exome_build         = $build->exome_build;
    my $tumor_rnaseq_build  = $build->tumor_rnaseq_build;
    my $normal_rnaseq_build = $build->normal_rnaseq_build;

    # set during build initialization
    my $common_name = $build->common_name;
    unless ($common_name) {
        die "no common name set on build during initialization?";
    }

    my $cancer_annotation_db = $build->cancer_annotation_db;
    my $misc_annotation_db   = $build->misc_annotation_db;
    my $cosmic_annotation_db = $build->cosmic_annotation_db;
    my ($mqs, $bqs) = $self->parse_qualities;

    # initial inputs used for various steps
    my @inputs = (
        model                   => $model,
        build                   => $build,
        wgs_build               => $wgs_build,
        exome_build             => $exome_build,
        tumor_rnaseq_build      => $tumor_rnaseq_build,
        normal_rnaseq_build     => $normal_rnaseq_build,
        working_dir             => $data_directory,
        common_name             => $common_name,
        verbose                 => 1,
        cancer_annotation_db    => $cancer_annotation_db,
        misc_annotation_db      => $misc_annotation_db,
        cosmic_annotation_db    => $cosmic_annotation_db,
        bam_readcount_version   => $self->bam_readcount_version,
        sireport_min_tumor_vaf  => $self->sireport_min_tumor_vaf,
        sireport_max_normal_vaf => $self->sireport_max_normal_vaf,
        sireport_min_coverage   => $self->sireport_min_coverage,
    );

    my $annotation_build = $self->_resolve_annotation;
    push @inputs, annotation_build => $annotation_build;

    #Identify test model
    if ($self->name =~ /^apipe\-test/) {
        push @inputs, test => 1;
    }

    my $patient_dir = $data_directory . "/" . $common_name;
    my @dirs        = ($patient_dir);

    #SummarizeBuilds
    my $input_summary_dir = $patient_dir . "/input";
    push @dirs,   $input_summary_dir;
    push @inputs, summarize_builds_outdir => $input_summary_dir;
    push @inputs, summarize_builds_log_file => $input_summary_dir . "/QC_Report.tsv";

    if ($build->id) {
        push @inputs, summarize_builds_skip_lims_reports => 0;
    }
    else {
        #Watch out for -ve build IDs which will occur when the ClinSeq.t test is run.  In that case, do not run the LIMS reports
        push @inputs, summarize_builds_skip_lims_reports => 1;
    }

    #DumpIgvXml
    my $igv_session_dir = $patient_dir . '/igv';
    push @dirs, $igv_session_dir;
    push @inputs, igv_session_dir => $igv_session_dir;

    #GetVariantSources
    if ($wgs_build or $exome_build) {
        my $variant_sources_dir = $patient_dir . '/variant_source_callers';
        push @dirs, $variant_sources_dir;
        if ($wgs_build) {
            my $wgs_variant_sources_dir = $variant_sources_dir . '/wgs';
            push @inputs, wgs_variant_sources_dir => $wgs_variant_sources_dir;
            push @dirs, $wgs_variant_sources_dir;
        }
        if ($exome_build) {
            my $exome_variant_sources_dir = $variant_sources_dir . '/exome';
            push @inputs, exome_variant_sources_dir => $exome_variant_sources_dir;
            push @dirs, $exome_variant_sources_dir;
        }
    }

    #ImportSnvsIndels
    push @inputs, import_snvs_indels_filter_mt => 1;
    push @inputs, import_snvs_indels_outdir    => $patient_dir;

    #CreateMutationDiagrams
    if ($wgs_build or $exome_build) {
        my $mutation_diagram_dir = $patient_dir . '/mutation_diagrams';
        push @dirs, $mutation_diagram_dir;
        push @inputs,
            (
            mutation_diagram_outdir              => $mutation_diagram_dir,
            mutation_diagram_collapse_variants   => 1,
            mutation_diagram_max_snvs_per_file   => 1500,
            mutation_diagram_max_indels_per_file => 1500,
            );
    }

    #rnaseq analysis steps according to what rna-seq builds are defined as inputs
    if ($normal_rnaseq_build or $tumor_rnaseq_build) {
        my $rnaseq_dir = $patient_dir . '/rnaseq';
        push @dirs, $rnaseq_dir;
        push @inputs, cufflinks_percent_cutoff => 1;

        #TophatJunctionsAbsolute and CufflinksExpressionAbsolute for 'normal' sample
        if ($normal_rnaseq_build) {
            my $normal_rnaseq_dir = $rnaseq_dir . '/normal';
            push @dirs, $normal_rnaseq_dir;
            my $normal_tophat_junctions_absolute_dir = $normal_rnaseq_dir . '/tophat_junctions_absolute';
            push @dirs, $normal_tophat_junctions_absolute_dir;
            push @inputs, normal_tophat_junctions_absolute_dir => $normal_tophat_junctions_absolute_dir;
            my $normal_cufflinks_expression_absolute_dir = $normal_rnaseq_dir . '/cufflinks_expression_absolute';
            push @dirs, $normal_cufflinks_expression_absolute_dir;
            push @inputs, normal_cufflinks_expression_absolute_dir => $normal_cufflinks_expression_absolute_dir;
        }

        #TophatJunctionsAbsolute and CufflinksExpressionAbsolute for 'tumor' sample
        if ($tumor_rnaseq_build) {
            my $tumor_rnaseq_dir = $rnaseq_dir . '/tumor';
            push @dirs, $tumor_rnaseq_dir;
            my $tumor_tophat_junctions_absolute_dir = $tumor_rnaseq_dir . '/tophat_junctions_absolute';
            push @dirs, $tumor_tophat_junctions_absolute_dir;
            push @inputs, tumor_tophat_junctions_absolute_dir => $tumor_tophat_junctions_absolute_dir;
            my $tumor_cufflinks_expression_absolute_dir = $tumor_rnaseq_dir . '/cufflinks_expression_absolute';
            push @dirs, $tumor_cufflinks_expression_absolute_dir;
            push @inputs, tumor_cufflinks_expression_absolute_dir => $tumor_cufflinks_expression_absolute_dir;
        }

        #CufflinksDifferentialExpression
        if ($normal_rnaseq_build and $tumor_rnaseq_build) {
            my $cufflinks_differential_expression_dir = $rnaseq_dir . '/cufflinks_differential_expression';
            push @dirs, $cufflinks_differential_expression_dir;
            push @inputs, cufflinks_differential_expression_dir => $cufflinks_differential_expression_dir;
        }

        #Filtered and Intersected ChimeraScan fusion output for tumor RNAseq
        #The following will only work if the rna-seq build used a processing profile that uses chimerascan.
        if ($tumor_rnaseq_build) {
            #Check for ChimeraScan fusion results
            if (-e $tumor_rnaseq_build->data_directory . '/fusions/filtered_chimeras.bedpe') {
                #copy over fusion files to this dir even if SV calls do not exist.
                my $tumor_filtered_fusion_dir = $patient_dir . '/rnaseq/tumor/fusions';
                push @dirs, $tumor_filtered_fusion_dir;
                if ($wgs_build) {
                    #Check for SV calls file
                    if (-e $wgs_build->data_directory . '/effects/svs.hq.annotated') {
                        my $ncbi_human_ensembl_build_id = $tumor_rnaseq_build->annotation_build->id;
                        my $tumor_filtered_fusion_file  = $tumor_filtered_fusion_dir . '/filtered_chimeras.bedpe';
                        my $wgs_sv_file = $build->wgs_build->data_directory . '/effects/svs.hq.annotated';
                        my $tumor_filtered_intersected_fusion_file =
                            $tumor_filtered_fusion_dir . '/chimeras.filtered.intersected.bedpe';
                        push @inputs, ncbi_human_ensembl_build_id            => $ncbi_human_ensembl_build_id;
                        push @inputs, wgs_sv_file                            => $wgs_sv_file;
                        push @inputs, tumor_filtered_fusion_file             => $tumor_filtered_fusion_file;
                        push @inputs, tumor_filtered_intersected_fusion_file => $tumor_filtered_intersected_fusion_file;
                    }
                }
            }
        }
    }

    #GenerateClonalityPlots
    if ($wgs_build) {
        my $clonality_dir = $patient_dir . "/clonality/";
        push @dirs, $clonality_dir;
        push @inputs, clonality_dir => $clonality_dir;
    }

    #Make base-dir for CNV's
    if ($wgs_build or $exome_build or $self->has_microarray_build()) {
        my $cnv_dir = $patient_dir . "/cnv/";
        push @dirs, $cnv_dir;
        push @inputs, cnv_dir => $cnv_dir;
    }

    #Make base-dir for WGS CNV's
    if ($wgs_build) {
        my $wgs_cnv_dir = $patient_dir . "/cnv/wgs_cnv";
        push @dirs, $wgs_cnv_dir;
        push @inputs, wgs_cnv_dir => $wgs_cnv_dir;
    }

    #RunMicroArrayCnView
    if ($self->has_microarray_build()) {
        my $microarray_cnv_dir = $patient_dir . "/cnv/microarray_cnv/";
        push @dirs, $microarray_cnv_dir;
        push @inputs, microarray_cnv_dir => $microarray_cnv_dir;
    }

    #ExomeCNV
    if ($build->should_run_exome_cnv) {
        my $exome_cnv_dir = $patient_dir . "/cnv/exome_cnv/";
        push @dirs, $exome_cnv_dir;
        push @inputs, exome_cnv_dir => $exome_cnv_dir;
    }

    my $docm_variants_file = $self->_get_docm_variants_file($self->cancer_annotation_db);
    #Make DOCM report
    if (($exome_build or $wgs_build) and $docm_variants_file) {
        my $docm_report_dir = $patient_dir . "/docm_report/";
        push @dirs,   $docm_report_dir;
        push @inputs, docm_report_dir => $docm_report_dir;
        push @inputs, docm_variants_file => $docm_variants_file;
        push @inputs, docmreport_min_bq => 0;
        push @inputs, docmreport_min_mq => 1;
    }

    #SummarizeSvs
    if ($wgs_build) {
        my $sv_dir = $patient_dir . "/sv/";
        push @dirs, $sv_dir;
        push @inputs, sv_dir => $sv_dir;
    }

    #Make base-dir for WGS CNV's
    if ($wgs_build) {
        my $wgs_cnv_summary_dir = $patient_dir . "/cnv/wgs_cnv/summary/";
        push @dirs, $wgs_cnv_summary_dir;
        push @inputs, wgs_cnv_summary_dir => $wgs_cnv_summary_dir;
    }

    #CreateMutationSpectrum
    my $iterator = List::MoreUtils::each_arrayref([1 .. @$mqs], $mqs, $bqs);
    while (my ($i, $mq, $bq) = $iterator->()) {
        if ($wgs_build) {
            push @inputs,
                'wgs_mutation_spectrum_outdir' . $i => $patient_dir . "/mutation-spectrum/b" . $bq . "_q" . $mq . "/";
            push @inputs, 'wgs_mutation_spectrum_datatype' . $i => 'wgs';
        }
        if ($exome_build) {
            push @inputs,
                'exome_mutation_spectrum_outdir' . $i => $patient_dir . "/mutation-spectrum/b" . $bq . "_q" . $mq . "/";
            push @inputs, 'exome_mutation_spectrum_datatype' . $i => 'exome';
        }
    }
    if ($wgs_build or $exome_build) {
        my $mutation_spectrum_dir = $patient_dir . "/mutation-spectrum";
        push @dirs, $mutation_spectrum_dir;
    }

    #AnnotateGenesByCategory
    push @inputs, 'gene_name_columns' => ['mapped_gene_name'];
    push @inputs, 'gene_name_regex'   => 'mapped_gene_name';

    #MakeCircosPlot
    if ($wgs_build) {
        my $circos_dir = $patient_dir . "/circos";
        push @dirs, $circos_dir;
        push @inputs, circos_outdir => $circos_dir;
    }

    #Converge SnvIndelReport and SciClone
    if ($exome_build || $wgs_build) {
        my $iterator = List::MoreUtils::each_arrayref([1 .. @$mqs], $mqs, $bqs);
        while (my ($i, $mq, $bq) = $iterator->()) {
            my $snv_indel_report_dir1 =
                File::Spec->join($patient_dir, "/snv_indel_report/" . "b" . $bq . "_" . "q" . $mq);
            push @dirs, $snv_indel_report_dir1;
            push @inputs, "snv_indel_report_dir" . $i => $snv_indel_report_dir1;
            if ($wgs_build or $build->should_run_exome_cnv) {
                my $sciclone_dir1 =
                    File::Spec->join($patient_dir, "/clonality/sciclone/" . "b" . $bq . "_" . "q" . $mq);
                push @dirs, $sciclone_dir1;
                push @inputs, "sciclone_dir" . $i => $sciclone_dir1;
            }
            push @inputs, "sireport_min_bq" . $i => $bq;
            push @inputs, "sireport_min_mq" . $i => $mq;
        }
        my $target_gene_list = $cancer_annotation_db->data_directory . "/CancerGeneCensus/cancer_gene_census_ensgs.tsv";
        my $target_gene_list_name = "CancerGeneCensus";

        push @inputs, snv_indel_report_clean                 => 1;
        push @inputs, snv_indel_report_tmp_space             => 1;
        push @inputs, snv_indel_report_target_gene_list      => $target_gene_list;
        push @inputs, snv_indel_report_target_gene_list_name => $target_gene_list_name;
        push @inputs, snv_indel_report_tiers                 => 'tier1';
    }

    #IdentifyLoh
    if ($exome_build || $wgs_build) {
        my $loh_dir = $patient_dir . "/loh/";
        push @dirs, $loh_dir;
        push @inputs, loh_output_dir => $loh_dir;
    }

    # For now it works to create directories here because the data_directory has been allocated.
    #It is possible that this would not happen until later, which would mess up assigning inputs to many of the commands.
    for my $dir (@dirs) {
        Genome::Sys->create_directory($dir);
    }

    return @inputs;
}


sub _resolve_workflow_for_build {
    # This is called by Genome::Model::Build::start()
    # Returns a Workflow::Operation
    my $self        = shift;
    my $build       = shift;
    my $lsf_queue   = shift;  # TODO: the workflow shouldn't need this yet
    my $lsf_project = shift;
    my ($mqs, $bqs) = $self->parse_qualities;
    if (!defined $lsf_queue || $lsf_queue eq '' || $lsf_queue eq 'inline') {
        $lsf_queue = Genome::Config::get('lsf_queue_build_worker_alt');
    }
    if (!defined $lsf_project || $lsf_project eq '') {
        $lsf_project = 'build' . $build->id;
    }

    #We need to call this so that the subdirectories get created for the test
    my %inputs           = $self->map_workflow_inputs($build);

    # Make the workflow
    my $workflow = Genome::WorkflowBuilder::DAG->create(
        name    => $build->workflow_name,
        log_dir => $build->log_directory,
    );

    #SummarizeBuilds - Summarize build inputs using SummarizeBuilds.pm
    my $summarize_builds_op = $self->summarize_builds_op($workflow);

    #ImportSnvsIndels - Import SNVs and Indels
    my $import_snvs_indels_op;
    if ($build->wgs_build or $build->exome_build) {
        $import_snvs_indels_op = $self->import_snvs_indels_op($workflow);
        if ($build->wgs_build) {
            $workflow->connect_input(
                input_property       => 'wgs_build',
                destination          => $import_snvs_indels_op,
                destination_property => 'wgs_build',
            );
        }
        if ($build->exome_build) {
            $workflow->connect_input(
                input_property       => 'exome_build',
                destination          => $import_snvs_indels_op,
                destination_property => 'exome_build',
            );
        }
    }

    #GetVariantSources - Determine source variant caller for SNVs and InDels for wgs data
    my $wgs_variant_sources_op;
    if ($build->wgs_build) {
        $wgs_variant_sources_op = $self->variant_sources_op($workflow, 'wgs');
    }
    #GetVariantSources - Determine source variant caller for SNVs and InDels for exome data
    my $exome_variant_sources_op;
    if ($build->exome_build) {
        $exome_variant_sources_op = $self->variant_sources_op($workflow, 'exome');
    }

    #CreateMutationDiagrams - Create mutation diagrams (lolliplots) for all Tier1 SNVs/Indels and compare to COSMIC SNVs/Indels
    if ($build->wgs_build or $build->exome_build) {
        my $mutation_diagram_op = $self->mutation_diagram_op($workflow);
        if ($build->wgs_build and $build->exome_build) {
            my $wgs_exome_build_converge_op = $self->wgs_exome_build_converge_op($workflow);
            $workflow->connect_input(
                input_property       => 'exome_variant_sources_dir',
                destination          => $exome_variant_sources_op,
                destination_property => 'outdir',
            );
            $workflow->create_link(
                source               => $wgs_exome_build_converge_op,
                source_property      => 'builds',
                destination          => $mutation_diagram_op,
                destination_property => 'builds',
            );
        }
        elsif ($build->wgs_build) {
            $workflow->connect_input(
                input_property       => 'wgs_build',
                destination          => $mutation_diagram_op,
                destination_property => 'builds',
            );
        }
        elsif ($build->exome_build) {
            $workflow->connect_input(
                input_property       => 'exome_build',
                destination          => $mutation_diagram_op,
                destination_property => 'builds',
            );
        }
    }

    #TophatJunctionsAbsolute - Run tophat junctions absolute analysis on normal
    #CufflinksExpressionAbsolute - Run cufflinks expression absolute analysis on normal
    if ($build->normal_rnaseq_build) {
        my $normal_tophat_junctions_absolute_op = $self->tophat_junctions_absolute_op($workflow, 'normal');
        my $normal_cufflinks_expression_absolute_op = $self->cufflinks_expression_absolute_op($workflow, 'normal');
    }
    #TophatJunctionsAbsolute - Run tophat junctions absolute analysis on tumor
    #CufflinksExpressionAbsolute - Run cufflinks expression absolute analysis on tumor
    #Intersect filtered fusion calls with WGS SV calls.
    my ($tumor_tophat_junctions_absolute_op, $tumor_cufflinks_expression_absolute_op);
    if ($build->tumor_rnaseq_build) {
        $tumor_tophat_junctions_absolute_op = $self->tophat_junctions_absolute_op($workflow, 'tumor');
        $tumor_cufflinks_expression_absolute_op = $self->cufflinks_expression_absolute_op($workflow, 'tumor');
        if (-e $build->tumor_rnaseq_build->data_directory . '/fusions/filtered_chimeras.bedpe') {
            #copy over fusion files
            $self->copy_fusion_files($build);
            if ($build->wgs_build) {
                if (-e $build->wgs_build->data_directory . '/effects/svs.hq.annotated') {
                    my $intersect_tumor_fusion_sv_op = $self->intersect_tumor_fusion_sv_op($workflow);
                }
            }
        }
    }

    #CufflinksDifferentialExpression - Run cufflinks differential expression
    my $cufflinks_differential_expression_op;
    if ($build->normal_rnaseq_build and $build->tumor_rnaseq_build) {
        $cufflinks_differential_expression_op = $self->cufflinks_differential_expression_op($workflow);
    }

    #DumpIgvXml - Create IGV xml session files with increasing numbers of tracks and store in a single (WGS and Exome BAM files, RNA-seq BAM files, junctions.bed, SNV bed files, etc.)
    #genome model clin-seq dump-igv-xml --outdir=/gscuser/mgriffit/ --builds=119971814
    my $igv_session_op = $self->igv_session_op($workflow);

    #GenerateClonalityPlots - Run clonality analysis and produce clonality plots
    #RunCnView - Produce copy number results with run-cn-view.
    #SummarizeCnvs - Generate a summary of CNV results, copy cnvs.hq, cnvs.png, single-bam copy number plot PDF, etc. to the cnv directory
    my ($clonality_op, $run_cn_view_op);
    if ($build->wgs_build) {
        $clonality_op = $self->clonality_op($workflow);
        $run_cn_view_op = $self->run_cn_view_op($workflow);
        for my $property (qw(cnv_hmm_file cnv_hq_file)) {
            $workflow->create_link(
                source               => $clonality_op,
                source_property      => $property,
                destination          => $run_cn_view_op,
                destination_property => $property,
            );
        }
        my $summarize_cnvs_op = $self->summarize_cnvs_op($workflow);
        for my $property (qw(cnv_hmm_file cnv_hq_file)) {
            $workflow->create_link(
                source               => $clonality_op,
                source_property      => $property,
                destination          => $summarize_cnvs_op,
                destination_property => $property,
            );
        }
        for my $property (qw(gene_amp_file gene_del_file)) {
            $workflow->create_link(
                source               => $run_cn_view_op,
                source_property      => $property,
                destination          => $summarize_cnvs_op,
                destination_property => $property,
            );
        }
    }

    #RunMicroarrayCNV - produce cnv plots using microarray results
    my $microarray_cnv_op;
    if ($self->has_microarray_build()) {
        $microarray_cnv_op = $self->microarray_cnv_op($workflow);
    }

    #RunExomeCNV - produce cnv plots using WEx results
    my $exome_cnv_op;
    if ($build->should_run_exome_cnv) {
        $exome_cnv_op = $self->exome_cnv_op($workflow);
    }

    #RunDOCMReport
    if (($build->exome_build or $build->wgs_build)
        and $self->_get_docm_variants_file($self->cancer_annotation_db))
    {
        my $docm_report_op = $self->docm_report_op($workflow);
    }

    #SummarizeSvs - Generate a summary of SV results from the WGS SV results
    my $summarize_svs_op;
    if ($build->wgs_build) {
        $summarize_svs_op = $self->summarize_svs_op($workflow);
    }

    #Add gene category annotations to some output files from steps above. (e.g. determine which SNV affected genes are kinases, ion channels, etc.)
    #AnnotateGenesByCategory - gene_category_<exome|wgs|wgs_exome>_<snv|indel>_result
    for my $sequencing_type (qw(exome wgs wgs_exome)) {
        my $build_accessor = "${sequencing_type}_build";
        if ($build->$build_accessor) {
            for my $variant_type (qw(snv indel)) {
                my $annotate_genes_by_category_op = $self->annotate_genes_by_category_op($workflow, $sequencing_type, $variant_type);
                $workflow->create_link(
                    source               => $import_snvs_indels_op,
                    source_property      => "${sequencing_type}_${variant_type}_file",
                    destination          => $annotate_genes_by_category_op,
                    destination_property => 'infile',
                );
            }
        }
    }
    #AnnotateGenesByCategory - gene_category_cnv_<amp|del|ampdel>_result
    if ($build->wgs_build) {
        for my $type (qw(amp del ampdel)) {
            my $annotate_genes_by_category_cnv_op = $self->annotate_genes_by_category_cnv_op($workflow, $type);
            $workflow->create_link(
                source               => $run_cn_view_op,
                source_property      => "gene_${type}_file",
                destination          => $annotate_genes_by_category_cnv_op,
                destination_property => 'infile',
            );
        }
    }
    #AnnotateGenesByCategory - gene_category_<cufflinks|tophat>_result
    if ($build->tumor_rnaseq_build) {
        my $annotate_genes_by_category_cufflinks_op = $self->annotate_genes_by_category_tumor_rna_op($workflow, 'cufflinks');
        $workflow->create_link(
            source               => $tumor_cufflinks_expression_absolute_op,
            source_property      => 'tumor_fpkm_topnpercent_file',
            destination          => $annotate_genes_by_category_cufflinks_op,
            destination_property => 'infile',
        );
        my $annotate_genes_by_category_tophat_op = $self->annotate_genes_by_category_tumor_rna_op($workflow, 'tophat');
        $workflow->create_link(
            source               => $tumor_tophat_junctions_absolute_op,
            source_property      => 'junction_topnpercent_file',
            destination          => $annotate_genes_by_category_tophat_op,
            destination_property => 'infile',
        );
    }
    #AnnotateGenesByCategory - gene_category_coding_<up|down|de>_result
    if ($build->normal_rnaseq_build && $build->tumor_rnaseq_build) {
        my %types = (
            up   => 'up-regulated',
            down => 'down-regulated',
            de   => 'DE',
        );
        while (my ($type_short, $type) = each %types) {
            my $annotate_genes_by_category_normal_rna_op = $self->annotate_genes_by_category_normal_rna_op($workflow, $type, $type_short);
            $workflow->create_link(
                source               => $cufflinks_differential_expression_op,
                source_property      => "coding_hq_${type_short}_file",
                destination          => $annotate_genes_by_category_normal_rna_op,
                destination_property => 'infile',
            );
        }
    }

    #DGIDB gene annotation
    for my $sequencing_type (qw(exome wgs wgs_exome)) {
        my $build_accessor = "${sequencing_type}_build";
        if ($build->$build_accessor) {
            for my $variant_type (qw(snv indel)) {
                my $annotate_genes_by_dgidb_op = $self->annotate_genes_by_dgidb_op($workflow, $sequencing_type, $variant_type);
                $workflow->create_link(
                    source               => $import_snvs_indels_op,
                    source_property      => "${sequencing_type}_${variant_type}_file",
                    destination          => $annotate_genes_by_dgidb_op,
                    destination_property => 'input_file',
                );
            }
        }
    }

    if ($build->wgs_build) {
        my $annotate_genes_by_dgidb_cnv_amp_op = $self->_annotate_genes_by_dgidb_op($workflow, 'gene_amp_file', 'dgidb_cnv_amp_result');
        $workflow->create_link(
            source               => $run_cn_view_op,
            source_property      => 'gene_amp_file',
            destination          => $annotate_genes_by_dgidb_cnv_amp_op,
            destination_property => 'input_file',
        );
        my $annotate_genes_by_dgidb_sv_fusion_op = $self->_annotate_genes_by_dgidb_op($workflow, 'fusion_output_file', 'dgidb_sv_fusion_result');
        $workflow->create_link(
            source               => $summarize_svs_op,
            source_property      => 'fusion_output_file',
            destination          => $annotate_genes_by_dgidb_sv_fusion_op,
            destination_property => 'input_file',
        );
    }

    if ($build->tumor_rnaseq_build) {
        my $annotate_genes_by_dgidb_cufflink_op = $self->_annotate_genes_by_dgidb_op($workflow, 'tumor_fpkm_topnpercent_file', 'dgidb_cufflinks_result');
        $workflow->create_link(
            source               => $tumor_cufflinks_expression_absolute_op,
            source_property      => 'tumor_fpkm_topnpercent_file',
            destination          => $annotate_genes_by_dgidb_cufflink_op,
            destination_property => 'input_file',
        );
        my $annotate_genes_by_dgidb_tophat_op = $self->_annotate_genes_by_dgidb_op($workflow, 'junction_topnpercent_file', 'dgidb_tophat_result');
        $workflow->create_link(
            source               => $tumor_tophat_junctions_absolute_op,
            source_property      => 'junction_topnpercent_file',
            destination          => $annotate_genes_by_dgidb_tophat_op,
            destination_property => 'input_file',
        );
    }

    #SummarizeTier1SnvSupport - For each of the following: WGS SNVs, Exome SNVs, and WGS+Exome SNVs, do the following:
    #Get BAM readcounts for WGS (tumor/normal), Exome (tumor/normal), RNAseq (tumor), RNAseq (normal) - as available of course
    #TODO: Break this down to do direct calls to GetBamReadCounts instead of wrapping it.
    for my $sequencing_type (qw(wgs exome wgs_exome)) {
        my $build_accessor = "${sequencing_type}_build";
        if ($build->$build_accessor) {
            my $summarize_tier1_snv_support_op = $self->summarize_tier1_snv_support_op($workflow, $sequencing_type);
            $workflow->create_link(
                source               => $import_snvs_indels_op,
                source_property      => "${sequencing_type}_snv_file",
                destination          => $summarize_tier1_snv_support_op,
                destination_property => "${sequencing_type}_positions_file",
            );
            if ($build->tumor_rnaseq_build) {
                $workflow->create_link(
                    source               => $tumor_cufflinks_expression_absolute_op,
                    source_property      => 'tumor_fpkm_file',
                    destination          => $summarize_tier1_snv_support_op,
                    destination_property => 'tumor_fpkm_file',
                );
            }
        }
    }

    #MakeCircosPlot - Creates a Circos plot to summarize the data using MakeCircosPlot.pm
    #Currently WGS data is a minimum requirement for Circos plot generation.
    if ($build->wgs_build) {
        my $make_circos_plot_op = $self->make_circos_plot_op($workflow);
        $workflow->create_link(
            source               => $summarize_svs_op,
            source_property      => 'fusion_output_file',
            destination          => $make_circos_plot_op,
            destination_property => 'candidate_fusion_infile',
        );
        $workflow->create_link(
            source               => $clonality_op,
            source_property      => 'cnv_hmm_file',
            destination          => $make_circos_plot_op,
            destination_property => 'cnv_hmm_file',
        );
        if ($build->normal_rnaseq_build) {
            $workflow->create_link(
                source               => $cufflinks_differential_expression_op,
                source_property      => 'coding_hq_de_file',
                destination          => $make_circos_plot_op,
                destination_property => 'coding_hq_de_file',
            );
        }
        elsif ($build->tumor_rnaseq_build) {
            $workflow->create_link(
                source               => $tumor_cufflinks_expression_absolute_op,
                source_property      => 'tumor_fpkm_topnpercent_file',
                destination          => $make_circos_plot_op,
                destination_property => 'tumor_fpkm_topnpercent_file',
            );
        }
        $workflow->create_link(
            source               => $import_snvs_indels_op,
            source_property      => 'result',
            destination          => $make_circos_plot_op,
            destination_property => 'import_snvs_indels_result',
        );
        $workflow->create_link(
            source               => $run_cn_view_op,
            source_property      => 'gene_ampdel_file',
            destination          => $make_circos_plot_op,
            destination_property => 'gene_ampdel_file',
        );
    }

    #Converge SnvIndelReport
    #CreateMutationDiagrams - Create mutation spectrum results for wgs and exome data
    #GenerateSciClonePlots - Run clonality analysis and produce clonality plots
    if ($build->wgs_build || $build->exome_build) {
        my %variant_sources_ops = (
            exome => $exome_variant_sources_op,
            wgs   => $wgs_variant_sources_op,
        );
        #Create a report for each $bq $mq combo.
        my $iterator = List::MoreUtils::each_arrayref([1 .. @$mqs], $mqs, $bqs);
        while (my ($i, $mq, $bq) = $iterator->()) {
            my $converge_snv_indel_report_op = $self->converge_snv_indel_report_op($workflow, $i);
            for my $sequencing_type (qw(exome wgs)) {
                my $build_accessor = "${sequencing_type}_build";
                if ($build->$build_accessor) {
                    for my $variant_type (qw(snv indel)) {
                        $workflow->create_link(
                            source               => $variant_sources_ops{$sequencing_type},
                            source_property      => "${variant_type}_variant_sources_file",
                            destination          => $converge_snv_indel_report_op,
                            destination_property => "_${sequencing_type}_${variant_type}_variant_sources_file",
                        );
                    }
                }
                my $create_mutation_spectrum_op = $self->create_mutation_spectrum_op($workflow, $sequencing_type, $i);
                $workflow->create_link(
                    source               => $converge_snv_indel_report_op,
                    source_property      => 'result',
                    destination          => $create_mutation_spectrum_op,
                    destination_property => 'converge_snv_indel_report_result',
                );
            }
            if ($build->wgs_build or $build->should_run_exome_cnv) {
                my $sciclone_op = $self->sciclone_op($workflow, $i);
                if ($build->wgs_build) {
                    $workflow->create_link(
                        source               => $run_cn_view_op,
                        source_property      => 'result',
                        destination          => $sciclone_op,
                        destination_property => 'wgs_cnv_result',
                    );
                }
                if ($build->should_run_exome_cnv) {
                    $workflow->create_link(
                        source               => $exome_cnv_op,
                        source_property      => 'result',
                        destination          => $sciclone_op,
                        destination_property => 'exome_cnv_result',
                    );
                }
                if ($self->has_microarray_build()) {
                    $workflow->create_link(
                        source               => $microarray_cnv_op,
                        source_property      => 'result',
                        destination          => $sciclone_op,
                        destination_property => 'microarray_cnv_result',
                    );
                }
                $workflow->create_link(
                    source               =>  $converge_snv_indel_report_op,
                    source_property      => 'result',
                    destination          => $sciclone_op,
                    destination_property => 'converge_snv_indel_report_result',
                );
            }
        }
    }

    #IdentifyLoh - Run identify-loh tool for exome or WGS data
    #genome model clin-seq identify-loh --clinseq-build=fafd219665d54462893fbacfe6639f70 --outdir=/Documents/GTB11/ --bamrc-version=0.7
    if ($build->wgs_build or $build->exome_build) {
        my $identify_loh_op = $self->identify_loh_op($workflow);
    }

    return $workflow;
}

sub summarize_builds_op {
    my $self = shift;
    my $workflow = shift;

    my $summarize_builds_op = Genome::WorkflowBuilder::Command->create(
        name => 'Creating a summary of input builds using summarize-builds',
        command => 'Genome::Model::ClinSeq::Command::SummarizeBuilds',
    );
    $workflow->add_operation($summarize_builds_op);
    $workflow->connect_input(
        input_property       => 'build',
        destination          => $summarize_builds_op,
        destination_property => 'builds',
    );
    $workflow->connect_input(
        input_property       => 'summarize_builds_outdir',
        destination          => $summarize_builds_op,
        destination_property => 'outdir',
    );
    $workflow->connect_input(
        input_property       => 'summarize_builds_skip_lims_reports',
        destination          => $summarize_builds_op,
        destination_property => 'skip_lims_reports',
    );
    $workflow->connect_input(
        input_property       => 'summarize_builds_log_file',
        destination          => $summarize_builds_op,
        destination_property => 'log_file',
    );
    $workflow->connect_output(
        output_property => 'summarize_builds_result',
        source          => $summarize_builds_op,
        source_property => 'result',
    );

    return $summarize_builds_op;
}

sub import_snvs_indels_op {
    my $self = shift;
    my $workflow = shift;

    my $import_snvs_indels_op = Genome::WorkflowBuilder::Command->create(
        name => 'Importing SNVs and Indels from somatic results, parsing, and merging exome/wgs if possible',
        command => 'Genome::Model::ClinSeq::Command::ImportSnvsIndels',
    );
    $workflow->add_operation($import_snvs_indels_op);
    $workflow->connect_input(
        input_property       => 'cancer_annotation_db',
        destination          => $import_snvs_indels_op,
        destination_property => 'cancer_annotation_db',
    );
    $workflow->connect_input(
        input_property       => 'import_snvs_indels_outdir',
        destination          => $import_snvs_indels_op,
        destination_property => 'outdir',
    );
    $workflow->connect_input(
        input_property       => 'import_snvs_indels_filter_mt',
        destination          => $import_snvs_indels_op,
        destination_property => 'filter_mt',
    );
    $workflow->connect_output(
        output_property => 'import_snvs_indels_result',
        source          => $import_snvs_indels_op,
        source_property => 'result',
    );

    return $import_snvs_indels_op;
}

sub variant_sources_op {
    my $self = shift;
    my $workflow = shift;
    my $type = shift;

    my $variant_sources_op = Genome::WorkflowBuilder::Command->create(
        name => "Determining source variant callers of all tier1-3 SNVs and InDels for $type data",
        command => 'Genome::Model::ClinSeq::Command::GetVariantSources',
    );
    $workflow->add_operation($variant_sources_op);
    $workflow->connect_input(
        input_property       => "${type}_build",
        destination          => $variant_sources_op,
        destination_property => 'builds',
    );
    $workflow->connect_input(
        input_property       => "${type}_variant_sources_dir",
        destination          => $variant_sources_op,
        destination_property => 'outdir',
    );
    $workflow->connect_output(
        output_property => "${type}_variant_sources_result",
        source          => $variant_sources_op,
        source_property => 'result',
    );

    return $variant_sources_op;
}

sub mutation_diagram_op {
    my $self = shift;
    my $workflow = shift;

    my $mutation_diagram_op = Genome::WorkflowBuilder::Command->create(
        name => 'Creating mutation-diagram plots',
        command => 'Genome::Model::ClinSeq::Command::CreateMutationDiagrams',
    );
    $workflow->add_operation($mutation_diagram_op);
    for my $property (qw(cancer_annotation_db cosmic_annotation_db)) {
        $workflow->connect_input(
            input_property => $property,
            destination => $mutation_diagram_op,
            destination_property => $property,
        );
    }
    for my $property (qw/outdir collapse_variants max_snvs_per_file max_indels_per_file/) {
        $workflow->connect_input(
            input_property       => "mutation_diagram_$property",
            destination          => $mutation_diagram_op,
            destination_property => $property,
        );
    }
    $workflow->connect_output(
        output_property => 'mutation_diagram_result',
        source          => $mutation_diagram_op,
        source_property => 'result',
    );

    return $mutation_diagram_op;
}

sub wgs_exome_build_converge_op {
    my $self = shift;
    my $workflow = shift;

    my $wgs_exome_build_converge_op = Genome::WorkflowBuilder::Converge->create(
        name => 'Converge wgs and exome builds',
        input_properties => ['wgs_build', 'exome_build'],
        output_properties => ['builds'],
    );
    $workflow->add_operation($wgs_exome_build_converge_op);
    for my $build (qw(wgs_build exome_build)) {
        $workflow->connect_input(
            input_property => $build,
            destination => $wgs_exome_build_converge_op,
            destination_property => $build,
        );
    }

    return $wgs_exome_build_converge_op;
}

sub tophat_junctions_absolute_op {
    my $self = shift;
    my $workflow = shift;
    my $type = shift;

    my $tophat_junctions_absolute_op = Genome::WorkflowBuilder::Command->create(
        name => "Performing tophat junction expression absolute analysis for $type sample",
        command => 'Genome::Model::ClinSeq::Command::TophatJunctionsAbsolute',
    );
    $workflow->add_operation($tophat_junctions_absolute_op);
    $workflow->connect_input(
        input_property       => 'cancer_annotation_db',
        destination          => $tophat_junctions_absolute_op,
        destination_property => 'cancer_annotation_db',
    );
    $workflow->connect_input(
        input_property       => "${type}_rnaseq_build",
        destination          => $tophat_junctions_absolute_op,
        destination_property => 'build',
    );
    $workflow->connect_input(
        input_property       => "${type}_tophat_junctions_absolute_dir",
        destination          => $tophat_junctions_absolute_op,
        destination_property => 'outdir',
    );
    $workflow->connect_output(
        output_property => "${type}tophat_junctions_absolute_result",
        source          => $tophat_junctions_absolute_op,
        source_property => 'result',
    );

    return $tophat_junctions_absolute_op;
}

sub cufflinks_expression_absolute_op {
    my $self = shift;
    my $workflow = shift;
    my $type = shift;

    my $cufflinks_expression_absolute_op = Genome::WorkflowBuilder::Command->create(
        name => "Performing cufflinks expression absolute analysis for $type sample",
        command => 'Genome::Model::ClinSeq::Command::CufflinksExpressionAbsolute',
    );
    $workflow->add_operation($cufflinks_expression_absolute_op);
    $workflow->connect_input(
        input_property       => 'cancer_annotation_db',
        destination          => $cufflinks_expression_absolute_op,
        destination_property => 'cancer_annotation_db',
    );
    $workflow->connect_input(
        input_property       => "${type}_rnaseq_build",
        destination          => $cufflinks_expression_absolute_op,
        destination_property => 'build',
    );
    $workflow->connect_input(
        input_property       => "${type}_cufflinks_expression_absolute_dir",
        destination          => $cufflinks_expression_absolute_op,
        destination_property => 'outdir',
    );
    $workflow->connect_input(
        input_property       => 'cufflinks_percent_cutoff',
        destination          => $cufflinks_expression_absolute_op,
        destination_property => 'percent_cutoff',
    );
    $workflow->connect_output(
        output_property => "${type}_cufflinks_expression_absolute_result",
        source          => $cufflinks_expression_absolute_op,
        source_property => 'result',
    );

    return $cufflinks_expression_absolute_op;
}

sub cufflinks_differential_expression_op {
    my $self = shift;
    my $workflow = shift;

    my $cufflinks_differential_expression_op = Genome::WorkflowBuilder::Command->create(
        name => 'Performing cufflinks differential expression analysis for case vs. control (e.g. tumor vs. normal)',
        command => 'Genome::Model::ClinSeq::Command::CufflinksDifferentialExpression',
    );
    $workflow->add_operation($cufflinks_differential_expression_op);
    $workflow->connect_input(
        input_property       => 'normal_rnaseq_build',
        destination          => $cufflinks_differential_expression_op,
        destination_property => 'control_build',
    );
    $workflow->connect_input(
        input_property       => 'tumor_rnaseq_build',
        destination          => $cufflinks_differential_expression_op,
        destination_property => 'case_build',
    );
    $workflow->connect_input(
        input_property       => 'cufflinks_differential_expression_dir',
        destination          => $cufflinks_differential_expression_op,
        destination_property => 'outdir',
    );
    $workflow->connect_output(
        output_property => 'cufflinks_differential_expression_result',
        source          => $cufflinks_differential_expression_op,
        source_property => 'result',
    );

    return $cufflinks_differential_expression_op;
}

sub intersect_tumor_fusion_sv_op {
    my $self = shift;
    my $workflow = shift;

    my $intersect_tumor_fusion_sv_op = Genome::WorkflowBuilder::Command->create(
        name => 'Intersecting filtered tumor ChimeraScan fusion calls with WGS SV calls',
        command => 'Genome::Model::Tools::ChimeraScan::IntersectSv',
    );
    $workflow->add_operation($intersect_tumor_fusion_sv_op);
    $workflow->connect_input(
        input_property       => 'ncbi_human_ensembl_build_id',
        destination          => $intersect_tumor_fusion_sv_op,
        destination_property => 'annotation_build_id',
    );
    $workflow->connect_input(
        input_property       => 'tumor_filtered_intersected_fusion_file',
        destination          => $intersect_tumor_fusion_sv_op,
        destination_property => 'output_file',
    );
    $workflow->connect_input(
        input_property       => 'wgs_sv_file',
        destination          => $intersect_tumor_fusion_sv_op,
        destination_property => 'sv_output_file',
    );
    $workflow->connect_input(
        input_property       => 'tumor_filtered_fusion_file',
        destination          => $intersect_tumor_fusion_sv_op,
        destination_property => 'filtered_bedpe_file',
    );
    $workflow->connect_output(
        output_property => 'intersect_tumor_fusion_sv_result',
        source          => $intersect_tumor_fusion_sv_op,
        source_property => 'result',
    );

    return $intersect_tumor_fusion_sv_op;
}

sub igv_session_op {
    my $self = shift;
    my $workflow = shift;

    my $igv_session_op = Genome::WorkflowBuilder::Command->create(
        name => 'Create IGV XML session files for varying levels of detail using the input builds',
        command => 'Genome::Model::ClinSeq::Command::DumpIgvXml',
    );
    $workflow->add_operation($igv_session_op);
    $workflow->connect_input(
        input_property       => 'build',
        destination          => $igv_session_op,
        destination_property => 'builds',
    );
    $workflow->connect_input(
        input_property       => 'igv_session_dir',
        destination          => $igv_session_op,
        destination_property => 'outdir',
    );
    $workflow->connect_output(
        output_property => 'igv_session_result',
        source          => $igv_session_op,
        source_property => 'result',
    );

    return $igv_session_op;
}

sub clonality_op {
    my $self = shift;
    my $workflow = shift;

    my $clonality_op = Genome::WorkflowBuilder::Command->create(
        name => 'Run clonality analysis and produce clonality plots',
        command => 'Genome::Model::ClinSeq::Command::GenerateClonalityPlots',
    );
    $workflow->add_operation($clonality_op);
    $workflow->connect_input(
        input_property       => 'misc_annotation_db',
        destination          => $clonality_op,
        destination_property => 'misc_annotation_db',
    );
    $workflow->connect_input(
        input_property       => 'wgs_build',
        destination          => $clonality_op,
        destination_property => 'somatic_var_build',
    );
    $workflow->connect_input(
        input_property       => 'clonality_dir',
        destination          => $clonality_op,
        destination_property => 'output_dir',
    );
    for my $property (qw(common_name bam_readcount_version)) {
        $workflow->connect_input(
            input_property       => $property,
            destination          => $clonality_op,
            destination_property => $property,
        );
    }
    $workflow->connect_output(
        output_property => 'clonality_result',
        source          => $clonality_op,
        source_property => 'result',
    );

    return $clonality_op;
}

sub run_cn_view_op {
    my $self = shift;
    my $workflow = shift;

    my $run_cn_view_op = Genome::WorkflowBuilder::Command->create(
        name => 'Use gmt copy-number cn-view to produce copy number tables and images',
        command => 'Genome::Model::ClinSeq::Command::RunCnView',
    );
    $workflow->add_operation($run_cn_view_op);
    $workflow->connect_input(
        input_property       => 'cancer_annotation_db',
        destination          => $run_cn_view_op,
        destination_property => 'cancer_annotation_db',
    );
    $workflow->connect_input(
        input_property       => 'wgs_build',
        destination          => $run_cn_view_op,
        destination_property => 'build',
    );
    $workflow->connect_input(
        input_property       => 'wgs_cnv_dir',
        destination          => $run_cn_view_op,
        destination_property => 'outdir',
    );
    $workflow->connect_output(
        output_property => 'run_cn_view_result',
        source          => $run_cn_view_op,
        source_property => 'result',
    );

    return $run_cn_view_op;
}

sub microarray_cnv_op {
    my $self = shift;
    my $workflow = shift;

    my $microarray_cnv_op = Genome::WorkflowBuilder::Command->create(
        name => 'Call somatic copy number changes using microarray calls',
        command => 'Genome::Model::ClinSeq::Command::MicroarrayCnv',
    );
    $workflow->add_operation($microarray_cnv_op);
    $workflow->connect_input(
        input_property       => 'microarray_cnv_dir',
        destination          => $microarray_cnv_op,
        destination_property => 'outdir',
    );
    $workflow->connect_input(
        input_property       => 'model',
        destination          => $microarray_cnv_op,
        destination_property => 'clinseq_model',
    );
    $workflow->connect_input(
        input_property       => 'annotation_build',
        destination          => $microarray_cnv_op,
        destination_property => 'annotation_build_id',
    );
    $workflow->connect_output(
        output_property => 'microarray_cnv_result',
        source          => $microarray_cnv_op,
        source_property => 'result',
    );

    return $microarray_cnv_op;
}

sub exome_cnv_op {
    my $self = shift;
    my $workflow = shift;

    my $exome_cnv_op = Genome::WorkflowBuilder::Command->create(
        name => 'Call somatic copy number changes using exome data',
        command => 'Genome::Model::Tools::CopyNumber::Cnmops',
    );
    $workflow->add_operation($exome_cnv_op);
    $workflow->connect_input(
        input_property       => 'exome_cnv_dir',
        destination          => $exome_cnv_op,
        destination_property => 'outdir',
    );
    $workflow->connect_input(
        input_property       => 'model',
        destination          => $exome_cnv_op,
        destination_property => 'clinseq_model',
    );
    $workflow->connect_input(
        input_property       => 'annotation_build',
        destination          => $exome_cnv_op,
        destination_property => 'annotation_build_id',
    );
    $workflow->connect_output(
        output_property => 'exome_cnv_result',
        source          => $exome_cnv_op,
        source_property => 'result',
    );

    return $exome_cnv_op;
}

sub docm_report_op {
    my $self = shift;
    my $workflow = shift;

    my $docm_report_op = Genome::WorkflowBuilder::Command->create(
        name => 'Produce a report using DOCM',
        command => 'Genome::Model::ClinSeq::Command::Converge::DocmReport',
    );
    $workflow->add_operation($docm_report_op);
    $workflow->connect_input(
        input_property       => 'docm_report_dir',
        destination          => $docm_report_op,
        destination_property => 'outdir',
    );
    $workflow->connect_input(
        input_property       => 'build',
        destination          => $docm_report_op,
        destination_property => 'builds',
    );
    for my $property (qw(docm_variants_file bam_readcount_version)) {
        $workflow->connect_input(
            input_property       => $property,
            destination          => $docm_report_op,
            destination_property => $property,
        );
    }
    for my $property (qw(bq mq)) {
        $workflow->connect_input(
            input_property       => "docmreport_min_$property",
            destination          => $docm_report_op,
            destination_property => $property,
        );
    }
    $workflow->connect_output(
        output_property => 'docm_report_result',
        source          => $docm_report_op,
        source_property => 'result',
    );

    return $docm_report_op;
}

sub summarize_cnvs_op {
    my $self = shift;
    my $workflow = shift;

    my $summarize_cnvs_op = Genome::WorkflowBuilder::Command->create(
        name => 'Summarize CNV results from WGS somatic variation',
        command => 'Genome::Model::ClinSeq::Command::SummarizeCnvs',
    );
    $workflow->add_operation($summarize_cnvs_op);
    $workflow->connect_input(
        input_property       => 'wgs_cnv_summary_dir',
        destination          => $summarize_cnvs_op,
        destination_property => 'outdir',
    );
    $workflow->connect_input(
        input_property       => 'wgs_build',
        destination          => $summarize_cnvs_op,
        destination_property => 'build',
    );
    $workflow->connect_output(
        output_property => 'summarize_cnvs_result',
        source          => $summarize_cnvs_op,
        source_property => 'result',
    );

    return $summarize_cnvs_op;
}

sub summarize_svs_op {
    my $self = shift;
    my $workflow = shift;

    my $summarize_svs_op = Genome::WorkflowBuilder::Command->create(
        name => 'Summarize SV results from WGS somatic variation',
        command => 'Genome::Model::ClinSeq::Command::SummarizeSvs',
    );
    $workflow->add_operation($summarize_svs_op);
    $workflow->connect_input(
        input_property       => 'cancer_annotation_db',
        destination          => $summarize_svs_op,
        destination_property => 'cancer_annotation_db',
    );
    $workflow->connect_input(
        input_property       => 'sv_dir',
        destination          => $summarize_svs_op,
        destination_property => 'outdir',
    );
    $workflow->connect_input(
        input_property       => 'wgs_build',
        destination          => $summarize_svs_op,
        destination_property => 'builds',
    );
    $workflow->connect_output(
        output_property => 'summarize_svs_result',
        source          => $summarize_svs_op,
        source_property => 'result',
    );

    return $summarize_svs_op;
}

sub annotate_genes_by_category_op {
    my $self = shift;
    my $workflow = shift;
    my $sequencing_type = shift;
    my $variant_type = shift;

    my $name = "Add gene category annotations to ${variant_type}s identified by $sequencing_type";
    my $output_property = "gene_category_${sequencing_type}_${variant_type}_result";

    return $self->_annotate_genes_by_category_op($workflow, $name, $output_property);
}

sub annotate_genes_by_category_cnv_op {
    my $self = shift;
    my $workflow = shift;
    my $type = shift;

    my $name = "Add gene category annotations to CNV $type genes";
    my $output_property = "gene_category_cnv_${type}_result";

    return $self->_annotate_genes_by_category_op($workflow, $name, $output_property);
}

sub annotate_genes_by_category_tumor_rna_op {
    my $self = shift;
    my $workflow = shift;
    my $software = shift;

    my $name = "Add gene category annotations to $software top1 percent genes";
    my $output_property = "gene_category_${software}_result";

    return $self->_annotate_genes_by_category_op($workflow, $name, $output_property);
}

sub annotate_genes_by_category_normal_rna_op {
    my $self = shift;
    my $workflow = shift;
    my $type = shift;
    my $type_short = shift;

    my $name = "Add gene category annotations to $type coding genes";
    my $output_property = "gene_category_coding_${type_short}_result";

    return $self->_annotate_genes_by_category_op($workflow, $name, $output_property);
}

sub _annotate_genes_by_category_op {
    my $self = shift;
    my $workflow = shift;
    my $name = shift;
    my $output_property = shift;

    my $annotate_genes_by_category_op = Genome::WorkflowBuilder::Command->create(
        name => $name,
        command => 'Genome::Model::ClinSeq::Command::AnnotateGenesByCategory',
    );
    $workflow->add_operation($annotate_genes_by_category_op);
    for my $property (qw(gene_name_columns cancer_annotation_db)) {
        $workflow->connect_input(
            input_property       => $property,
            destination          => $annotate_genes_by_category_op,
            destination_property => $property,
        );
    }
    $workflow->connect_output(
        output_property => $output_property,
        source          => $annotate_genes_by_category_op,
        source_property => 'result',
    );

    return $annotate_genes_by_category_op;
}

sub annotate_genes_by_dgidb_op {
    my $self = shift;
    my $workflow = shift;
    my $sequencing_type = shift;
    my $variant_type = shift;

    return $self->_annotate_genes_by_dgidb_op($workflow, "${sequencing_type}_${variant_type}_file", "dgidb_${sequencing_type}_${variant_type}_result");
}

sub _annotate_genes_by_dgidb_op {
    my $self = shift;
    my $workflow = shift;
    my $file_name = shift;
    my $output_property = shift;

    my $annotate_genes_by_dgidb_op = Genome::WorkflowBuilder::Command->create(
        name => "Add dgidb gene annotations to $file_name",
        command => 'Genome::Model::ClinSeq::Command::AnnotateGenesByDgidb',
    );
    $workflow->add_operation($annotate_genes_by_dgidb_op);
    $workflow->connect_input(
        input_property       => 'gene_name_regex',
        destination          => $annotate_genes_by_dgidb_op,
        destination_property => 'gene_name_regex',
    );
    $workflow->connect_output(
        output_property => $output_property,
        source          => $annotate_genes_by_dgidb_op,
        source_property => 'result',
    );

    return $annotate_genes_by_dgidb_op;
}

sub summarize_tier1_snv_support_op {
    my $self = shift;
    my $workflow = shift;
    my $sequencing_type = shift;

    my %sequencing_types = (
        wgs       => 'WGS',
        exome     => 'Exome',
        wgs_exome => 'WGS and Exome'
    );
    my $txt_name = $sequencing_types{$sequencing_type};
    my $summarize_tier1_snv_support_op = Genome::WorkflowBuilder::Command->create(
        name => "$txt_name Summarize Tier 1 SNV Support (BAM read counts)",
        command => 'Genome::Model::ClinSeq::Command::SummarizeTier1SnvSupport',
    );
    $workflow->add_operation($summarize_tier1_snv_support_op);
    for my $property (qw(wgs_build exome_build tumor_rnaseq_build normal_rnaseq_build verbose cancer_annotation_db)) {
        $workflow->connect_input(
            input_property       => $property,
            destination          => $summarize_tier1_snv_support_op,
            destination_property => $property,
        );
    }
    $workflow->connect_output(
        output_property => "summarize_${sequencing_type}_tier1_snv_support_result",
        source          => $summarize_tier1_snv_support_op,
        source_property => 'result',
    );

    return $summarize_tier1_snv_support_op;
}

sub make_circos_plot_op {
    my $self = shift;
    my $workflow = shift;

    my $make_circos_plot_op = Genome::WorkflowBuilder::Command->create(
        name => 'Create a Circos plot using MakeCircosPlot',
        command => 'Genome::Model::ClinSeq::Command::MakeCircosPlot',
    );
    $workflow->add_operation($make_circos_plot_op);
    $workflow->connect_input(
        input_property       => 'build',
        destination          => $make_circos_plot_op,
        destination_property => 'build',
    );
    $workflow->connect_input(
        input_property       => 'circos_outdir',
        destination          => $make_circos_plot_op,
        destination_property => 'output_directory',
    );
    $workflow->connect_output(
        output_property => 'circos_result',
        source          => $make_circos_plot_op,
        source_property => 'result',
    );

    return $make_circos_plot_op;
}

sub converge_snv_indel_report_op {
    my $self = shift;
    my $workflow = shift;
    my $iterator = shift;

    my $converge_snv_indel_report_op = Genome::WorkflowBuilder::Command->create(
        name => "Generate SnvIndel Report$iterator",
        command => 'Genome::Model::ClinSeq::Command::Converge::SnvIndelReport',
    );
    $workflow->add_operation($converge_snv_indel_report_op);
    $workflow->connect_input(
        input_property       => 'build',
        destination          => $converge_snv_indel_report_op,
        destination_property => 'builds',
    );
    for my $property (qw(bam_readcount_version annotation_build)) {
        $workflow->connect_input(
            input_property       => $property,
            destination          => $converge_snv_indel_report_op,
            destination_property => $property,
        );
    }
    $workflow->connect_input(
        input_property       => "snv_indel_report_dir$iterator",
        destination          => $converge_snv_indel_report_op,
        destination_property => 'outdir',
    );
    for my $property (qw(clean tmp_space target_gene_list target_gene_list_name)) {
        $workflow->connect_input(
            input_property       => "snv_indel_report_$property",
            destination          => $converge_snv_indel_report_op,
            destination_property => $property,
        );
    }
    for my $property (qw(min_tumor_vaf max_normal_vaf min_coverage)) {
        $workflow->connect_input(
            input_property       => "sireport_$property",
            destination          => $converge_snv_indel_report_op,
            destination_property => $property,
        );
    }
    for my $property (qw(bq mq)) {
        $workflow->connect_input(
            input_property       => "sireport_min_$property$iterator",
            destination          => $converge_snv_indel_report_op,
            destination_property => $property,
        );
    }
    #If this is a build of a test model, perform a faster analysis (e.g. apipe-test-clinseq-wer)
    if ($self->name =~ /^apipe\-test/) {
        $workflow->connect_input(
            input_property       => 'snv_indel_report_tiers',
            destination          => $converge_snv_indel_report_op,
            destination_property => 'tiers',
        );
    }
    $workflow->connect_output(
        output_property => "converge_snv_indel_report_result$iterator",
        source          => $converge_snv_indel_report_op,
        source_property => 'result',
    );

    return $converge_snv_indel_report_op;
}

sub create_mutation_spectrum_op {
    my $self = shift;
    my $workflow = shift;
    my $sequencing_type = shift;
    my $iterator = shift;

    my $create_mutation_spectrum_op = Genome::WorkflowBuilder::Command->create(
        name => "Creating mutation spectrum results for $sequencing_type snvs using create-mutation-spectrum$iterator",
        command => 'Genome::Model::ClinSeq::Command::CreateMutationSpectrum',
    );
    $workflow->add_operation($create_mutation_spectrum_op);
    $workflow->connect_input(
        input_property       => 'build',
        destination          => $create_mutation_spectrum_op,
        destination_property => 'clinseq_build',
    );
    $workflow->connect_input(
        input_property       => "${sequencing_type}_build",
        destination          => $create_mutation_spectrum_op,
        destination_property => 'somvar_build',
    );
    $workflow->connect_input(
        input_property       => "${sequencing_type}_mutation_spectrum_outdir$iterator",
        destination          => $create_mutation_spectrum_op,
        destination_property => 'outdir',
    );
    $workflow->connect_input(
        input_property       => "${sequencing_type}_mutation_spectrum_datatype$iterator",
        destination          => $create_mutation_spectrum_op,
        destination_property => 'datatype',
    );
    $workflow->connect_input(
        input_property       => "sireport_min_bq$iterator",
        destination          => $create_mutation_spectrum_op,
        destination_property => 'min_base_quality',
    );
    $workflow->connect_input(
        input_property       => "sireport_min_mq$iterator",
        destination          => $create_mutation_spectrum_op,
        destination_property => 'min_quality_score',
    );
    $workflow->connect_output(
        output_property => "${sequencing_type}_mutation_spectrum_result$iterator",
        source          => $create_mutation_spectrum_op,
        source_property => 'result',
    );

    return $create_mutation_spectrum_op;
}

sub sciclone_op {
    my $self = shift;
    my $workflow = shift;
    my $iterator = shift;

    my $sciclone_op = Genome::WorkflowBuilder::Command->create(
        name => "Run clonality analysis and produce clonality plots using SciClone $iterator",
        command => 'Genome::Model::ClinSeq::Command::GenerateSciclonePlots',
    );
    $workflow->add_operation($sciclone_op);
    $workflow->connect_input(
        input_property       => "sciclone_dir$iterator",
        destination          => $sciclone_op,
        destination_property => 'outdir',
    );
    $workflow->connect_input(
        input_property       => 'build',
        destination          => $sciclone_op,
        destination_property => 'clinseq_build',
    );
    $workflow->connect_input(
        input_property       => 'sireport_min_coverage',
        destination          => $sciclone_op,
        destination_property => 'min_coverage',
    );
    for my $property (qw(min_mq min_bq)) {
        $workflow->connect_input(
            input_property       => "sireport_$property$iterator",
            destination          => $sciclone_op,
            destination_property => $property,
        );
    }
    $workflow->connect_output(
        output_property => "sciclone_result$iterator",
        source          => $sciclone_op,
        source_property => 'result',
    );

    return $sciclone_op;
}

sub identify_loh_op {
    my $self = shift;
    my $workflow = shift;

    my $identify_loh_op = Genome::WorkflowBuilder::Command->create(
        name => 'Identify regions of LOH and create plots',
        command => 'Genome::Model::ClinSeq::Command::IdentifyLoh',
    );
    $workflow->add_operation($identify_loh_op);
    $workflow->connect_input(
        input_property       => 'build',
        destination          => $identify_loh_op,
        destination_property => 'clinseq_build',
    );
    $workflow->connect_input(
        input_property       => 'loh_output_dir',
        destination          => $identify_loh_op,
        destination_property => 'outdir',
    );
    $workflow->connect_input(
        input_property       => 'bam_readcount_version',
        destination          => $identify_loh_op,
        destination_property => 'bamrc_version',
    );

    if ($self->name =~ /^apipe\-test/) {
        $workflow->connect_input(
            input_property       => 'test',
            destination          => $identify_loh_op,
            destination_property => 'test',
        );
    }
    $workflow->connect_output(
        output_property => 'loh_result',
        source          => $identify_loh_op,
        source_property => 'result',
    );

    return $identify_loh_op;
}

sub _infer_candidate_subjects_from_input_models {
    my $self = shift;
    my %subjects;
    my @input_models = ($self->wgs_model, $self->exome_model, $self->tumor_rnaseq_model, $self->normal_rnaseq_model);
    foreach my $input_model (@input_models) {
        next unless $input_model;
        my $patient;
        if ($input_model->subject->isa("Genome::Individual")) {
            $patient = $input_model->subject;
        }
        else {
            $patient = $input_model->subject->individual;
        }
        $subjects{$patient->id} = $patient;

        if ($input_model->can("tumor_model")) {
            $subjects{$input_model->tumor_model->subject->individual->id} =
                $input_model->tumor_model->subject->individual;
        }
        if ($input_model->can("normal_model")) {
            $subjects{$input_model->normal_model->subject->individual->id} =
                $input_model->normal_model->subject->individual;
        }
    }
    my @subjects = sort {$a->id cmp $b->id} values %subjects;
    return @subjects;
}

sub _infer_annotations_from_input_models {
    my $self = shift;
    my %annotations;
    my @input_models = ($self->wgs_model, $self->exome_model, $self->tumor_rnaseq_model, $self->normal_rnaseq_model);
    push(@input_models, $self->wgs_model->normal_model)   if $self->wgs_model;
    push(@input_models, $self->wgs_model->tumor_model)    if $self->wgs_model;
    push(@input_models, $self->exome_model->normal_model) if $self->exome_model;
    push(@input_models, $self->exome_model->tumor_model)  if $self->exome_model;

    foreach my $input_model (@input_models) {
        next unless $input_model;
        if ($input_model->can("annotation_build")) {
            my $annotation_build = $input_model->annotation_build;
            $annotations{$annotation_build->id} = $annotation_build;
        }
        if ($input_model->can("annotation_reference_build")) {
            my $annotation_build = $input_model->annotation_reference_build;
            $annotations{$annotation_build->id} = $annotation_build;
        }
    }
    my @annotations = sort {$a->id cmp $b->id} values %annotations;
    return @annotations;
}

sub _infer_references_from_input_models {
    my $self = shift;
    my %references;
    my @input_models = ($self->tumor_rnaseq_model, $self->normal_rnaseq_model);
    push(@input_models, $self->wgs_model->normal_model)   if $self->wgs_model;
    push(@input_models, $self->wgs_model->tumor_model)    if $self->wgs_model;
    push(@input_models, $self->exome_model->normal_model) if $self->exome_model;
    push(@input_models, $self->exome_model->tumor_model)  if $self->exome_model;

    foreach my $input_model (@input_models) {
        next unless $input_model;
        if ($input_model->can("reference_sequence_build")) {
            my $reference_build = $input_model->reference_sequence_build;
            $references{$reference_build->id} = $reference_build;
        }
    }
    my @references = sort {$a->id cmp $b->id} values %references;
    return @references;
}

sub _resolve_resource_requirements_for_build {
    #Set LSF parameters for the ClinSeq job
    my $lsf_resource_string = "-R 'rusage[tmp=10000:mem=4000]' -M 4000000";
    return ($lsf_resource_string);
}

# This is implemented here until refactoring is done on the model/build API
# It ensures that when the CI server compares two clinseq builds it ignores the correct files.
# Keep it in sync with the diff conditions in ClinSeq.t.
sub files_ignored_by_build_diff {
    return qw(
        build.xml
        reports/Build_Initialized/report.xml
        reports/Build_Succeeded/report.xml
        logs/.*
        .*.R$
        .*.Robject$
        .*.pdf$
        .*.jpg$
        .*.jpeg$
        .*.png$
        .*.svg$
        .*.xls$
        .*.bam$
        .*/summary/rc_summary.stderr$
        .*._COSMIC.svg$
        .*.clustered.data.tsv$
        .*QC_Report.tsv$
        .*.DumpIgvXml.log.txt$
        .*/mutation_diagrams/cosmic.mutation-diagram.stderr$
        .*/mutation_diagrams/somatic.mutation-diagram.stderr$
        .*.sequence-context.stderr$
        .*.annotate.stderr$
        .*/mutation-spectrum/exome/summarize_mutation_spectrum/summarize-mutation-spectrum.stderr$
        .*/mutation-spectrum/wgs/summarize_mutation_spectrum/summarize-mutation-spectrum.stderr$
        .*LIMS_Sample_Sequence_QC_library.tsv$
        .*.R.stderr$
        .*/snv_indel_report/subjects_legend.txt$
        .*sciclone.*.clusters.txt$
        .*mutation_spectrum.input.tsv$
        .*/kinase_only.tsv$
        .*/all_interactions.tsv$
        .*/expert_antineoplastic.tsv$
    );
}

# Below are some methods to retrieve files from a build
# Ideally they would be tied to creation of these file paths, but they currently
# throw an exception if the files don't exist. Consider another approach...
sub patient_dir {
    my ($self, $build) = @_;
    unless ($build->common_name) {
        die $self->error_message("Common name is not defined.");
    }
    my $patient_dir = $build->data_directory . "/" . $build->common_name;
    unless (-d $patient_dir) {
        die $self->error_message("ClinSeq patient directory not found. Expected: $patient_dir");
    }
    return $patient_dir;
}

sub snv_variant_source_file {
    my ($self, $build, $data_type,) = @_;

    my $patient_dir = $self->patient_dir($build);
    my $source      = $patient_dir . "/variant_source_callers";
    my $exception;
    if (-d $source) {
        my $dir = $source . "/$data_type";
        if (-d $dir) {
            my $file = $dir . "/snv_sources.tsv";
            if (-e $file) {
                return $file;
            }
            else {
                $exception = $self->error_message("Expected $file inside $dir and it did not exist.");
            }
        }
        else {
            $exception = $self->error_message("$data_type sub-directory not found in $source.");
        }
    }
    else {
        $exception = $self->error_message("$source directory not found");
    }
    die $exception;
}

sub copy_fusion_files {
    my ($self, $build) = @_;
    my $rnaseq_build_dir             = $build->tumor_rnaseq_build->data_directory;
    my $tumor_unfiltered_fusion_file = $rnaseq_build_dir
        . '/fusions/Genome_Model_RnaSeq_DetectFusionsResult_Chimerascan_VariableReadLength_Result/chimeras.bedpe';
    my $tumor_filtered_fusion_file           = $rnaseq_build_dir . '/fusions/filtered_chimeras.bedpe';
    my $tumor_filtered_annotated_fusion_file = $rnaseq_build_dir . '/fusions/filtered_chimeras.catanno.bedpe';
    my $clinseq_tumor_unfiltered_fusion_file = $self->patient_dir($build) . '/rnaseq/tumor/fusions/chimeras.bedpe';
    my $clinseq_tumor_filtered_fusion_file =
        $self->patient_dir($build) . '/rnaseq/tumor/fusions/filtered_chimeras.bedpe';
    my $clinseq_tumor_filtered_annotated_fusion_file =
        $self->patient_dir($build) . '/rnaseq/tumor/fusions/filtered_chimeras.catanno.bedpe';

    if (-e $tumor_unfiltered_fusion_file) {
        unless (Genome::Sys->copy_file($tumor_unfiltered_fusion_file, $clinseq_tumor_unfiltered_fusion_file)) {
            die "unable to copy $tumor_unfiltered_fusion_file";
        }
    }
    if (-e $tumor_filtered_fusion_file) {
        unless (Genome::Sys->copy_file($tumor_filtered_fusion_file, $clinseq_tumor_filtered_fusion_file)) {
            die "unable to copy $tumor_filtered_fusion_file";
        }
    }
    if (-e $tumor_filtered_annotated_fusion_file) {
        unless (
            Genome::Sys->copy_file(
                $tumor_filtered_annotated_fusion_file, $clinseq_tumor_filtered_annotated_fusion_file
            )
            )
        {
            die "unable to copy $tumor_filtered_annotated_fusion_file";
        }
    }
}

sub clonality_dir {
    my ($self, $build) = @_;

    my $patient_dir   = $self->patient_dir($build);
    my $clonality_dir = $patient_dir . "/clonality";

    unless (-d $clonality_dir) {
        die $self->error_message("Clonality directory does not exist. Expected: $clonality_dir");
    }
    return $clonality_dir;
}

sub varscan_formatted_readcount_file {
    my ($self, $build) = @_;
    my $clonality_dir  = $self->clonality_dir($build);
    my $readcount_file = $clonality_dir . "/allsnvs.hq.novel.tier123.v2.bed.adapted.readcounts.varscan";
    unless (-e $readcount_file) {
        die $self->error_message("Unable to find varscan formatted readcount file. Expected: $readcount_file");
    }
    return $readcount_file;
}

sub cnaseq_hmm_file {
    my ($self, $build) = @_;
    my $clonality_dir = $self->clonality_dir($build);
    my $hmm_file      = $clonality_dir . "/cnaseq.cnvhmm";
    unless (-e $hmm_file) {
        die $self->error_message("Unable to find cnaseq hmm file. Expected: $hmm_file");
    }
    return $hmm_file;
}

sub has_microarray_build {
    my $self = shift;
    my $base_model;
    if ($self->exome_model) {
        $base_model = $self->exome_model;
    }
    elsif ($self->wgs_model) {
        $base_model = $self->wgs_model;
    }
    else {
        return 0;
    }
    if ($base_model->tumor_model->genotype_microarray && $base_model->normal_model->genotype_microarray) {
        if (   $base_model->tumor_model->genotype_microarray->last_succeeded_build
            && $base_model->normal_model->genotype_microarray->last_succeeded_build)
        {
            return 1;
        }
        else {
            return 0;
        }
    }
    else {
        return 0;
    }
}

sub _get_docm_variants_file {
    my $self                 = shift;
    my $cancer_annotation_db = shift;
    my $cad_data_directory   = $cancer_annotation_db->data_directory;
    my $docm_variants_file   = $cad_data_directory . "/DOCM/DOCM_v0.1.tsv";
    if (-e $docm_variants_file) {
        return $docm_variants_file;
    }
    else {
        $self->status_message(
            "Unable to find DOCM variants file in cancer annotation db directory $docm_variants_file");
        return 0;
    }
}

1;
