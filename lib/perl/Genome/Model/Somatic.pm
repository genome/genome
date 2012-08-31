package Genome::Model::Somatic;

use strict;
use warnings;

use Genome;

class Genome::Model::Somatic {
    is  => 'Genome::ModelDeprecated',
    has_param => [
        only_tier_1 => {
            doc => "If set to true, the pipeline will skip ucsc annotation and produce only tier 1 snps",
        },
        min_mapping_quality => {
            doc => "minimum average mapping quality threshold for high confidence call",
        },
        min_somatic_quality => {
            doc => "minimum somatic quality threshold for high confidence call",
        },
        skip_sv => {
            doc => "If set to true, the pipeline will skip structural variation detection",
        },
        sv_detector_version => {
            doc => "Version of the SV detector to use. Reference the Genome::Model::Tools::DetectVariants::SomaticBreakdancer module for currently available versions.",
        },
        sv_detector_params => {
            doc => "Parameters to pass to the SV detector.  For breakdancer, separate params for bam2cfg & BreakDancerMax with a colon. If no parameters are desired, just provide ':'.",
        },
        bam_window_version => {
            doc => "Version to use for bam-window in the copy number variation step.",
        },
        bam_window_params => {
            doc => "Parameters to pass to bam-window in the copy number variation step.",
        },
        sniper_version => {
            doc => "Version to use for bam-somaticsniper for detecting snps and indels.",
        },
        sniper_params => {
            doc => "Parameters to pass to bam-somaticsniper for detecting snps and indels",
        },
        snv_detector_name => {
            doc => "The name of the variant detector to use for snv detection",
        },
        snv_detector_version => {
            doc => "The version of the variant detector to use for snv detection",
        },
        snv_detector_params => {
            doc => "The params to pass to the variant detector to use for snv detection",
        },
        indel_detector_name => {
            doc => "The name of the variant detector to use for indel detection",
        },
        indel_detector_version => {
            doc => "The version of the variant detector to use for indel detection",
        },
        indel_detector_params => {
            doc => "The params to pass to the variant detector to use for indel detection",
        },
        bam_readcount_version => {
            doc => "Version to use for bam-readcount in the high confidence step.",
        },
        bam_readcount_params=> {
            doc => "Parameters to pass to bam-readcount in the high confidence step",
        },
        require_dbsnp_allele_match => {
            doc => "If set to true, the pipeline will require the allele to match during Lookup Variants"  
        },
        transcript_variant_annotator_version => {
            doc => 'Version of the "annotate transcript-variants" tool to run during the annotation step',
            is_optional => 1,
            default_value => Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version,
            valid_values => [ 0,1,2,3],#Genome::Model::Tools::Annotate::TranscriptVariants->available_versions ],
        },
    ],
    has_optional => [
        tumor_model_links => {
            is => 'Genome::Model::Link',
            reverse_as => 'to_model',
            where => [ role => 'tumor'],
            is_many => 1,
        },
        tumor_model => {
            is => 'Genome::Model',
            via => 'tumor_model_links',
            to => 'from_model',
        },
        tumor_model_id => {
            is => 'Integer',
            via => 'tumor_model',
            to => 'id',
        },
        normal_model_links => {
            is => 'Genome::Model::Link',
            reverse_as => 'to_model',
            where => [ role => 'normal'],
            is_many => 1,
        },
        normal_model => {
            is => 'Genome::Model',
            via => 'normal_model_links',
            to => 'from_model',
        },
        normal_model_id => {
            is => 'Integer',
            via => 'normal_model',
            to => 'id',
        },
    ],
};

sub create {
    my $class = shift;
    my %params = @_;

    my $tumor_model_id = delete $params{tumor_model_id};
    my $tumor_model = delete $params{tumor_model};
    my $normal_model_id = delete $params{normal_model_id};
    my $normal_model = delete $params{normal_model};

    unless($tumor_model) {
        $tumor_model = Genome::Model->get($tumor_model_id);

        unless($tumor_model) {
            $class->error_message('Could not find tumor model.' );
            return;
        }
    }

    unless($normal_model) {
        $normal_model = Genome::Model->get($normal_model_id);

        unless($normal_model) {
            $class->error_message('Could not find normal model.');
            return;
        }
    }

    my $tumor_subject = $tumor_model->subject;
    my $normal_subject = $normal_model->subject;

    if($tumor_subject->can('source') and $normal_subject->can('source')) {
        my $tumor_source = $tumor_subject->source;
        my $normal_source = $normal_subject->source;

        if($tumor_source eq $normal_source) {
            my $subject = $tumor_source;

            #Set up other parameters for call to parent execute()
            $params{subject_id} = $subject->id;
            $params{subject_class_name} = $subject->class;
        } else {
            $class->error_message('Tumor and normal samples are not from same source!');
            return;
        }
    } else {
        $class->error_message('Unexpected subject for tumor or normal model!');
        return;
    }

    my $self = $class->SUPER::create(%params);

    $self->add_from_model(from_model => $normal_model, role => 'normal');
    $self->add_from_model(from_model => $tumor_model, role => 'tumor');

    return $self;
}

sub get_all_objects {
    my $self = shift;

    my @objects = $self->SUPER::get_all_objects(@_);
    my @validations = Genome::Model::VariantValidation->get(model_id=>$self->id);
    push @objects, @validations;
    return @objects;
}

sub _initialize_build {
    my($self,$build) = @_;

    # Check that the annotator version param is sane before doing the build
    my $annotator_version;
    my $worked = eval {
        $annotator_version = $self->transcript_variant_annotator_version;
        # When all processing profiles have a param for this, remove this unless block so
        # they'll fail if it's missing
        unless (defined $annotator_version) {
            $annotator_version = Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version;
        }

        my %available_versions = map { $_ => 1 } Genome::Model::Tools::Annotate::TranscriptVariants->available_versions;
        unless ($available_versions{$annotator_version}) {
            die "Requested annotator version ($annotator_version) is not in the list of available versions: "
                . join(', ',keys(%available_versions));
        }
        1;
    };
    unless ($worked) {
        $self->error_message("Could not determine which version of the Transcript Variants annotator to use: $@");
        return;
    }
    return 1;
}

sub _resolve_workflow_for_build {
    my $self = shift;
    my $build = shift;

    my $operation = Workflow::Operation->create_from_xml(__FILE__ . '.xml');

    my $log_directory = $build->log_directory;
    $operation->log_dir($log_directory);
    $operation->name($build->workflow_name);

    return $operation;
}

sub map_workflow_inputs {
    my $self = shift;
    my $build = shift;

    my @inputs = ();

    my $tumor_build = $build->tumor_build;
    my $normal_build = $build->normal_build;

    unless ($tumor_build) {
        $self->error_message("Failed to get a tumor_build associated with this somatic capture build!");
        die $self->error_message;
    }

    unless ($normal_build) {
        $self->error_message("Failed to get a normal_build associated with this somatic capture build!");
        die $self->error_message;
    }

    my $data_directory = $build->data_directory;
    unless ($data_directory) {
        $self->error_message("Failed to get a data_directory for this build!");
        die $self->error_message;
    }

    my $tumor_bam = $tumor_build->whole_rmdup_bam_file;
    unless (-e $tumor_bam) {
        $self->error_message("Tumor bam file $tumor_bam does not exist!");
        die $self->error_message;
    }

    my $normal_bam = $normal_build->whole_rmdup_bam_file;
    unless (-e $normal_bam) {
        $self->error_message("Normal bam file $normal_bam does not exist!");
        die $self->error_message;
    }

    # Get the snp file from the tumor and normal models
    my $tumor_snp_file = $tumor_build->snv_file;
    unless (-e $tumor_snp_file) {
        $self->error_message("Tumor snp file $tumor_snp_file does not exist!");
        die $self->error_message;
    }
    my $normal_snp_file = $normal_build->snv_file;
    unless (-e $normal_snp_file) {
        $self->error_message("Normal snp file $normal_snp_file does not exist!");
        die $self->error_message;
    }

    push @inputs,
        build_id => $build->id,
        normal_bam_file => $normal_bam,
        tumor_bam_file => $tumor_bam,
        normal_filtered_snp_file => $normal_snp_file,
        tumor_filtered_snp_file => $tumor_snp_file,
        data_directory => $data_directory;


    my %default_filenames = $self->default_filenames;
    for my $param (keys %default_filenames) {
        # set a default param if one has not been specified
        my $default_filename = $default_filenames{$param};
        push @inputs,
            $param => join('/', $data_directory, $default_filename);
    }

    #These filenames don't take a directory prepended
    push @inputs,
        snp_base_name => 'snp_output.csv',
        indel_base_name => 'indel_output.csv',
        sv_base_name => 'sv_output.csv';

    # Set (hardcoded) defaults for tools that have defaults that do not agree with somatic pipeline
    push @inputs,
        skip_if_output_present => 0,
        imported_bams => 0,
        lookup_variants_report_mode => "novel-only",
        lookup_variants_filter_out_submitters => "SNP500CANCER,OMIMSNP,CANCER-GENOME,CGAP-GAI,LCEISEN,ICRCG,DEVINE_LAB",
        annotate_no_headers => 1,
        transcript_annotation_filter => "top",
        annotation_reference_transcripts => 'NCBI-human.combined-annotation/54_36p_v2',
        #FIXME Get the reference from the constituent builds instead of hard-coding
        reference_fasta => (Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa');

    # Set values from processing profile parameters
    push @inputs,
        only_tier_1 => (defined($self->only_tier_1)? $self->only_tier_1 : 0),
        skip_sv => (defined($self->skip_sv)? $self->skip_sv : 0),
        breakdancer_version => (defined($self->sv_detector_version)? $self->sv_detector_version : ''),
        sv_params => (defined($self->sv_detector_params)? $self->sv_detector_params : ':'),
        transcript_variant_annotator_version => (defined($self->transcript_variant_annotator_version)? $self->transcript_variant_annotator_version : ':'),
        min_mapping_quality => (defined($self->min_mapping_quality)? $self->min_mapping_quality : 70),
        min_somatic_quality => (defined($self->min_somatic_quality)? $self->min_somatic_quality : 40),
        require_dbsnp_allele_match => (defined($self->require_dbsnp_allele_match)? $self->require_dbsnp_allele_match : 1),
        sniper_params => (defined($self->sniper_params)? $self->sniper_params : ''),
        sniper_version => (defined($self->sniper_version)? $self->sniper_version : ''),
        bam_window_params => (defined($self->bam_window_params)? $self->bam_window_params : ''),
        bam_window_version => (defined($self->bam_window_version)? $self->bam_window_version : ''),
        bam_readcount_params => (defined($self->bam_readcount_params)? $self->bam_readcount_params : ''),
        bam_readcount_version => (defined($self->bam_readcount_version)? $self->bam_readcount_version : '');

    return @inputs;
}

sub default_filenames{
    my $self = shift;

    my %default_filenames = (
        breakdancer_working_directory       => 'breakdancer/',
        sniper_working_directory            => 'sniper/',
        snp_filter_output                   => 'sfo_snp_filtered.csv',
        loh_fail_output_file                => 'loh.csv',
        loh_output_file                     => 'noloh.csv',
        filter_ceu_yri_output               => 'ceu_yri_filtered.csv',
        adaptor_output_snp                  => 'adv_adapted_snp.csv',
        dbsnp_output                        => 'dbsnp_filtered.csv',
        annotate_output_snp                 => 'anv_annotated_snp.csv',
        ucsc_output_snp                     => 'uca_ucsc_annotated_snp.csv',
        ucsc_unannotated_output_snp         => 'ucu_ucsc_unannotated_snp.csv',
        tier_1_snp_file                     => 't1v_tier1_snp.csv',
        tier_2_snp_file                     => 't2v_tier2_snp.csv',
        tier_3_snp_file                     => 't3v_tier3_snp.csv',
        tier_4_snp_file                     => 't4v_tier4_snp.csv',
        tier_1_snp_high_confidence_file     => 'hc1_tier1_snp_high_confidence.csv',
        tier_2_snp_high_confidence_file     => 'hc2_tier2_snp_high_confidence.csv',
        tier_3_snp_high_confidence_file     => 'hc3_tier3_snp_high_confidence.csv',
        tier_4_snp_high_confidence_file     => 'hc4_tier4_snp_high_confidence.csv',
        upload_variants_snp_1_output        => 'uv1_uploaded_tier1_snp.csv',
        upload_variants_snp_2_output        => 'uv2_uploaded_tier2_snp.csv',
        indel_lib_filter_multi_output       => 'iml_indel_multi_lib_filtered.csv',
        indel_lib_filter_single_output      => 'isl_indel_single_lib_filtered.csv',
        adaptor_output_indel                => 'adi_adaptor_output_indel.csv',
        annotate_output_indel               => 'ani_annotated_indel.csv',
        ucsc_output_indel                   => 'uci_ucsc_annotated_indel.csv',
        ucsc_unannotated_output_indel       => 'ucn_ucsc_unannotated_indel.csv',
        tier_1_indel_file                   => 't1i_tier1_indel.csv',
        tier_2_indel_file                   => 't2i_tier2_indel.csv',
        tier_3_indel_file                   => 't3i_tier3_indel.csv',
        tier_4_indel_file                   => 't4i_tier4_indel.csv',
        upload_variants_indel_1_output      => 'ui1_uploaded_tier1_indel.csv',
        upload_variants_indel_2_output      => 'ui2_uploaded_tier2_indel.csv',
        copy_number_output                  => 'cno_copy_number.csv',
        circos_graph                        => 'circos_graph.png',
        variant_report_output               => 'cancer_report.html', 
        file_summary_report_output          => 'file_summary_report.html', 
        indel_lib_filter_preferred_output   => 'NULL',
    );

    return %default_filenames;
}

1;
