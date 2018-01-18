package Genome::Model::Build::SomaticValidation;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::SomaticValidation {
    is => [
        'Genome::Model::Build',
        'Genome::Model::Build::RunsDV2',
        'Genome::Model::Build::HasFeatureLists',
        'Genome::Model::Build::SomaticInterface',
    ],
    has_optional => [
        reference_sequence_build => {
            is => 'Genome::Model::Build::ReferenceSequence', via => 'inputs', to => 'value', where => [name => 'reference_sequence_build'],
        },
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            via => 'inputs', to => 'value', where => [ name => 'annotation_build' ],
            is_mutable => 1,
        },
        previously_discovered_variations_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            via => 'inputs', to => 'value', where => [ name => 'previously_discovered_variations_build' ],
            is_mutable => 1,
            doc => 'build of variants to screen out from consideration (such as from dbSNP)',
        },
        snv_variant_list => {
            is => 'Genome::SoftwareResult',
            via => 'inputs', to => 'value', where => [ name => 'snv_variant_list' ],
            is_mutable => 1,
        },
        indel_variant_list => {
            is => 'Genome::SoftwareResult',
            via => 'inputs', to => 'value', where => [ name => 'indel_variant_list' ],
            is_mutable => 1,
        },
        sv_variant_list => {
            is => 'Genome::SoftwareResult',
            via => 'inputs', to => 'value', where => [ name => 'sv_variant_list' ],
            is_mutable => 1,
        },
        alignment_strategy => {
            is => 'Text',
            via => 'model',
        },
        snv_detection_strategy => {
            is => 'Text',
            via => 'model',
        },
        sv_detection_strategy => {
            is => 'Text',
            via => 'model',
        },
        indel_detection_strategy => {
            is => 'Text',
            via => 'model',
        },
        cnv_detection_strategy => {
            is => 'Text',
            via => 'model',
        },
        target_region_set_name => {
            is => 'Text',
            via => 'target_region_set',
            to => 'name',
        },
        target_region_set => {
            is => 'Genome::FeatureList',
            via => 'inputs', to => 'value', where => [ name => 'target_region_set' ],
            is_mutable => 1,
        },
        region_of_interest_set => {
            is => 'Genome::FeatureList',
            via => 'inputs', to => 'value', where => [ name => 'region_of_interest_set' ],
            is_mutable => 1,
        },
        design_set => {
            is => 'Genome::FeatureList',
            via => 'inputs', to => 'value', where => [ name => 'design_set' ],
            is_mutable => 1,
        },
        tumor_sample => {
            is => 'Genome::Sample',
            via => 'inputs', to => 'value', where => [ name => 'tumor_sample' ],
            is_mutable => 1,
        },
        normal_sample => {
            is => 'Genome::Sample',
            via => 'inputs', to => 'value', where => [ name => 'normal_sample' ],
            is_mutable => 1,
        },

        merged_alignment_result => {
            is => 'Genome::InstrumentData::AlignedBamResult::Merged',
            via => 'result_users',
            to => 'software_result',
            where => [label => 'merged_alignment'],
        },
        control_merged_alignment_result => {
            is => 'Genome::InstrumentData::AlignedBamResult::Merged',
            via => 'result_users',
            to => 'software_result',
            where => [label => 'control_merged_alignment'],
        },

        coverage_stats_result => {
            is => 'Genome::InstrumentData::AlignmentResult::Merged::CoverageStats',
            via => 'result_users',
            to => 'software_result',
            where => [label => 'coverage_stats_tumor'],
        },
        control_coverage_stats_result  => {
            is => 'Genome::InstrumentData::AlignmentResult::Merged::CoverageStats',
            via => 'result_users',
            to => 'software_result',
            where => [label => 'coverage_stats_normal'],
        },

        loh_version => {
            via => 'model',
        },
        loh_snv_detection_strategy => {
            via => 'model',
        },
        tiering_version => {
            via => 'model',
        },

        map { $_ => { via => 'model' } } ('run_snv_validation', 'run_sv_validation', 'run_indel_validation'),
    ],
};

sub post_allocation_initialization {
    my $self = shift;

    my @result_subfolders;
    for my $subdir ('alignments', 'variants', 'coverage') {
        push @result_subfolders, $self->data_directory."/".$subdir;
    }

    for my $subdir (@result_subfolders){
        Genome::Sys->create_directory($subdir) unless -d $subdir;
    }

    for my $variant_type ('snv', 'indel', 'sv') {
        my $property = $variant_type . '_variant_list';
        my $variant_list = $self->$property;

        if($variant_list) {
            $variant_list->add_user(label => 'uses', user => $self);
        }
    }

    return 1;
}

sub tumor_bam {
    my $self = shift;

    my $result = $self->merged_alignment_result;
    return unless $result;
    return $result->bam_file;
}

sub normal_bam {
    my $self = shift;

    my $result = $self->control_merged_alignment_result;
    return unless $result;
    return $result->bam_file;
}

sub workflow_name {
    my $self = shift;
    return $self->build_id . ' Somatic Variant Validation Pipeline';
}

sub calculate_estimated_kb_usage {
    my $self = shift;

    # FIXME find out how much we probably really need
    return 15_728_640;
}

sub regex_files_for_diff {
    return qw(
        ^alignments/normal/[[:xdigit:]]+\.(bam.*)$
        ^alignments/tumor/[[:xdigit:]]+\.(bam.*)$
        sv/alignments/tumor/[[:xdigit:]]+\.(bam.*)$
        sv/alignments/normal/[[:xdigit:]]+\.(bam.*)$
        validation/large_indel/alignments/tumor/[[:xdigit:]]+\.(bam.*)$
        validation/large_indel/alignments/normal/[[:xdigit:]]+\.(bam.*)$
        coverage/(tumor|normal)/wingspan_(\d+)/[[:xdigit:]]+_(\w+)_STATS.t(sv|xt)
        coverage/(tumor|normal)/[[:xdigit:]]+-(\w+)-wingspan_(\d+)-alignment_summary.tsv
        coverage/(tumor|normal)/[[:xdigit:]]+-(\w+)-wingspan_(\d+)-alignment_summary-v2.tsv
    );
}

sub files_ignored_by_diff {
    return qw(
        reports/Build_Initialized/report.xml
        reports/Build_Succeeded/report.xml
        variants/dispatcher.cmd
        \.vcf$
        \.vcf.idx$
        workflow\.xml$
        build\.xml$
        \d+\.out$
        \d+\.err$
        \.png$
        readcounts$
        variants/sv/breakdancer
        variants/sv/squaredancer
        variants/sv/union
        variants/svs.hq
        svs\.merge\.index$
        sv/assembly_output\.csv\.index$
        sv/assembly_input
        sv/assembly_output.csv
        cnv_graph\.pdf$
        pindel\.config$
        variants/cnv/plot-cnv-[^/]+/[[:xdigit:]]+\.bam\.cnvseg$
        variants/cnv/plot-cnv-[^/]+/[[:xdigit:]]+.bam.bamtocn$
        variants/indel/pindel-[^/]+/pindel-somatic-calls-[^/]+/pindel-read-support-[^/]+/indels.hq.read_support.bed$
        alignments/.*\.log$
        alignments/.*\.metrics$
        alignments/.*\.bam\.md5$
        alignments/.*\.bam$
        alignments/.*\.bam\.bai$
        control_variants_for_loh/dispatcher.cmd
        validation/review/newcalls.xml
        validation/small_indel/indel_files_to_validate
        validation/large_indel/tumor.csv
        validation/large_indel/normal.csv
        variants/(.*)\.tbi$
        control_variants_for_loh/(.*)\.tbi$
        variants/snv/strelka.*/output/.*
        variants/indel/strelka.*/output/.*
    );
}
sub dirs_ignored_by_diff {
    return qw(
        logs/
        reports/
        variants/\d+/
        variants/sv/breakdancer
        variants/sv/squaredancer
        variants/sv/union
        validation/small_indel/realigned_bams
        /validation/small_indel
    );
}

sub validate_for_start_methods {
    my $self = shift;
    my @methods = $self->SUPER::validate_for_start_methods;
    push @methods, qw(
        _validate_required_for_start_properties
    );
    return @methods;
}

sub _validate_required_for_start_properties {
    my $model_method = $_[0]->model->class . '::_validate_required_for_start_properties';
    return (\&$model_method)->(@_);
}

sub reference_being_replaced_for_input {
    my ($self, $input) = @_;

    if($input->name eq "target_region_set" or $input->name eq "design_set"){
        return 1;
    }

    if($input->name eq "region_of_interest_set"){
        my $rsb = $self->reference_sequence_build;
        my $roi_reference = $input->value->reference;
        unless ($roi_reference) {
            return;
        }

        if ($roi_reference and !$rsb->is_compatible_with($roi_reference)) {
            my $converter_exists =  Genome::Model::Build::ReferenceSequence::Converter->exists_for_references(
                $roi_reference, $rsb,
            );

            if ($converter_exists) {
                return 1;
            }
        }
    }

    if($input->name eq 'instrument_data') {
        return ($self->processing_profile->alignment_strategy !~ 'imported'); #as long as it's not "imported", we'll re-align
    }

    return;
}

sub whole_rmdup_bam_file {
    my $self = shift;
    return $self->tumor_bam;
}

sub get_target_region_feature_list {
    my $self = shift;

    if (defined($self->target_region_set_name)) {
        return Genome::FeatureList->get(
            name => $self->target_region_set_name);
    } else {
        return;
    }
}

sub get_feature_list_from_reference {
    my ($self, $feature_list_accessor) = @_;
    return $self->reference_sequence_build->get_feature_list($feature_list_accessor);
}

sub indels_effects_file {
    my $self = shift;
    my $tier = shift;
    return File::Spec->join($self->data_directory, 'effects', "indels.hq.$tier.bed");
}

sub has_microarray_build {
    my $self = shift;
    return 0;
}

sub has_snvs_annotated_variants_vcf_file {
    my $self = shift;
    return 0;
}

sub build_for_normal_bam {
    my $self = shift;
    return $self;
}

sub build_for_tumor_bam {
    my $self = shift;
    return $self;
}

sub _disk_usage_result_subclass_names {
    my $self = shift;

    my @classes = ($self->_dv2_result_subclass_names, qw(
        Genome::InstrumentData::AlignmentResult::Bwa
        Genome::InstrumentData::AlignmentResult::Bwamem
        Genome::InstrumentData::AlignmentResult::BwamemStream
        Genome::InstrumentData::AlignmentResult::Merged
        Genome::InstrumentData::AlignmentResult::Merged::BamQc
        Genome::InstrumentData::AlignmentResult::Merged::CoverageStats
        Genome::InstrumentData::AlignmentResult::Speedseq
        Genome::InstrumentData::AlignmentResult::Merged::Speedseq
        Genome::Qc::Result
    ));

    return \@classes;
}

1;
