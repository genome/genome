package Genome::Model::Build::SomaticVariation::CalcCovgResult;

use strict;
use warnings;

use Sys::Hostname;
use Genome;

class Genome::Model::Build::SomaticVariation::CalcCovgResult {
    is => 'Genome::SoftwareResult::Stageable',
    has_input => [
        somatic_variation_build => {
            is => 'Genome::Model::Build::SomaticVariation',
            doc => 'Build containing normal and tumor bams to analyze',
        },
    ],
    has_param => [
        reference_sequence => {
            is => 'Text',
            doc => "Path to reference sequence in FASTA format",
        },
        normal_min_depth => {
            is => 'Text',
            doc => "The minimum read depth to consider a Normal BAM base as covered",
        },
        tumor_min_depth => {
            is => 'Text',
            doc => "The minimum read depth to consider a Tumor BAM base as covered",
        },
        min_mapq => {
            is => 'Text',
            doc => "The minimum mapping quality of reads to consider towards read depth counts",
        },
        roi_file => {
            is => 'Text',
            doc => "Tab delimited list of ROIs [chr start stop gene_name] (See Description)",
        },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    $self->_prepare_staging_directory;
    $self->_prepare_output_directory;

    $self->status_message("CalcCovg for build ".$self->somatic_variation_build->id);

    my $rv = $self->_run_calc_covg_helper;
    unless ($rv) {
        $self->error_message("CalcCovg failed for build ".$self->somatic_variation_build->id);
        $self->delete;
        return;
    }

    my $status = "CalcCovg done";
    $self->status_message($status);

    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;

}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $hostname = hostname;

    my $user = $ENV{'USER'};
    my $base_dir = sprintf("calc_covg_result-%s-%s-%s-%s", $hostname, $user, $$, $self->id);
    my $directory = join('/', 'build_merged_alignments',$self->id,$base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    return 'info_genome_models';
}

sub output_file {
    my $self = shift;

    my $sample_name = $self->somatic_variation_build->tumor_build->subject->extraction_label;

    return $sample_name.".covg";
}

sub _run_calc_covg_helper {
    my $self = shift;

    my $sample_name = $self->somatic_variation_build->tumor_build->subject->extraction_label;
    my $normal_bam = $self->somatic_variation_build->normal_bam;
    my $tumor_bam = $self->somatic_variation_build->tumor_bam;

    my $cmd = Genome::Model::Tools::Music::Bmr::CalcCovgHelper->create (
        roi_file => $self->roi_file,
        reference_sequence => $self->reference_sequence,
        normal_tumor_bam_pair => "$sample_name\t$normal_bam\t$tumor_bam",
        output_file => $self->temp_staging_directory."/".$self->output_file,
        normal_min_depth => $self->normal_min_depth,
        tumor_min_depth => $self->tumor_min_depth,
        min_mapq => $self->min_mapq,
    );

    my $rv = $cmd->execute;

    return $rv;
}

sub _gather_params_for_get_or_create {
    my $class = shift;

    my $bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);

    my %params = $bx->params_list;
    my %is_input;
    my %is_param;
    my $class_object = $class->__meta__;
    for my $key ($class->property_names) {
        my $meta = $class_object->property_meta_for_name($key);
        if ($meta->{is_input} && exists $params{$key}) {
            $is_input{$key} = $params{$key};
        } elsif ($meta->{is_param} && exists $params{$key}) {
            $is_param{$key} = $params{$key};
        }
    }

    my $inputs_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_input);
    my $params_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_param);

    my %software_result_params = (
        params_id=>$params_bx->id,
        inputs_id=>$inputs_bx->id,
        subclass_name=>$class
    );

    return {
        software_result_params => \%software_result_params,
        subclass => $class,
        inputs=>\%is_input,
        params=>\%is_param,
    };
}

1;
