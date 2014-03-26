package Genome::Model::SomaticValidation::Command::CoverageStats;

use strict;
use warnings;
use Genome;

class Genome::Model::SomaticValidation::Command::CoverageStats {
    is => ['Command::V2'],
    has_input => [
        merged_alignment_result_id => {
            is => 'Text',
        },
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
        },
    ],
};

sub shortcut {
    my $self = shift;
    my $params = $self->params_for_result;
    my $result = Genome::InstrumentData::VerifyBamIdResult->get_with_lock(
        %{$params}
    );
    if ($result) {
        $self->debug_message("Using existing result ".$result->__display_name__);
        return $self->link_result_to_build($result);
    }
    return;
}

sub execute {
    my $self = shift;
    my $params = $self->params_for_result;
    my $result = Genome::InstrumentData::VerifyBamIdResult->get_or_create(
        %{$params}
    );
    if ($result) {
        $self->debug_message("Created result ".$result->__display_name__);
        return $self->link_result_to_build($result);
    }
    else {
        $self->error_message("Failed to create result");
        return;
    }
}

sub params_for_result {
    my $self = shift;
    my $aligned_bam_result = Genome::SoftwareResult->get($self->merged_alignment_result_id);
    my $genotype_vcf = Genome::InstrumentData::GenotypeVcf->get_or_create(
        sample => $aligned_bam_result->instrument_data->sample,
        known_sites_build => $self->build->previously_discovered_variations_build,
        filters => ["chromosome:exclude=".$self->build->previously_discovered_variations_build->reference_sequence_build->allosome_names],
    );
    my %params = (
        aligned_bam_result_id => $self->merged_alignment_result_id,
        max_depth => 1000,
        precise => 1,
        version => $self->build->verify_bam_id_version,
    );
    if (defined $self->build->target_region_set) {
        $params{on_target_list} = $self->build->target_region_set;
    }
    return \%params;
}

sub link_result_to_build {
    my $self = shift;
    my $result = shift;
    my $build = $self->build;

    my $link = join('/', $build->data_directory, 'verifyBamId', $self->merged_alignment_result_id);
    my $label = join('_', 'verify_bam_id');
    Genome::Sys->create_symlink($result->output_dir, $link);
    $result->add_user(label => $label, user => $build);
    $self->add_metrics_to_build;
    return 1;
}

sub add_metrics_to_build {
    my $self = shift;
    return 1;
}

1;

