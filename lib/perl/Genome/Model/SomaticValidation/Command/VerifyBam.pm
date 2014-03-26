package Genome::Model::SomaticValidation::Command::VerifyBam;

use strict;
use warnings;
use Genome;

class Genome::Model::SomaticValidation::Command::VerifyBam {
    is => ['Command::V2'],
    has_transient_optional_output => [
        result => {
            is => 'Genome::InstrumentData::VerifyBamIdResult',
        },
    ],
};

sub shortcut {
    my $self = shift;

    unless ($self->should_run) {
        return 1;
    }

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

    unless($self->should_run) {
        return 1;
    }

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

sub should_run {
    my $self = shift;

    unless ($self->SUPER::should_run) {
        return 0;
    }

    unless (defined $self->sample_for_mode->default_genotype_data) {
        $self->debug_message('No default genotype data for sample '.$self->sample_for_mode->__display_name__.' Skipping VerifyBamId');
        return 0;
    }
}

sub params_for_result {
    my $self = shift;
    my $genotype_vcf = Genome::InstrumentData::GenotypeVcf->get_or_create(
        sample => $self->sample_for_mode,
        known_sites_build => $self->build->previously_discovered_variations_build,
        filters => ["chromosome:exclude=".$self->build->previously_discovered_variations_build->reference_sequence_build->allosome_names],
    );
    my %params = (
        aligned_bam_result_id => $self->alignment_result_for_mode,
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
    $self->result($result);
    $self->add_metrics_to_build;
    return $self->SUPER::link_result_to_build($result, "verifyBamId", "verify_bam_id");
}

sub add_metrics_to_build {
    my $self = shift;
    return 1;
}

1;

