package Genome::InstrumentData::Command::RefineReads::GatkBestPractices;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Gatk::IndelRealignerResult;
use Genome::InstrumentData::Gatk::BaseRecalibratorBamResult;
require Genome::Utility::Text;

class Genome::InstrumentData::Command::RefineReads::GatkBestPractices {
    is => 'Command::V2',
    has => [
        version => { is => 'Text', },
        bam_source => { is => 'Genome::InstrumentData::AlignedBamResult', },
    ],
    has_many => [
        known_sites => { is => 'Genome::Model::Build::ImportedVariationList', },
    ],
    has_optional => [
        params => { is => 'Text', }, # not used
    ],
    has_optional_transient => [
        indel_realigner_result => { is => 'Genome::InstrumentData::Gatk::IndelRealignerResult', },
        base_recalibrator_bam_result => { is => 'Genome::InstrumentData::Gatk::BaseRecalibratorBamResult', },
    ],
};

sub shortcut {
    my $self = shift;
    $self->debug_message('Attempting to shortcut...');

    # Get indel aligner result
    my $indel_realigner_result = $self->_get_indel_realigner_result;
    if ( not $indel_realigner_result ) {
        $self->debug_message('Failed to find indel realigner result, cannot shortcut!');
        return;
    }

    # Get base recalibrator result
    my $base_recalibrator_result = $self->_get_base_recalibrator_result;
    if ( not $base_recalibrator_result ) {
        $self->debug_message('Failed to find base recalibrator result, cannot shortcut!');
        return;
    }

    $self->debug_message('Shortcut...OK');
    return $base_recalibrator_result;
}

sub execute {
    my $self = shift;

    # Try to shortcut
    my $shortcut = $self->shortcut;
    return 1 if $shortcut;

    # [Get or] Create indel aligner result
    my $indel_realigner_result = $self->_get_or_create_indel_realigner_result;
    if ( not $indel_realigner_result ) {
        $self->error_message('Failed to create indel realigner result!');
        return;
    }

    # [Get or] Create base recalibrator result
    my $base_recalibrator_result = $self->_get_or_create_base_recalibrator_result;
    if ( not $base_recalibrator_result ) {
        $self->error_message('Failed to create base recalibrator result!');
        return;
    }

    $self->debug_message('Execute...OK');
    return $base_recalibrator_result;
}

sub _load_result {
    my ($self, $result_name, $retrieval_method) = @_;

    # Check accessor
    my $result_method = $result_name;
    $result_method =~ s/\s/_/g;
    $result_method .= '_result';
    return $self->$result_method if $self->$result_method;

    $self->debug_message("Looking for $result_name result...");
    my $result_class = 'Genome::InstrumentData::Gatk::'.Genome::Utility::Text::string_to_camel_case($result_name).'Result';
    $self->debug_message('Class: '.$result_class);
    $self->debug_message('Method: '.$retrieval_method);

    my $params_method = '_params_for_'.$result_method;
    my %params = $self->$params_method;
    $self->debug_message("Params: \n".$self->_params_display_name(\%params));
    my $result = $result_class->$retrieval_method(%params);
    return if not $result; # let caller handle error/status
    $self->$result_method($result);

    $self->debug_message(ucfirst($result_name).': '.$result->__display_name__);
    $self->debug_message(ucfirst($result_name).' output directory: '.$result->output_dir);
    return $result;
}

sub _get_indel_realigner_result {
    return $_[0]->_load_result('indel realigner', 'get_with_lock');
}

sub _get_or_create_indel_realigner_result {
    return $_[0]->_load_result('indel realigner', 'get_or_create');
}
sub _get_base_recalibrator_result {
    return $_[0]->_load_result('base recalibrator bam', 'get_with_lock');
}

sub _get_or_create_base_recalibrator_result {
    return $_[0]->_load_result('base recalibrator bam', 'get_or_create');
}

sub _common_params_for_gatk_results {
    my $self = shift;

    my @known_sites = $self->known_sites;
    my %params = (
        version => $self->version,
        reference_build => $self->bam_source->reference_build,
        known_sites => \@known_sites,
    );

    return %params;
}

sub _params_for_indel_realigner_result {
    my $self = shift;

    my %indel_realigner_params = $self->_common_params_for_gatk_results;
    return if not %indel_realigner_params;

    $indel_realigner_params{bam_source} = $self->bam_source;

    return %indel_realigner_params;
}

sub _params_for_base_recalibrator_bam_result {
    my $self = shift;

    my %base_recalibrator_params = $self->_common_params_for_gatk_results;
    return if not %base_recalibrator_params;

    my $indel_realigner_result = $self->indel_realigner_result;
    if ( not $indel_realigner_result ) {
        $self->error_message('No indel aligner result set!');
        return;
    }
    $base_recalibrator_params{bam_source} = $indel_realigner_result;

    return %base_recalibrator_params;
}

sub _params_display_name {
    my ($self, $params) = @_;

    my $display_name;
    for my $key ( keys %$params ) {
        $display_name .= $key.': ';
        if ( my $ref = ref $params->{$key} ) {
            my @params = ( $ref eq 'ARRAY' ) ? @{$params->{$key}} : $params->{$key};
            for my $param ( @params ) {
                $display_name .= $param->id."\n";
            }
        }
        else {
            $display_name .=  $params->{$key}."\n";
        }
    }

    return $display_name;
}

1;

