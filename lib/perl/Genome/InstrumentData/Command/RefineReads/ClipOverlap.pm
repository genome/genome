package Genome::InstrumentData::Command::RefineReads::ClipOverlap;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::BamUtil::ClipOverlapResult;

class Genome::InstrumentData::Command::RefineReads::ClipOverlap {
    is => 'Command::V2',
    has => [
        version => { is => 'Text', },
        bam_source => { is => 'Genome::InstrumentData::AlignedBamResult', },
    ],
    has_optional => [
        params => { is => 'Text', },
    ],
    has_many_optional => [
        known_sites => { is => 'Genome::Model::Build::ImportedVariationList', }, # FIXME not needed... add stuff to a base class and make this optional (but required for gatk)
    ],
    has_optional_transient => [
        result => { is => 'Genome::InstrumentData::BamUtil::ClipOverlapResult', },
    ],
};

sub shortcut {
    my $self = shift;
    $self->debug_message('Attempting to shortcut...');

    my $result_class = 'Genome::InstrumentData::BamUtil::ClipOverlapResult';
    my %params = $self->_params_for_result;
    my $result = $result_class->get_with_lock(%params);

    if ( not $result ) {
        $self->debug_message('Failed to find clip overlap result, cannot shortcut!');
        return;
    }

    $self->debug_message('Shortcut...OK');
    return $result;
}

sub execute {
    my $self = shift;

    # Try to shortcut
    my $shortcut = $self->shortcut;
    return $shortcut if $shortcut;

    # [Get or] Create clip overlap result
    my $clip_overlap_result = $self->_get_or_create_clip_overlap_result;
    if ( not $clip_overlap_result ) {
        $self->error_message('Failed to create clip overlap result!');
        return;
    }

    $self->debug_message('Execute...OK');
    return $clip_overlap_result;
}

sub _get_or_create_clip_overlap_result {
    my $self = shift;

    # Check accessor
    return $self->result if $self->result;

    my $result_class = 'Genome::InstrumentData::BamUtil::ClipOverlapResult';
    my %params = $self->_params_for_result;
    my $result = $result_class->get_or_create(%params);
    return if not $result; # let caller handle error/status

    $self->result($result);

    return $result;
}

sub _params_for_result {
    my $self = shift;

    my %params = (
        version => $self->version,
        bam_source => $self->bam_source,
        reference_build => $self->bam_source->reference_build,
        # FIXME more needed here probably, params?
    );

    return %params;
}

1;

