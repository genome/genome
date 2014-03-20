package Genome::InstrumentData::Command::RefineReads::ClipOverlap;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::BamUtils::ClipOverlapResult;
require Genome::Utility::Text;

class Genome::InstrumentData::Command::RefineReads::ClipOverlap {
    is => 'Command::V2',
    has => [
        version => { is => 'Text', },
        bam_source => { is => 'Genome::InstrumentData::AlignedBamResult', },
    ],
    has_optional => [
        params => { is => 'Text', },
    ],
    has_optional_transient => [
        clip_overlap_result => { is => 'Genome::InstrumentData::BamUtils::ClipOverlapResult', },
    ],
};

sub shortcut {
    my $self = shift;
    $self->debug_message('Attempting to shortcut...');

    my $clip_overlap_result = $self->_get_clip_overlap_result;
    if ( not $clip_overlap_result ) {
        $self->debug_message('Failed to find clip overlap result, cannot shortcut!');
        return;
    }

    $self->debug_message('Shortcut...OK');
    return $clip_overlap_result;
}

sub execute {
    my $self = shift;

    # Try to shortcut
    my $shortcut = $self->shortcut;
    return 1 if $shortcut;

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

};

sub _load_result {
    my ($self, $retrieval_method) = @_;

    # Check accessor
    return $self->clip_overlap_result if $self->clip_overlap_result;

    my $result_class = 'Genome::InstrumentData::BamUtils::ClipOverlapResult';
    $self->debug_message('Class: '.$result_class);
    $self->debug_message('Method: '.$retrieval_method);

    my %params = $self->_params_for_result;
    $self->debug_message("Params: \n".$self->_params_display_name(\%params));
    my $result = $result_class->$retrieval_method(%params);
    return if not $result; # let caller handle error/status
    $self->clip_overlap_result($result);

    $self->debug_message('Clip overlap: '.$result->__display_name__);
    $self->debug_message('Clip overlap output directory: '.$result->output_dir);
    return $result;
}

sub _params_for_result {
    my $self = shift;

    my %params = (
        version => $self->version,
        bam_source => $self->bam_source,
        # FIXME more needed here probably, params?
    );

    return %params;
}

#FIXME is this needed? Put in a Base class somewhere since this is pasted from gatk?
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

