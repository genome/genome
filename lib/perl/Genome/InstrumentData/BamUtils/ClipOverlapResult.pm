package Genome::InstrumentData::BamUtils::ClipOverlapResult;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::BamUtils::ClipOverlapResult {
    is => 'Genome::InstrumentData::AlignedBamResult',
    has_input => [
        bam_source => { # PROVIDES bam_path SHOULD be in aligned bam result, but would be incompatible with AR::Merged
            is => 'Genome::InstrumentData::AlignedBamResult',
        },
    ],
    has_param => [
        version => {
            is => 'Text',
            doc => 'Version of BamUtils to use.',
            valid_values => [qw/ /], #TODO
        },
    ],
    has_constant => [
        _tmpdir => {  calculate => q| return File::Temp::tempdir(CLEANUP => 1); |, },
        # from inputs
        input_bam_path => { 
            via => 'bam_source',
            to => 'bam_path', 
        },
    ],
};

sub instrument_data { # SHOULD be in AlignedBamResult, but would be incompatible with AR::Merged
    my $original_source = shift;

    do {
        $original_source = $original_source->bam_source
    } until not $original_source->can('bam_source');

    if ( $original_source->isa('Genome::InstrumentData') ) {
        return $original_source;
    }
    elsif ( $original_source->can('instrument_data') ) {
        return $original_source->instrument_data;
    }
    return;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $class = $self->class;
    $class =~ s/^Genome::InstrumentData::BamUtils:://;
    $class =~ s/Result$//;
    return sprintf(
        "model_data/gatk/%s-%s-%s-%s-%s", 
        Genome::Utility::Text::camel_case_to_string($class, '_'), 
        Sys::Hostname::hostname(),
        $ENV{USER}, $$, $self->id,
    );
}

sub resolve_allocation_kilobytes_requested {
    my $self = shift;
    my $kb_requested = -s $self->input_bam_path;
    return int($kb_requested / 1024 * 1.5);
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $prepare_output_directory = eval{ $self->_prepare_output_directory; };
    if ( not $prepare_output_directory ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to prepare output directory!') if $@;
        $self->delete;
        return;
    }

    my $bam_source = $self->bam_source;
    $self->status_message('Bam source: '.$bam_source->id);
    $bam_source->add_user(user => $self, label => 'bam source');

    $self->status_message('Reference: '.$self->reference_build->id);

    return $self;
}

1;

