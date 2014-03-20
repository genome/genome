package Genome::InstrumentData::BamUtil::ClipOverlapResult;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::BamUtil::ClipOverlapResult {
    is => 'Genome::InstrumentData::AlignedBamResult',
    has_input => [
        bam_source => { # PROVIDES bam_path SHOULD be in aligned bam result, but would be incompatible with AR::Merged
            is => 'Genome::InstrumentData::AlignedBamResult',
        },
    ],
    has_param => [
        version => {
            is => 'Text',
            doc => 'Version of BamUtil to use.',
            valid_values => [qw/ /], #TODO
        },
    ],
    has_constant => [
        output_bam_path => {
            calculate_from => ['temp_staging_directory'],
            calculate => q| File::Spec->join($temp_staging_directory, 'clipped.bam'); |,
        },
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
    $class =~ s/^Genome::InstrumentData::BamUtil:://;
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

    my $prepare_staging_directory = $self->_prepare_staging_directory;

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

    unless ($self->_run_clip_overlap) {
        $self->delete;
        return;
    }

    my $allocation = $self->disk_allocations;
    eval { $allocation->reallocate }; # FIXME pasted from GATk... I don't think this is needed? Should happen automatically

    return $self;
}

sub _run_clip_overlap {
    my $self = shift;

    my %clip_overlap_params = (
        version => $self->version,
        input_bam => $self->input_bam_path,
        output_bam => $self->output_bam_path,
    );
    $self->status_message('Params: '.Data::Dumper::Dumper(\%clip_overlap_params));

    my $clip_overlap = Genome::Model::Tools::BamUtil::ClipOverlap->create(%clip_overlap_params);
    if ( not $clip_overlap ) {
        $self->error_message('Failed to create clip overlap command!');
        return;
    }

    if ( not eval{ $clip_overlap->execute; } ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to execute clip overlap command!');
        return;
    }

    $self->_validate_results;
    return 1;
}

sub _validate_results {
    my $self = shift;

    if ( not -s $self->output_bam_path) {
        $self->error_message('Ran clip overlap command, but failed to make a output bam with size!');
        return;
    }

    $self->status_message('Output bam: ' . $self->output_bam_path);

    return 1;
}

1;

