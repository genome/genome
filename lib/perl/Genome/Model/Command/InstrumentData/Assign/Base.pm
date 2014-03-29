package Genome::Model::Command::InstrumentData::Assign::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::InstrumentData::Assign::Base {
    is => 'Command::V2',
    is_abstract => 1,
    has_input => [
        model => {
            is => 'Genome::Model',
            doc => 'model to which to assign instrument data',
        },
    ],
    has_optional_input => [
        force => {
            is => 'Boolean',
            default_value => 0,
            doc => 'unconditionally assign instrument data',
        },
    ],
    has_transient_optional => [
        _assigned_instrument_data => {
            is => 'HASH',
            default_value => {},
            doc => 'used internally by the command to hold already-assigned IDs for the model',
        },
    ],
    has_optional_output => [
        newly_added_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    my @model_instrument_data = $self->model->instrument_data;
    my $assigned_instrument_data = $self->_assigned_instrument_data;
    for my $instrument_data ( @model_instrument_data ) {
        $assigned_instrument_data->{ $instrument_data->id }++;
    }

    my @instrument_data = $self->_resolve_instrument_data;

    unless(@instrument_data) {
        die $self->error_message("No instrument data to assign.");
    }

    my $count = 0;
    for my $instrument_data (@instrument_data) {
        my $filter = $self->_resolve_filter_for_instrument_data($instrument_data);
        $count += $self->_assign_instrument_data_if_acceptable($instrument_data, $filter);
    }

    $self->status_message('Assigned %s instrument data to model.', $count);
    return 1;
}

sub _resolve_instrument_data {
    my $self = shift;
    my $class = ref $self || $self;

    die $self->error_message('The command %s must implement _resolve_instrument_data.', $class);
}

sub _resolve_filter_for_instrument_data {
    my $self = shift;
    my $instrument_data = shift;

    return;
}

sub _assign_instrument_data {
    my ($self, $instrument_data, $filter) = @_;

    # Check if already assigned
    my $model = $self->model;
    my $assigned_instrument_data = $self->_assigned_instrument_data;
    if ( exists $assigned_instrument_data->{ $instrument_data->id } ) {
        $self->status_message(
            'Instrument data (%s) already assigned to model (%s). Skipping.',
            $instrument_data->id,
            $model->__display_name__,
        );
        return 0;
    }
    $assigned_instrument_data->{ $instrument_data->id }++;

    my $add = $model->add_instrument_data(
        value => $instrument_data,
        filter_desc => $filter,
    );
    if ( not $add ) {
        $self->error_message(
            'Failed to add instrument data (%s) to model (%s).',
            $instrument_data->id,
            $model->__display_name__,
        );
        return;
    }

    $self->status_message(
        'Instrument data (%s) assigned to model (%s)%s.',
        $instrument_data->id,
        $model->__display_name__,
        ( $filter ? (' with filter ' . $filter) : '')
    );

    $self->add_newly_added_data($instrument_data);

    return 1;
}

sub _assign_instrument_data_if_acceptable {
    my $self = shift;
    my $instrument_data = shift;
    my $filter = shift;

    my @issues = @_; #if this method is extended by a subclass, pass other issues here


    if($instrument_data->ignored) {
        push @issues,
            'Ignored flag is set.';
    }

    unless($self->_instrument_data_matches_model_subject($instrument_data)) {
        push @issues,
            sprintf('Model subject (%s) does not match.', $self->model->subject->__display_name__);
    }

    unless($self->_instrument_data_matches_model_target_region_set($instrument_data)) {
        push @issues, 'Model target region set does not match.';
    }

    if(@issues and not $self->force) {
        $self->warning_message('Skipping instrument data %s because: %s',
            $instrument_data->__display_name__,
            join(" ", @issues)
        );
        return;
    }

    if(@issues and $self->force) {
        $self->status_message('Forcing assignment of instrument data %s despite: %s',
            $instrument_data->__display_name__,
            join(" ", @issues),
        );
    }

    return $self->_assign_instrument_data($instrument_data, $filter);
}

sub _instrument_data_matches_model_subject {
    my ($self, $instrument_data) = @_;

    my $subject = $self->model->subject;
    if($subject->isa('Genome::Taxon')) {
        return $instrument_data->taxon eq $subject;
    } elsif ($subject->isa('Genome::SampleSource')) {
        return $instrument_data->sample->source eq $subject;
    } elsif ($subject->isa('Genome::Sample')) {
        return $instrument_data->sample eq $subject;
    } else {
        die $self->error_message('Unexpected subject type for model: ' . $subject->class);
    }
}

sub _instrument_data_matches_model_target_region_set {
    my ($self, $instrument_data) = @_;

    my $model = $self->model;
    my @inputs = Genome::Model::Input->get(model_id => $model->id, name => 'target_region_set_name');
    my %model_capture_targets = map { $_->value_id() => 1 } @inputs;

    my @set_inputs = Genome::Model::Input->get(model_id => $model->id, name => 'target_region_set');
    for my $si (@set_inputs) { $model_capture_targets{$si->value->name} = 1; }

    my $id_capture_target = eval { $instrument_data->target_region_set_name };

    #if model has no TRS but instrument data does
    if(not %model_capture_targets and defined $id_capture_target) {
        return;
    }

    #if model has TRS but not what instrument data does
    if(%model_capture_targets and (not defined $id_capture_target or not exists $model_capture_targets{$id_capture_target})) {
        return;
    }

    return 1;
}

1;
