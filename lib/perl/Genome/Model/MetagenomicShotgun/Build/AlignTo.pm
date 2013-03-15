package Genome::Model::MetagenomicShotgun::Build::AlignTo;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicShotgun::Build::AlignTo {
    is => 'Command::V2',
    has_input => [
        input_build => {
            is => 'Genome::Model::Build::MetagenomicShotgun',
            is_many => 1,
            doc => 'The meta shot build to work with.',
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            doc => 'The instrument data to work with.',
        },
        sub_model_label => { 
            is => 'Text',
            valid_values => [ Genome::Model::MetagenomicShotgun->sub_model_labels ],
        },
    ],
    has_output => [
        build => {
            is => 'Genome::Model::Build::MetagenomicShotgun',
            calculate_from => ['input_build'],
            calculate => sub{ return $_[0]; },
        },
    ]
};

sub execute {
    my $self = shift;

    my $sub_model = $self->build->model->sub_model_for_label($self->sub_model_label);
    my @instrument_data = $self->instrument_data;

    my $assign_instrument_data_ok = $self->_assign_instrument_data($sub_model, @instrument_data);
    return if not $assign_instrument_data_ok;

    my $sub_build = $self->_build_if_necessary($sub_model);
    return if not $sub_build;

    my $sub_build_ok = $self->_wait_for_build($sub_build);
    return if not $sub_build_ok;

    return 1;
}

sub _assign_instrument_data  {
    my ($self, $model, @instrument_data) = @_;
    $self->status_message('Ensure correct assigned to sub model...');
    $self->status_message('Model: '.$model->__display_name__);
    $self->status_message('Instrument data: '.join(' ', map { $_->id } @instrument_data));

    # Ensure correct inst data on model
    my %assigned_instrument_data = map { $_->id => $_ } $model->instrument_data;
    my @instrument_data_ids_to_assign;
    for my $instrument_data ( @instrument_data ) {
        if ( $assigned_instrument_data{$instrument_data->id} ) {
            delete $assigned_instrument_data{$instrument_data->id};
        }
        else {
            push @instrument_data_ids_to_assign, $instrument_data->id;
        }
    }

    # Unassign incorrect inst data
    my $model_id = $model->id;
    my @instrument_data_ids_to_unassign = keys %assigned_instrument_data;
    if ( @instrument_data_ids_to_unassign ) {
        $self->status_message("Unassign incorrect instrument data ".join(' ', @instrument_data_ids_to_unassign)." from model ".$model->__display_name__);
        my $instrument_data_expression = 'id'.( @instrument_data_ids_to_unassign > 1 ) 
        ? '='.$instrument_data_ids_to_unassign[0]
        : ':'.join(@instrument_data_ids_to_unassign);
        my $cmd = "genome model instrument-data unassign $model_id --instrument-data $instrument_data_expression --force";
        my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
        if ( not $rv ) {
            $self->error_message($@) if $@;
            $self->error_message('Failed to unassign instrument data from model!');
            return;
        }
        $self->status_message('Unassign incorrect instrument data...done');
    }

    # Add correct inst data
    if ( @instrument_data_ids_to_assign ) {
        $self->status_message("Assign correct instrument data ".join(' ', @instrument_data_ids_to_assign)." to model ".$model->__display_name__);
        my $instrument_data_expression = 'id'.( @instrument_data_ids_to_assign > 1 ) 
        ? '='.$instrument_data_ids_to_assign[0]
        : ':'.join(@instrument_data_ids_to_assign);
        my $cmd = "genome model instrument-data unassign $model_id --instrument-data $instrument_data_expression --force";
        my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
        if ( not $rv ) {
            $self->error_message($@) if $@;
            $self->error_message('Failed to assign instrument data to model!');
            return;
        }
        $self->status_message('Assign correct instrument data...done');
    }

    # Reload the model and inputs
    UR::Context->reload('Genome::Model', id => $model->id);
    UR::Context->reload('Genome::Model::Input', model_id => $model->id);

    # Verify correct instrument data was added
    %assigned_instrument_data = map { $_->id => $_ } $model->instrument_data;
    for my $instrument_data ( @instrument_data ) {
        next if delete $assigned_instrument_data{$instrument_data->id};
        $self->error_message('Failed to assign correct instrument data!');
        return;
    }

    $self->status_message('Ensure correct assigned to sub model...OK');
    return 1;
}

sub _build_if_necessary {
    my ($self, @models) = @_;

    my (@succeeded_builds, @watched_builds);
    for my $model ( @models ) {
        $self->status_message('Model: '. $model->__display_name__);
        $self->status_message('Search for succeeded build');
        my $succeeded_build = $model->last_succeeded_build;
        if ( $succeeded_build and $self->_verify_model_and_build_instrument_data_match($model, $succeeded_build) ) {
            $self->status_message('Found succeeded build: '.$succeeded_build->__display_name__);
            push @succeeded_builds, $succeeded_build;
            next;
        }
        $self->status_message('No succeeded build');
        $self->status_message('Search for scheduled or running build');
        my $watched_build = $self->_find_scheduled_or_running_build_for_model($model);
        if ( not $watched_build ) {
            $self->status_message('No scheduled or running build');
            $self->status_message('Start build');
            $watched_build = $self->_start_build_for_model($model);
            return if not $watched_build;
        }
        $self->status_message('Watching build: '.$watched_build->__display_name__);
        push @watched_builds, $watched_build;
    }

    my @builds = (@succeeded_builds, @watched_builds);
    if ( not @builds ) {
        $self->error_message('Failed to find or start any builds');
        return;
    }

    if ( @models != @builds ) {
        $self->error_message('Failed to find or start a build for each model');
        return;
    }

    return ( @builds > 1 ? @builds : $builds[0] );
}

sub _verify_model_and_build_instrument_data_match {
    my ($self, $model, $build) = @_;

    Carp::confess('No model to verify instrument data') if not $model;
    Carp::confess('No build to verify instrument data') if not $build;

    my @build_instrument_data = sort {$a->id <=> $b->id} $build->instrument_data;
    my @model_instrument_data = sort {$a->id <=> $b->id} $model->instrument_data;

    $self->status_message('Model: '.$model->__display_name__);
    $self->status_message('Model instrument data: '.join(' ', map { $_->id } @model_instrument_data));
    $self->status_message('Build: '.$build->__display_name__);
    $self->status_message('Build instrument data: '.join(' ', map { $_->id } @build_instrument_data));

    if ( @build_instrument_data != @model_instrument_data ) {
        $self->status_message('Model and build instrument data count does not match');
        return;
    }

    for ( my $i = 0; $i < @model_instrument_data; $i++ ) {
        my $build_instrument_data = $build_instrument_data[$i];
        my $model_instrument_data = $model_instrument_data[$i];

        if ($build_instrument_data->id ne $model_instrument_data->id) {
            $self->status_message("Missing instrument data.");
            return;
        }
    }

    return 1;
}

sub _find_scheduled_or_running_build_for_model {
    my ($self, $model) = @_;

    Carp::confess('No model sent to find running or scheduled build') if not $model;

    $self->status_message('Looking for running or scheduled build for model: '.$model->__display_name__);

    UR::Context->reload('Genome::Model', id => $model->id);
    UR::Context->reload('Genome::Model::Build', model_id => $model->id);
    UR::Context->reload('Genome::Model::Event', model_id => $model->id);

    my $build = $model->latest_build;
    return if not $build;

    $self->status_message( sprintf('Build: %s %s', $build->id, $build->status) );
    if ( grep { $build->status eq $_ } (qw/ Scheduled Running /) ) {
        return $build;
    }

    return;
}

sub _start_build_for_model {
    my ($self, $model) = @_;

    Carp::confess('no model sent to start build') if not $model;

    my $cmd = 'genome model build start '.$model->id.' --job-dispatch apipe --server-dispatch workflow'; # these are defaults
    $self->status_message('cmd: '.$cmd);

    UR::Context->commit();
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        die $self->error_message('failed to execute build start command');
    }

    my $build = $self->_find_scheduled_or_running_build_for_model($model);
    if ( not $build ) {
        die $self->error_message('executed build start command, but cannot find build.');
    }

    return $build;
}

sub _wait_for_build {
    my ($self, $build) = @_;

    if ( not $build ) {
        $self->status_message("No build to wait!");
        return;
    }
    $self->status_message('Watching build: '.$build->__display_name__);

    my $last_status = '';
    my $time = 0;
    my $inc = 30;
    while (1) {
        UR::Context->current->reload($build->the_master_event);
        my $status = $build->status;
        if ($status and !($status eq 'Running' or $status eq 'Scheduled')){
            return 1;
        }

        if ($last_status ne $status or !($time % 300)){
            $self->status_message("Waiting for build(~$time sec) ".$build->id.", status: $status");
        }
        sleep $inc;
        $time += $inc;
        $last_status = $status;
    }

    my $status = $build->status;
    if ( $status eq 'Succeeded' ) {
        $self->status_message($status.'! '.$build->__display_name__);
        return 1;
    }
    else {
        $self->error_message($status.'! '.$build->__display_name__);
        return;
    }
}

1;

