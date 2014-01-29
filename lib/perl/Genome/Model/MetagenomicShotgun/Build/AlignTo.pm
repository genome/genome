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
    ],
    has_optional_transient => [
        sub_model => {
            is => 'Genome::Model',
            calculate_from => [qw/ build sub_model_label /],
            calculate => sub{
                my ($build, $sub_model_label) = @_;
                my $sub_model = $build->model->sub_model_for_label($sub_model_label);
                return $sub_model;
            },
        },
    ],
};

sub execute {
    my $self = shift;

    my $assign_instrument_data_ok = $self->_assign_instrument_data;
    return if not $assign_instrument_data_ok;

    my $sub_build = $self->_build_if_necessary;
    return if not $sub_build;

    my $sub_build_ok = $self->_wait_for_build($sub_build);
    return if not $sub_build_ok;

    return 1;
}

sub _assign_instrument_data  {
    my $self = shift;
    $self->debug_message('Ensure correct assigned to sub model...');

    my $model = $self->sub_model;
    if ( not $model ) {
        $self->error_message('Failed to get sub model!');
        return;
    }
    $self->debug_message('Model: '.$model->__display_name__);

    my @instrument_data = $self->instrument_data;
    $self->debug_message('Instrument data: '.join(' ', map { $_->id } @instrument_data));

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
        $self->debug_message("Unassign incorrect instrument data ".join(' ', @instrument_data_ids_to_unassign)." from model ".$model->__display_name__);
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
        $self->debug_message('Unassign incorrect instrument data...done');
    }

    # Add correct inst data
    if ( @instrument_data_ids_to_assign ) {
        $self->debug_message("Assign correct instrument data ".join(' ', @instrument_data_ids_to_assign)." to model ".$model->__display_name__);
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
        $self->debug_message('Assign correct instrument data...done');
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
    if ( %assigned_instrument_data ) {
        $self->error_message('Found incorrect instrument data still assigned to model!');
        return;
    }

    $self->debug_message('Ensure correct assigned to sub model...OK');
    return 1;
}

sub _build_if_necessary {
    my ($self, $model) = @_;
    $self->debug_message('Build if necessary...');

    $self->debug_message('Model: '. $model->__display_name__);
    $self->debug_message('Search for succeeded build...');
    my $sub_build = $model->build_needed;
    if ( $sub_build ) {
        $self->debug_message('Found current build! '.$sub_build->__display_name__);
        return;
    }

    $self->debug_message('No succeeded build. Look for a build that is scheduled or running...');
    $sub_build = $self->_find_scheduled_or_running_build_for_model($model);
    if ( not $sub_build ) {
        $self->debug_message('None found. Creating a build...');
        my $cmd = 'genome model build start '.$model->id.' --job-dispatch apipe --server-dispatch workflow'; # these are defaults
        $self->debug_message('cmd: '.$cmd);
        my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
        if ( not $rv ) {
            $self->error_message($@);
            $self->error_message('Failed to execute build start command!');
            return;
        }
        $sub_build = $self->_find_scheduled_or_running_build_for_model($model);
        if ( not $sub_build ) {
            $self->error_message('Executed the build start command, but cannot the build!');
            return;
        }
    }
    else {
        $self->debug_message('Found scheduled/running build!');
    }
    $self->debug_message('Build id: '.$sub_build->id);
    $self->debug_message('Build data directory '.$sub_build->data_directory);
    $self->debug_message('Watching build...');
    my $time = 0;
    my $inc = 30;
    my $status = $sub_build->status;
    while ( $status eq 'Running' or $status eq 'Scheduled' ) {
        if ( $time % 150 ){ # Only report every 5 minutes
            $self->debug_message("Watching build. Time: $time Status: $status");
        }
        sleep $inc;
        $time += $inc;
        UR::Context->current->reload($sub_build->the_master_event);
        $status = $sub_build->status;
    }

    $self->debug_message('Build has stopped running, checking status...');
    $self->debug_message('Status: '.$status);
    return ( $status eq 'Succeeded' ? 1 : 0 );
}

sub _find_scheduled_or_running_build_for_model {
    my ($self, $model) = @_;

    Carp::confess('No model sent to find running or scheduled build') if not $model;

    $self->debug_message('Looking for running or scheduled build for model: '.$model->__display_name__);

    UR::Context->reload('Genome::Model', id => $model->id);
    UR::Context->reload('Genome::Model::Build', model_id => $model->id);
    UR::Context->reload('Genome::Model::Event', model_id => $model->id);

    my $build = $model->latest_build;
    return if not $build;

    $self->debug_message( sprintf('Build: %s %s', $build->id, $build->status) );
    if ( grep { $build->status eq $_ } (qw/ Scheduled Running /) ) {
        return $build;
    }

    return;
}

1;

