package Genome::Model::MetagenomicShotgun::Build::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicShotgun {
    is => 'Command::V2',
    has => [
        build => {
            is => 'Genome::Model::Build::MetagenomicShotgun',
            doc => 'The meta shot build to work with.',
        },
    ],
};

sub _start_build  {
    my ($self, $model, @instrument_data) = @_;

    my %existing;
    for my $inst_data($model->instrument_data){
        $existing{$inst_data->id} = $inst_data;
    }
    my @to_add;
    for my $inst_data(@instrument_data){
        if ($existing{$inst_data->id}) {
            delete $existing{$inst_data->id};
        }
        else {
            push @to_add, $inst_data;
        }
    }
    for my $inst_data (values %existing){
        $self->status_message("Removing Instrument Data " . $inst_data->id . " from model " . $model->__display_name__);
        $model->remove_instrument_data($inst_data)
    }
    for my $inst_data(@to_add){
        $self->status_message("Adding Instrument Data " . $inst_data->id . " to model " . $model->__display_name__);
        $model->add_instrument_data($inst_data);
    }
    if (@instrument_data == 0){
        $self->status_message("No instrument data for model ".$model->__display_name__.", skipping build");
        return;
    }

    my $build = $self->_build_if_necessary($model);
    return $build;
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

sub _extract_data {
    my ($self, $from_build, $extraction_type) = @_;


    my $extract_from_alignment = Genome::Model::MetagenomicShotgun::Build::ExtractFromAlignment->create(
        build => $self,
        sub_build => $from_build,
        type => $extraction_type,
    );
    if ( not $extract_from_alignment ) {
        $self->error_message("Failed to create extract from alignment! Tried to extract $extraction_type from ".$from_build->__display_name__);
        return;
    }

    my $execute_ok = $extract_from_alignment->execute;
    if ( not $execute_ok ) {
        $self->error_message("Failed to execute extract from alignment! Tried to extract $extraction_type from ".$from_build->__display_name__);
        return;
    }

    return 1;
}

sub _link_sub_build_alignments_to_build {
    my ($self, %params) = @_;

    my $build = delete $params{build};
    Carp::confess('No build given to link alignments!') if not $build;
    my $sub_build = delete $params{sub_build};
    Carp::confess('No sub-build given to link alignments!') if not $sub_build;
    my $sub_model_name = delete $params{sub_model_name};
    Carp::confess('No sub model name given to link alignments!') if not $sub_model_name;
    Carp::confess('Unknown params given to _link_sub_build_alignments_to_build! '.Data::Dumper::Dumper(\%params)) if %params;

    my $dir = $build->data_directory;
    my $sub_dir = $dir.'/'.$sub_model_name;
    my $create_ok = eval{ Genome::Sys->create_directory($sub_dir); };
    if ( not $create_ok ) {
        $self->error_message($@) if $@;
        $self->error_message("Failed to create $sub_model_name sub dir! ".$sub_dir);
        return;
    }

    for my $instrument_data ( $sub_build->instrument_data ) {
        my @alignments = $sub_build->alignment_results_for_instrument_data($instrument_data); # This should only be one.
        for my $alignment ( @alignments ) {
            my $target = $alignment->output_dir;
            my $link = $sub_dir.'/'.$instrument_data->id.'-'.$alignment->id;
            unlink $link;
            my $link_ok = eval{ Genome::Sys->create_symlink($target, $link); };
            if ( not $link_ok ) {
                $self->error_message($@) if $@;
                $self->error_message("Failed to create symlink! From $link to $target");
                return;
            }
        }
    }

    return 1;
}

1;

