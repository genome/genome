package Genome::Model::Command::Services::BuildQueuedModels;

use strict;
use warnings;

use Genome;
use Genome::Model::Build::Command::Start;

use POSIX qw(ceil);

class Genome::Model::Command::Services::BuildQueuedModels {
    is => 'Genome::Command::Base',
    doc => "Build queued models.",
    has_optional => [
        max_scheduled_builds => {
            is => 'Integer',
            default => 50,
        },
        channels => {
            is => 'Integer',
            default => 1,
            doc => 'number of "channels" to parallelize models by',
        },
        channel => {
            is => 'Integer',
            default => 0,
            doc => 'zero-based channel to use',
        },
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
            is_output => 1,
        },
        _builds_started => {
            is => 'Integer',
            default => 0,
        },
        _total_count => {
            is => 'Integer',
            default => 0,
        },
        _create_params => {
            is => 'Hash',
            default => {},
        },
        _start_params => {
            is => 'Hash',
            default => {},
        },
    ],
};

sub help_synopsis {
    return <<EOS;
genome model services build-queued-models
EOS
}

sub help_detail {
    return <<EOS;
Builds queued models.
EOS
}

# We had to debug a failure on the cron server where a downstream pipe crashed
# which meant that we got no mail due to having no SIGPIPE handler.
$SIG{PIPE} = sub {
    die "Cannot write to pipe (downstream pipe closed)";
};

sub execute {
    my $self = shift;

    unless ($self->channel < $self->channels) {
        die $self->error_message('--channel must be less than --channels');
    }

    my $lock_resource = $ENV{GENOME_LOCK_DIR} . '/genome_model_services_builed_queued_models_' . $self->channel . '_' . $self->channels;

    my $lock = Genome::Sys->lock_resource(resource_lock => $lock_resource, max_try => 1);
    unless ($lock) {
        $self->error_message("Could not lock, another instance of BQM must be running.");
        return;
    }

    my $context = UR::Context->current;
    $context->add_observer(
        aspect => 'commit',
        callback => sub{ Genome::Sys->unlock_resource(resource_lock => $lock) },
    );
    $context->add_observer(
        # this observer is so we can filter emails that were able to commit
        aspect => 'commit',
        callback => sub{ print "\nCOMMITTED\n" },
    );

    my $max_builds_to_start = $self->num_builds_to_start;
    unless ($max_builds_to_start) {
        $self->status_message("There are already " . $self->max_scheduled_builds . " builds scheduled.");
        return 1;
    }
    $self->status_message("Will try to start up to $max_builds_to_start builds.");

    my @iterator_params = (
        # prioritize genotype microarray over other builds because their
        # runtime is so short and are frequently prerequisite for other builds
        {build_requested => '1', type_name => 'genotype microarray', -order_by => 'subject_id'}, 
        {build_requested => '1', -order_by => 'subject_id'},
    );

    ITERATOR:
    while (my $iterator_params = shift @iterator_params) {
        my $models = Genome::Model->create_iterator(%{$iterator_params});

        MODEL:
        while (my $model = $models->next) {
            next MODEL unless ($model->id % $self->channels == $self->channel);

            if ($self->_builds_started >= $max_builds_to_start){
                $self->status_message("Already started max builds (" . $self->_builds_started . "), quitting...");
                last ITERATOR; 
            }

            $self->_total_command_count($self->_total_command_count + 1);
            Genome::Model::Build::Command::Start::create_and_start_build($self, $model);
        }
    }

    $self->display_command_summary_report();

    $self->status_message("Tried to start up to $max_builds_to_start builds.");

    if(not $self->_total_command_count) {
        $self->status_message("No queued models found to start.");
    }

    return !scalar(keys %{$self->_command_errors});
}


sub num_builds_to_start {
    my $self = shift;
    
    my $scheduled_builds = Genome::Model::Build->create_iterator(
        run_by => Genome::Sys->username,
        status => 'Scheduled',
    );
    
    my $scheduled_build_count = 0;
    while ($scheduled_builds->next && ++$scheduled_build_count <= $self->max_scheduled_builds) { 1; }
    
    my $max_per_channel = int($self->max_scheduled_builds / $self->channels);
    if ($scheduled_build_count >= $self->max_scheduled_builds) {
        return 0;
    }
    elsif (($scheduled_build_count + $max_per_channel) > $self->max_scheduled_builds) {
        return ceil($max_per_channel / $self->channels);
    }
    else {
        return $max_per_channel;
    }
}


1;
